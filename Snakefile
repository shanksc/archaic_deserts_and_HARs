configfile: "config.yaml"

INTERMEDIATE = f"{config['output_dir']}/intermediate"
RESULTS = f"{config['output_dir']}/final_results"

rule all:
    input:
        expand("{results}/archaic_enrichment_results.txt", results=RESULTS),
        expand("{results}/har_enrichment_results.txt", results=RESULTS) if config.get("HAR_regions") else [],
rule parse_gencode_gtf:
    input:
        gtf=config["gencode_gtf"]
    output:
        bed=f"{INTERMEDIATE}/gencode_genes_all.bed"
    shell:
        "python3 scripts/parse_gencode_gtf.py {input.gtf} {output.bed}"

rule filter_protein_coding:
    input:
        bed=f"{INTERMEDIATE}/gencode_genes_all.bed"
    output:
        bed=f"{INTERMEDIATE}/gencode_genes.bed"
    shell:
        "awk '$7 == \"protein_coding\"' {input.bed} | sort -k1,1 -k2,2n > {output.bed}"

rule filter_gwas_by_pvalue:
    input:
        gwas=config["gwas_file"]
    output:
        filtered=f"{INTERMEDIATE}/gwas_filtered.txt"
    params:
        pval=config["p_threshold"],
        pcol=config["p_column"]
    run:
        import pandas as pd
        df = pd.read_csv(input.gwas, sep=r"\s+", na_values=['-nan', 'nan', 'NA'])
        filtered = df[df[params.pcol] < params.pval]
        filtered.to_csv(output.filtered, sep="\t", index=False)
        print(f"Filtered {len(filtered)} variants with p < {params.pval}")

rule convert_to_bed:
    input:
        filtered=f"{INTERMEDIATE}/gwas_filtered.txt"
    output:
        bed=f"{INTERMEDIATE}/gwas_variants.bed"
    run:
        import pandas as pd
        df = pd.read_csv(input.filtered, sep="\t")
        with open(output.bed, 'w') as f:
            for _, row in df.iterrows():
                chrom = str(row['chr']) if not str(row['chr']).startswith('chr') else str(row['chr']).replace('chr', '')
                pos = int(row['ps'])
                rsid = row['rs']
                allele1 = row['allele1']
                allele0 = row['allele0']
                beta = row['beta']
                pval = row[config['p_column']]
                f.write(f"chr{chrom}\t{pos}\t{pos+1}\t{rsid}\t{beta}\t.\t{allele0}\t{allele1}\t{pval}\n")

rule sort_bed:
    input:
        bed=f"{INTERMEDIATE}/gwas_variants.bed"
    output:
        sorted_bed=f"{INTERMEDIATE}/gwas_significant_variants.bed"
    shell:
        "sort -k1,1 -k2,2n {input.bed} > {output.sorted_bed}"

rule assign_nearest_gene:
    input:
        variants=f"{INTERMEDIATE}/gwas_significant_variants.bed",
        genes=f"{INTERMEDIATE}/gencode_genes.bed"
    output:
        annotated=f"{INTERMEDIATE}/gwas_variants_annotated.bed"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools closest -a {input.variants} -b {input.genes} -d > {output.annotated}"

rule extract_gene_list:
    input:
        annotated=f"{INTERMEDIATE}/gwas_variants_annotated.bed"
    output:
        genes=f"{INTERMEDIATE}/gwas_genes.txt"
    shell:
        "cut -f13 {input.annotated} | awk -F'.' '{{print $1}}' | sort | uniq > {output.genes}"

rule calculate_gene_density:
    input:
        genes=f"{INTERMEDIATE}/gencode_genes.bed"
    output:
        density=f"{INTERMEDIATE}/gene_density.bed"
    params:
        window=config["gene_density_window"]
    run:
        import pandas as pd
        genes_df = pd.read_csv(input.genes, sep="\t", header=None,
                               names=['chrom', 'start', 'end', 'gene_id', 'gene_name', 'strand', 'gene_type'])
        
        results = []
        for idx, gene in genes_df.iterrows():
            chrom = gene['chrom']
            start = int(gene['start'])
            end = int(gene['end'])
            center = (start + end) // 2
            window_start = max(0, center - params.window)
            window_end = center + params.window
            
            nearby = genes_df[
                (genes_df['chrom'] == chrom) &
                (genes_df['end'].astype(int) >= window_start) &
                (genes_df['start'].astype(int) <= window_end)
            ]
            gene_density = len(nearby)
            gene_length = end - start
            
            results.append([chrom, start, end, gene['gene_id'], gene['gene_name'], 
                          gene_density, gene_length])
        
        result_df = pd.DataFrame(results, columns=['chrom', 'start', 'end', 'gene_id', 
                                                   'gene_name', 'gene_density', 'gene_length'])
        result_df.to_csv(output.density, sep="\t", index=False, header=True)

rule overlap_with_archaic_deserts:
    input:
        density=f"{INTERMEDIATE}/gene_density.bed",
        deserts=config["archaic_deserts"]
    output:
        overlap=f"{INTERMEDIATE}/genes_archaic_overlap.bed"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.density} -b {input.deserts} -c > {output.overlap}
        """

rule filter_human_derived_variants:
    input:
        catalog=config["archaic_catalog"]
    output:
        bed=f"{INTERMEDIATE}/human_derived_variants.bed"
    run:
        import gzip
        import pandas as pd
        
        results = []
        with gzip.open(input.catalog, 'rt') as f:
            header = next(f).strip().split('\t')
            col_idx = {name: i for i, name in enumerate(header)}
            
            for line in f:
                vals = line.strip().split('\t')
                try:
                    daf = float(vals[col_idx['human_DAF']])
                    arch_freq = float(vals[col_idx['Arch_freq']])
                except:
                    continue
                
                canc = vals[col_idx['CAnc']]
                if canc in ['NA', 'N', '.', '*']:
                    continue
                
                # Human-DERIVED: high-frequency derived in humans (DAF>=0.9), ancestral in archaics
                # DAF is relative to derived allele
                altai_anc = vals[col_idx['Altai_allele']] == canc if vals[col_idx['Altai_allele']] != 'NA' else False
                vindija_anc = vals[col_idx['Vindija_allele']] == canc if vals[col_idx['Vindija_allele']] != 'NA' else False
                denisova_anc = vals[col_idx['Denisova_allele']] == canc if vals[col_idx['Denisova_allele']] != 'NA' else False
                
                # Human high-frequency derived (HHMC) = DAF>=0.9 AND archaics carry ancestral
                if daf >= 0.9 and sum([altai_anc, vindija_anc, denisova_anc]) >= 2:
                    pos = vals[col_idx['POS']]
                    chrom, coord = pos.split(':')
                    coord = int(coord)
                    results.append([f"chr{chrom}", coord, coord+1, pos, vals[col_idx['REF']], vals[col_idx['ALT']]])
        
        df = pd.DataFrame(results, columns=['chrom', 'start', 'end', 'id', 'ref', 'alt'])
        df = df.sort_values(['chrom', 'start'])
        df.to_csv(output.bed, sep="\t", index=False, header=False)

rule filter_archaic_derived_variants:
    input:
        catalog=config["archaic_catalog"]
    output:
        bed=f"{INTERMEDIATE}/archaic_derived_variants.bed"
    run:
        import gzip
        import pandas as pd
        
        results = []
        with gzip.open(input.catalog, 'rt') as f:
            header = next(f).strip().split('\t')
            col_idx = {name: i for i, name in enumerate(header)}
            
            for line in f:
                vals = line.strip().split('\t')
                try:
                    daf = float(vals[col_idx['human_DAF']])
                    arch_freq = float(vals[col_idx['Arch_freq']])
                except:
                    continue
                
                canc = vals[col_idx['CAnc']]
                if canc in ['NA', 'N', '.', '*']:
                    continue
                
                # Archaic-DERIVED: all 3 archaics carry derived (not ancestral), rare/absent in humans
                # DAF is derived allele frequency
                if arch_freq == 3.0 and daf <= 0.1:
                    pos = vals[col_idx['POS']]
                    chrom, coord = pos.split(':')
                    coord = int(coord)
                    results.append([f"chr{chrom}", coord, coord+1, pos, vals[col_idx['REF']], vals[col_idx['ALT']]])
        
        df = pd.DataFrame(results, columns=['chrom', 'start', 'end', 'id', 'ref', 'alt'])
        df = df.sort_values(['chrom', 'start'])
        df.to_csv(output.bed, sep="\t", index=False, header=False)

rule assign_derived_variants_to_genes:
    input:
        human_vars=f"{INTERMEDIATE}/human_derived_variants.bed",
        archaic_vars=f"{INTERMEDIATE}/archaic_derived_variants.bed",
        genes=f"{INTERMEDIATE}/gencode_genes.bed"
    output:
        human_annotated=f"{INTERMEDIATE}/human_derived_variants_annotated.bed",
        archaic_annotated=f"{INTERMEDIATE}/archaic_derived_variants_annotated.bed"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools closest -a {input.human_vars} -b {input.genes} -d > {output.human_annotated}
        bedtools closest -a {input.archaic_vars} -b {input.genes} -d > {output.archaic_annotated}
        """

rule overlap_genes_with_hars:
    input:
        genes=f"{INTERMEDIATE}/gencode_genes.bed",
        har=config["HAR_regions"]
    output:
        overlap=f"{INTERMEDIATE}/genes_har_overlap.bed"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        tail -n +2 {input.har} | bedtools intersect -a {input.genes} -b - -c > {output.overlap}
        """

rule test_har_enrichment:
    input:
        overlap=f"{INTERMEDIATE}/genes_har_overlap.bed",
        endurance_genes=f"{INTERMEDIATE}/gwas_genes.txt",
        human_vars=f"{INTERMEDIATE}/human_derived_variants.bed",
        har=config["HAR_regions"]
    output:
        results=f"{RESULTS}/har_enrichment_results.txt"
    shell:
        """
        python3 scripts/test_har_enrichment.py \
            {input.overlap} \
            {input.endurance_genes} \
            {input.human_vars} \
            {input.har} \
            {output.results}
        """

rule test_archaic_enrichment:
    input:
        overlap=f"{INTERMEDIATE}/genes_archaic_overlap.bed",
        endurance_genes=f"{INTERMEDIATE}/gwas_genes.txt"
    output:
        results=f"{RESULTS}/archaic_enrichment_results.txt"
    run:
        import pandas as pd
        import numpy as np
        from scipy import stats
        
        overlap_df = pd.read_csv(input.overlap, sep="\t")
        overlap_df.columns = ['chrom', 'start', 'end', 'gene_id', 'gene_name', 
                              'gene_density', 'gene_length', 'desert_overlap']
        
        endurance = set()
        with open(input.endurance_genes) as f:
            for line in f:
                endurance.add(line.strip())
        
        overlap_df['is_endurance'] = overlap_df['gene_id'].apply(
            lambda x: 1 if x.split('.')[0] in endurance else 0
        )
        overlap_df['in_desert'] = (overlap_df['desert_overlap'] > 0).astype(int)
        
        overlap_df['log_gene_length'] = np.log10(overlap_df['gene_length'] + 1)
        overlap_df['log_gene_density'] = np.log10(overlap_df['gene_density'] + 1)
        
        with open(output.results, 'w') as out:
            out.write("Archaic Desert Enrichment\n\n")
            out.write(f"Total protein-coding genes: {len(overlap_df)}\n")
            out.write(f"Endurance genes (GWAS): {overlap_df['is_endurance'].sum()}\n")
            out.write(f"Genes in archaic deserts: {overlap_df['in_desert'].sum()}\n\n")
            
            endurance_in_desert = overlap_df[overlap_df['is_endurance']==1]['in_desert'].sum()
            endurance_total = overlap_df['is_endurance'].sum()
            background_in_desert = overlap_df[overlap_df['is_endurance']==0]['in_desert'].sum()
            background_total = len(overlap_df) - endurance_total
            
            out.write("Raw Counts\n")
            out.write(f"Endurance genes in deserts: {endurance_in_desert}/{endurance_total} ({100*endurance_in_desert/endurance_total:.1f}%)\n")
            out.write(f"Background genes in deserts: {background_in_desert}/{background_total} ({100*background_in_desert/background_total:.1f}%)\n\n")
            
            contingency = [[endurance_in_desert, endurance_total - endurance_in_desert],
                          [background_in_desert, background_total - background_in_desert]]
            odds_ratio, p_fisher = stats.fisher_exact(contingency)
            out.write(f"Fisher's exact test p-value: {p_fisher:.3e}\n")
            out.write(f"Odds ratio: {odds_ratio:.3f}\n\n")
            
            import statsmodels.api as sm
                
            X = overlap_df[['is_endurance', 'log_gene_density', 'log_gene_length']]
            y = overlap_df['in_desert']
            X = sm.add_constant(X)
                
            model = sm.Logit(y, X).fit(disp=0)
            
            out.write("Logistic Regression\n")
            out.write(model.summary().as_text())
            out.write("\n\n")
            
            coef = model.params['is_endurance']
            pval = model.pvalues['is_endurance']
            or_adj = np.exp(coef)
            ci = np.exp(model.conf_int().loc['is_endurance'])
            
            out.write(f"Adjusted odds ratio: {or_adj:.3f} (95% CI: {ci[0]:.3f}-{ci[1]:.3f})\n")
            out.write(f"P-value for endurance enrichment: {pval:.3e}\n")
            