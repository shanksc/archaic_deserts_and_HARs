import pandas as pd
import numpy as np
from scipy import stats
import sys

def run_har_enrichment(overlap_file, endurance_file, human_vars_file, har_bed, output_file):
    """Test for HAR enrichment in endurance genes."""
    
    overlap = pd.read_csv(overlap_file, sep="\t", header=None)
    overlap.columns = ['chrom', 'start', 'end', 'gene_id', 'gene_name', 'strand', 'gene_type', 'har_count']
    
    endurance = set()
    with open(endurance_file) as f:
        for line in f:
            endurance.add(line.strip())
    
    overlap['gene_id_short'] = overlap['gene_id'].str.split('.').str[0]
    overlap['is_endurance'] = overlap['gene_id_short'].isin(endurance).astype(int)
    overlap['has_har'] = (overlap['har_count'] > 0).astype(int)
    
    # Count overlaps
    endurance_genes = overlap[overlap['is_endurance'] == 1]
    background_genes = overlap[overlap['is_endurance'] == 0]
    
    endurance_har = endurance_genes['has_har'].sum()
    endurance_total = len(endurance_genes)
    background_har = background_genes['has_har'].sum()
    background_total = len(background_genes)
    
    # Fisher's exact test
    contingency = [[endurance_har, endurance_total - endurance_har],
                   [background_har, background_total - background_har]]
    odds_ratio, p_fisher = stats.fisher_exact(contingency)
    
    # Write TSV with HAR overlap for endurance genes
    tsv_file = output_file.replace('.txt', '_genes.tsv')
    endurance_har_data = endurance_genes[['gene_name', 'gene_id_short', 'har_count', 'chrom', 'start', 'end']].copy()
    endurance_har_data = endurance_har_data.sort_values('har_count', ascending=False)
    endurance_har_data.to_csv(tsv_file, sep="\t", index=False)
    
    # Write results
    with open(output_file, 'w') as out:
        out.write("HAR Enrichment Analysis\n\n")
        
        out.write(f"Total genes: {len(overlap)}\n")
        out.write(f"Genes with HARs: {overlap['has_har'].sum()}\n")
        out.write(f"Endurance genes: {endurance_total}\n\n")
        
        out.write("HAR OVERLAP\n\n")
        
        out.write(f"Endurance genes with HARs: {endurance_har}/{endurance_total} ({100*endurance_har/endurance_total:.1f}%)\n")
        out.write(f"Background genes with HARs: {background_har}/{background_total} ({100*background_har/background_total:.1f}%)\n\n")
        
        out.write(f"Fisher's exact test:\n")
        out.write(f"  Odds ratio: {odds_ratio:.3f}\n")
        out.write(f"  P-value: {p_fisher:.3e}\n\n")
        
        out.write("TOP ENDURANCE GENES BY HAR COUNT\n\n")
        
        endurance_har_genes = endurance_genes[endurance_genes['har_count'] > 0].copy()
        if len(endurance_har_genes) > 0:
            endurance_har_genes = endurance_har_genes.sort_values('har_count', ascending=False)
            top = endurance_har_genes.head(20)[['gene_name', 'gene_id_short', 'har_count']]
            out.write(top.to_string(index=False))
            out.write(f"\n\n({len(endurance_har_genes)} total endurance genes with HARs)")
            out.write(f"\nFull gene list saved to: {tsv_file}\n")
        else:
            out.write("No endurance genes overlap HARs\n")

if __name__ == '__main__':
    run_har_enrichment(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])