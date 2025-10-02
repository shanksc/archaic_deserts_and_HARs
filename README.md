## Dependencies 

```bash
conda create -n archaic_analysis python=3.10
conda activate archaic_analysis
conda install -c bioconda -c conda-forge snakemake bedtools pandas scipy statsmodels scikit-learn
```

### Running pipeline 
```
snakemake --use-conda -c 1
```


## Configuration

Edit `config.yaml` to modify:
```yaml
p_threshold: 1.0e-5              # GWAS p-value threshold
gene_density_window: 5000000     # for gene density calculation
protein_coding_only: true         # Filter non-coding genes
```

## Results

Final results are in `results/final_results/`:

- `archaic_enrichment_results.txt` - Tests whether endurance genes are enriched in archaic deserts (regions depleted of Neanderthal ancestry). Includes Fisher's exact test and logistic regression controlling for gene density and length.

- `har_enrichment_results.txt` - Tests whether endurance genes overlap Human Accelerated Regions (HARs) more than expected. 

- `har_enrichment_results_genes.tsv` - Complete list of all endurance genes with HAR overlap counts (gene name, ENSG ID, HAR count, coordinates).


