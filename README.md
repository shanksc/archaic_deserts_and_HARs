## Dependencies 

```bash
mamba create -n archaic_analysis python=3.10
mamba activate archaic_analysis
mamba install -c bioconda -c conda-forge snakemake bedtools pandas scipy statsmodels scikit-learn
```

### Running pipeline 
```
snakemake -c1
```


## Configuration

Edit `config.yaml` to modify parameters 


## Results

Final results are in `results/final_results/`:

- `archaic_enrichment_results.txt` - Tests whether endurance genes are enriched in archaic deserts (regions depleted of Neanderthal ancestry). Includes Fisher's exact test and logistic regression controlling for gene density and length.

- `har_enrichment_results.txt` - Tests whether endurance genes overlap Human Accelerated Regions (HARs) more than expected. 

- `har_enrichment_results_genes.tsv` - Complete list of all endurance genes with HAR overlap counts (gene name, ENSG ID, HAR count, coordinates).


