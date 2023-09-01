# methyl_longRead_wf
a snakemake workflow to compute methylation activity of genes with long read sequencing and associate with gene expression

## how to install pipeline
```
$ git clone git@github.com:sahuno/methyl_longRead_wf.git
$ cd methyl_longRead_wf
```

option: install snakemake conda environment
## Run pipeline

option: install snakemake conda environment
```
# conda activate snakemake #activate conda env
#MUST: create `config/samples.yaml` to match paths to per_read_modified_base.txt of each sample
#MUST: edit `config/config.yaml` to set ref genomes

$ sh run_snakefile.sh
```


## inputs
1. `samples.yaml` file specifying sample name and file paths to megalodon per read modified bases file
2. 

## outputs

## QC 
there are a number of qc available
1. calculate depth at each CpG site
```
$ sh samtools_depth_qc.sh '/path/to/megalodon_results_dir' '/path/to/CpG_sites_0based_mm10.bed'
```

### notes on scripts
scripts/R:
aggregate_per_read_5mC_5hmc_main.r  - gather reads per site, calc mean/median methylation for 5mc
methyl_activity_X_gene_expression_and_promoters_pipeline.r - associate gene expressin with promoters
dseq2_analysis.r - dseq2 analysis                    
methylation_activity_gene_expression_and_promoters.r
files_mod_bases2yaml.yaml           
plot_methylation_TransposableElements.r
gene_promoters_with_encode.r        
plot_samtools_depth_coverage.r
make_sample_yaml.r                  
plugNplay_any_granges.r
mergeStats_ONT.r                    
promoter_methylation_knownGenes.r