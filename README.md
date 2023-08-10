# methyl_longRead_wf
a snakemake workflow to compute methylation activity of genes with long read sequencing and associate with gene expression

## how to install pipeline
```
$ git clone git@github.com:sahuno/methyl_longRead_wf.git
$ cd methyl_longRead_wf
# option: install snakemake conda environment 
$ conda activate snakemake #activate conda env

#MUST: edit `config/samples.yaml` to match paths to per_read_modified_base.txt of each sample

$ sh run_snakefile.sh
```
## inputs
1. `samples.yaml` file specifying sample name and file paths to megalodon per read modified bases file
2. 

## QC 
there are a number of qc available
1. depth of sequencing per CpG site
