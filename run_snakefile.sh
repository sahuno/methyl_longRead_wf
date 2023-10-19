#!/bin/bash

# Run snakemake
snakemake --jobname 's.{jobid}.{rulename}' \
	--snakefile $1 \
	--keep-going \
	--reason \
	--printshellcmds \
	--latency-wait 10 \
	--rerun-incomplete \
	--stats logs/stats/snakemake_$(date +"%Y%m%d_%H%M%S").stats \
	-j 10 \
	--cluster-config config/cluster.json \
	--cluster "bsub -q {cluster.queue} -n {cluster.threads} -W {cluster.time} -M{cluster.mem} -R\"span[hosts=1] select[mem>{cluster.mem}] rusage[mem={cluster.mem}]\" {cluster.extra} -o logs/cluster/cluster_$(date +"%Y%m%d_%H%M%S").out -e logs/cluster/cluster_$(date +"%Y%m%d_%H%M%S").err"
# --unlock

#use `-R` to run specific rules
#	-R gather_files \
# -R per_read_aggregate \
#	-R plugNplay_plots_TE \
#	-R gene_promoters \
#	-R per_read_aggregate \
# sh run_snakefile.sh Snakemake_methyl_analysis_from_dorado_modkit.smk

# snakemake -s Snakemake_split_input_by_chrom.smk -np 
# sh run_snakefile.sh Snakemake_split_input_by_chrom.smk
