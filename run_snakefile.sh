#!/bin/bash

# Run snakemake
snakemake --jobname 's.{jobid}.{rulename}' \
	--snakefile Snakefile \
	--keep-going \
	--reason \
	--printshellcmds \
	--latency-wait 10 \
	--rerun-incomplete \
	--stats snakemake_$(date +"%Y%m%d_%H%M%S").stats \
	-j 500 \
	--cluster-config config/cluster.json \
	--cluster "bsub -q {cluster.queue} -n {cluster.threads} -W {cluster.time} -M{cluster.mem} -R\"span[hosts=1] select[mem>{cluster.mem}] rusage[mem={cluster.mem}]\" {cluster.extra} -o out.txt -e err.txt"
# --unlock

#use `-R` to run specific rules
#	-R gather_files \
# -R per_read_aggregate \
#	-R plugNplay_plots_TE \
#	-R gene_promoters \
#	-R per_read_aggregate \

