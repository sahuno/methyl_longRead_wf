configfile: "config/config_split_chroms.yaml"
configfile: "config/samples_automatic.yaml"

import yaml

# Load the samples from the provided YAML configuration file
with open("config/samples_automatic.yaml", 'r') as stream:
    try:
        data = yaml.safe_load(stream)
        SAMPLES = list(data['samples'].values())
    except yaml.YAMLError as exc:
        print(exc)

rule all:
    input:
        "results/build_consensus/consensus_HighQual_motif_sites.bed", 
        expand("results/splitChrom/{samples}/data/{samples}.chr{chr}.bed", samples=config["samples"], chr=config["chrs"]),
        expand("results/sortBed/{samples}/data/{samples}.chr{chr}.sorted.bed", samples=config["samples"], chr=config["chrs"]),
        expand("results/filter_consensus/{samples}/data/{samples}.chr{chr}.filtered_Consensus.bed", samples=config["samples"], chr=config["chrs"]),
        expand("results/promoters/{samples}/data/{samples}.chr{chr}.filtered_Consensus.overlaps.bed", samples=config["samples"], chr=config["chrs"]),
        expand("results/promoter_stats/{samples}/data/{samples}.chr{chr}.methyl_metrics_promoter_data.tsv", samples=config["samples"], chr=config["chrs"]),
        expand("results/gather_files/{samples}/{samples}.methyl_metrics_promoter_data_all_chroms.tsv", samples=config["samples"])


rule build_consensus:
    input:
        SAMPLES
    output:
        consensus_sites="results/build_consensus/consensus_HighQual_motif_sites.bed",
        dir_out=directory("results/build_consensus")
    params:
        R_script_cons=config["consensus_sites"],
        Perc_samples_conse=config["P_samples_conse"],
        nReads_build_consen=config["nReads_conse"],
        dir="results/build_consensus/"
    log:
        "logs/build_consensus/build_consensus.log"
    shell:
        """
        Rscript {params.R_script_cons} \
        {params.nReads_build_consen} \
        {params.Perc_samples_conse} \
        {output.dir_out} \
        {input:q} \
        &> {log}
        """

rule splitChrom:
    input:
        input_file_from_config=lambda wildcards: config["samples"][wildcards.samples]
    output:
        unsorted_beds="results/splitChrom/{samples}/data/{samples}.chr{chr}.bed"
    params:
        split_chr_script_mod=config["splitChr_script"],
        dir="results/splitChrom/{samples}/data"
    # log:
    #     "logs/splitChrom/{samples}/{samples}.chr{chr}.log"
    shell:
        """
        bash {params.split_chr_script_mod} {wildcards.samples} {input.input_file_from_config} {params.dir}
        """

rule sortBed:
    input:
        inp_bed="results/splitChrom/{samples}/data/{samples}.chr{chr}.bed"
    output:
        sorted_beds="results/sortBed/{samples}/data/{samples}.chr{chr}.sorted.bed"
    log:
       "logs/sortBed/{samples}.chr{chr}.log"
    shell:
        """
        /bin/bash -c 'sort -k1,1 -k2,2n {input.inp_bed} > {output.sorted_beds} 2> {log}'
        """

# #filter good regions only 
rule filter_consensus:
    input:
        methyl_chrom_bed="results/splitChrom/{samples}/data/{samples}.chr{chr}.bed",
        consensus_motif_sites="results/build_consensus/consensus_HighQual_motif_sites.bed"
    output:
        sample_consensus_motif_sites="results/filter_consensus/{samples}/data/{samples}.chr{chr}.filtered_Consensus.bed"
    # params:
    #     consensus_motif_sites="results/build_consensus/consensus_HighQual_motif_sites.bed"
    log:
        "logs/filter_consensus/{samples}/{samples}.chr{chr}.log"
    shell:
        """
        bedtools intersect -a {input.methyl_chrom_bed} -b {input.consensus_motif_sites} > {output.sample_consensus_motif_sites} 2> {log[0]}
        """

rule promoters:
    input:
        input_bed="results/filter_consensus/{samples}/data/{samples}.chr{chr}.filtered_Consensus.bed"
    output:
        output_bed="results/promoters/{samples}/data/{samples}.chr{chr}.filtered_Consensus.overlaps.bed"
    params:
        Bedtools_ov=config["intersect_script"],
        prom_bed=config["path_promoter_bed"]
    log:
      "logs/promoters/{samples}/{samples}.chr{chr}.log"
    shell:
        """
bedtools intersect -a {input.input_bed} -b {params.prom_bed} -wb > {output.output_bed} 2> {log}
          """

# bash {params.Bedtools_ov} {input.input_bed} {params.prom_bed} {output.output_bed} 1> {log}


rule promoter_stats:
    input:
        ov_promoters="results/promoters/{samples}/data/{samples}.chr{chr}.filtered_Consensus.overlaps.bed"
    output:
        ov_stats="results/promoter_stats/{samples}/data/{samples}.chr{chr}.methyl_metrics_promoter_data.tsv"
    params:
        R_script=config["prom_stats"],
        methyl_cutoff=config["methyl_thresh"],
        nreads_thresh=config["nReads_conse"]
    log:
        "logs/promoter_stats/{samples}.chr{chr}.log"
    shell:
        """
        Rscript {params.R_script} \
        --input_file {input.ov_promoters} \
        --nReads {params.nreads_thresh} \
        --methyl_percent {params.methyl_cutoff} \
        --methyl_metrics_promoter_data {output.ov_stats} &> {log}
        """

rule gather_files:
    input:
        prom_ovs=expand("results/promoter_stats/{samples}/data/{samples}.chr{chr}.methyl_metrics_promoter_data.tsv", chr=config["chrs"], allow_missing=True)
        # methyl_rate_prom_data=expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_metrics.txt", chr=config["chrs"], allow_missing=True)
    output:
        prom_stats_merge="results/gather_files/{samples}/{samples}.methyl_metrics_promoter_data_all_chroms.tsv",
    params:
        gatherFilesScript=config["combineStats_files"],
        awk_arg=r"""'BEGIN {{FS=OFS="\t"}} NR==1 || FNR!=1'"""
    log:
       "logs/gather_files/{samples}.log"
    run:
        shell("awk {params.awk_arg} {input.prom_ovs} > {output.prom_stats_merge}")

    # shell:
    #     """
    #     HEADER_prom=$(head -n 1 {input.prom_ovs[0]})  # Use the header from the first input file
    #     echo -e $HEADER_prom > {output.prom_stats_merge}
    #     for file in {input.prom_ovs}; do
    #         if [ "$file" != "{input.prom_ovs[0]}" ]; then
    #             awk {params.awk_arg} "$file" >> {output.prom_stats_merge}
    #         fi
    #     done 2> {log}
    #     """



# snakemake -s Snakemake_split_input_by_chrom.smk -np
#sh run_snakefile.sh Snakemake_split_input_by_chrom.smk



# rule makeBedFiles:
#     input:
#         in_bedMethyl="results/splitChrom/{samples}/data/{samples}.chr{chr}.bed"
#     output:
#         out_bed="results/makeBedFiles/{samples}/data/{samples}.chr{chr}.bed"
#     log:
#        "logs/makeBedFiles/{samples}.chr{chr}.log"
#     shell:
#         """
#         awk 'BEGIN {{ OFS="\t" }}{{print $1, $2, $3, $6, $4, $5,$7,$8,$9, $10, $11, $12, $13, $14, $15, $16, $17, $18}}' {input.in_bedMethyl} > {output.out_bed}
#         """



        # # expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters.pdf", samples=config["samples"], chr=config["chrs"]),
        # # expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_plots.rds", samples=config["samples"], chr=config["chrs"]),
        # #expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.data_methylation_rate_promo_stats.txt", samples=config["samples"], chr=config["chrs"]),
        # expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_metrics.txt", samples=config["samples"], chr=config["chrs"]),
        # expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_data.rds", samples=config["samples"], chr=config["chrs"])
        # expand("results/makeBedFiles/{samples}/data/{samples}.chr{chr}.bed", samples=config["samples"], chr=config["chrs"]),


                # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-1/D-0-1_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-3/D-A-3_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-2/D-A-2_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-2/D-0-2_modpileup_5mC.bed"
        # aggregate_stat_data=expand("results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt", chr=config["samples"], allow_missing=True),


##QC
# remove any spaces in the bed file
# tr ' ' '\t' < results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted.bed > results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted_spacesRemov.bed
# grep ' ' results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted.bed # check where there are possible spaces

#give me non integers values in the second and column, all passes
# awk -F'\t' '($2 !~ /^[0-9]+$/) {print $2}' results/build_consensus/consensus_HighQual_motif_sites.bed
# awk -F'\t' '($2 !~ /^[0-9]+$/) {print $2}' results/build_consensus/consensus_HighQual_motif_sites.bed
# awk -F'\t' '($2 !~ /^[0-9]+$/) {print $2}' results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted.bed
# awk -F'\t' '($3 !~ /^[0-9]+$/) {print $3}' results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted.bed


# results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted_spacesRemov.bed
# bedtools intersect -a results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted_spacesRemov.bed -b results/build_consensus/consensus_HighQual_motif_sites.bed > results/filter_consensus/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.filtered_Consensus.bed
# awk '{print NF}' file | sort -nu | tail -n 1
# awk '{print NF}' /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr1.sorted.bed | sort -nu | tail -n 1
# 
# awk '{print NF}' results/build_consensus/consensus_HighQual_motif_sites.bed | sort -nu | tail -n 1
# cat results/build_consensus/consensus_HighQual_motif_sites.bed | cut -f 1| sort | uniq -c | sort -nr | head -n 1
