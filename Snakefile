configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input: 
        expand("results/splitChrom/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/addHeader/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.header.txt", samples=config["samples"], chr=config["chrs"]),
        # expand("results/per_read_aggregate/{samples}/data/{samples}.chr{chr}.methyl_UMH_Stats.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.density_plot_mod_prob.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.histogram_mod_prob.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.ecdf_plot_mod_prob.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.density_plot_nReads.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.histogram_nReads.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.ecdf_plot_nReads.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.dot_plot_nReadsVrsMedianProb.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_plots.rds", samples=config["samples"], chr=config["chrs"]),
        #expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.data_methylation_rate_promo_stats.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_metrics.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_data.rds", samples=config["samples"], chr=config["chrs"]),
        expand("results/TEs/{samples}/chr{chr}/{samples}.chr{chr}.methylation_TEs_data.rds", samples=config["samples"], chr=config["chrs"]),
        expand("results/TEs/{samples}/chr{chr}/{samples}.chr{chr}.methylation_TEs.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/TEs/{samples}/chr{chr}/{samples}.chr{chr}.methylation_TEs.rds", samples=config["samples"], chr=config["chrs"]),
        expand("results/gather_files/{samples}/{samples}.data_aggregate_stats_all_chroms.txt", samples=config["samples"]),
        expand("results/gather_files/{samples}/{samples}.data_methylation_gene_promoters_metrics_all_chroms.txt", samples=config["samples"]),
        expand("results/plugNplay_plots_TE/{samples}/{samples}.methylation_TEs.pdf", samples=config["samples"]),
        expand("results/plugNplay_plots_TE/{samples}/{samples}.methylation_TEs.rds",samples=config["samples"])

rule splitChrom:
    input:
        input_file_from_config=lambda wildcards: config["samples"][wildcards.samples]
    output:
        #directory("results/splitChrom/{samples}"),
        output_file_user_sets=expand("results/splitChrom/{{samples}}/data/{{samples}}.chr{chr}.per_read_modified_base_calls.txt", chr=config["chrs"])
    params:
        split_chr_script=config["splitChr_script"],
        # chrs=config["chrs"],
        dir="results/splitChrom/{samples}/data"
    log:
        "logs/splitChrom/{samples}/{samples}.log"
    shell:
        "bash {params.split_chr_script} {wildcards.samples} {input.input_file_from_config} {params.dir} &> {log}"

rule addHeader:
    input:
        "results/splitChrom/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.txt"
    output:
        "results/addHeader/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.header.txt"
    log:
        "logs/addHeader/{samples}/{samples}.chr{chr}.log"
    shell:
       """
       (printf 'read_id\tchrm\tstrand\tpos\tmod_log_prob\tcan_log_prob\tmod_base\n' && cat {input}) > {output} #&& mv {output} {input}
       echo "done adding header to methylation file"
       """

rule per_read_aggregate:
    input:
        infile="results/addHeader/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.header.txt"
    output:
        MStats="results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt",
        densPlot_meth="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.density_plot_mod_prob.pdf",
        hist_meth="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.histogram_mod_prob.pdf",
        ecdf_meth="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.ecdf_plot_mod_prob.pdf",
        densPlot_nReads="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.density_plot_nReads.pdf",
        hist_nReads="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.histogram_nReads.pdf",
        ecdf_nReads="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.ecdf_plot_nReads.pdf",
        plot_nreads_meth="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.dot_plot_nReadsVrsMedianProb.pdf"#,
        #genomewide_meth="results/per_read_aggregate/{samples}/chr{chr}/figures/{samples}.chr{chr}.GenomeWide_methylation.pdf"

    params:
        R_script=config["qc_Script"]
    log:
        "logs/per_read_aggregate/{samples}/{samples}.chr{chr}.log"
    shell:
       """
Rscript {params.R_script} --input_file {input.infile} \
--data_aggregate_stats {output.MStats} \
--plot_density_mod_prob {output.densPlot_meth} \
--plot_hist_mod_prob {output.hist_meth} \
--plot_ecdf_mod_prob {output.ecdf_meth} \
--plot_density_nReads {output.densPlot_nReads} \
--plot_hist_nReads {output.hist_nReads} \
--plot_ecdf_nReads {output.ecdf_nReads} \
--plot_dot_nReadsVrsMedianProb {output.plot_nreads_meth} \
2> {log}
        """
#--plot_GenomeWide {output.genomewide_meth} \
# --methyl_UMH_Stats {output.MStats_UMH} \

rule gene_promoters:
    input:
        infile="results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt"
    output:
        data_overlaps="results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_data.rds",
        data_methyl_metrics_prom="results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_metrics.txt",
        plot_pdf="results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters.pdf",
        plot_rds="results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_plots.rds"
    params:
        promoter_scripts=config["promoter_plotsR"],
        #HgVersion="latest",
        stats="median",
        minReads=3,
        Methyl_cutoff=0.90,
        results_dir="results/gene_promoters/{samples}/chr{chr}"
        # chromosomes=lambda wildcards: config["chrs"][wildcards.chrs],
    log:
       "logs/gene_promoters/{samples}/{samples}.chr{chr}.log"
    shell:
       """
Rscript {params.promoter_scripts} \
        --stats2Use {params.stats} \
        --nReads {params.minReads} \
        --methyl_percent {params.Methyl_cutoff} \
        --input_file {input.infile} \
        --plots_pdf {output.plot_pdf} \
        --plots_rds {output.plot_rds} \
        --chrom {wildcards.chr} \
        --overlaps_rds {output.data_overlaps} \
        --methyl_metrics_promo_data {output.data_methyl_metrics_prom} &> {log}
        """

rule TEs:
    input:
        infile="results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt"
    output:
        data_rds="results/TEs/{samples}/chr{chr}/{samples}.chr{chr}.methylation_TEs_data.rds",
        plot_pdf="results/TEs/{samples}/chr{chr}/{samples}.chr{chr}.methylation_TEs.pdf",
        plot_rds="results/TEs/{samples}/chr{chr}/{samples}.chr{chr}.methylation_TEs.rds"
    params:
        TE_scripts=config["TE_plotsR"],
        # chromosomes=lambda wildcards: config["chrs"][wildcards.chrs],
        results_dir="results/TEs/{samples}/chr{chr}",
        TE_file="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/hsflil1_8438.bed"
    log:
       "logs/TEs/{samples}/{samples}.chr{chr}.log"
    shell:
       """
Rscript {params.TE_scripts} --nReads 3 --input_file {input.infile} --plots_pdf {output.plot_pdf} --plots_rds {output.plot_rds} --chrom {wildcards.chr} --output_rds {output.data_rds} --TE_ref {params.TE_file} &> {log}
        """

#from subprocess import run
##############
##############
##TODOS: gather all files
rule gather_files:
    input:
        aggregate_stat_data=expand("results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt", chr=config["chrs"], allow_missing=True),
        methyl_rate_prom_data=expand("results/gene_promoters/{samples}/chr{chr}/{samples}.chr{chr}.methylation_gene_promoters_metrics.txt", chr=config["chrs"], allow_missing=True)
    output:
        aggregate_stat_merged="results/gather_files/{samples}/{samples}.data_aggregate_stats_all_chroms.txt",
        merged_methyl_rate_data="results/gather_files/{samples}/{samples}.data_methylation_gene_promoters_metrics_all_chroms.txt"
    params:
        gatherFilesScript=config["combineStats_files"],
        awk_arg=r"""-vFS="\t" -vOFS="\t" FNR!=1"""
    log:
       "logs/gather_files/{samples}.log"
    run:
       shell("awk {params.awk_arg} {input.aggregate_stat_data} > {output.aggregate_stat_merged}")
       shell("awk {params.awk_arg} {input.methyl_rate_prom_data} > {output.merged_methyl_rate_data}")
#    shell: "cat {input} > {output} 2> {log}" #this works but you end up with 


rule plugNplay_plots_TE:
    input:
        infile="results/gather_files/{samples}/{samples}.data_aggregate_stats_all_chroms.txt"
    output:
        plot_pdf="results/plugNplay_plots_TE/{samples}/{samples}.methylation_TEs.pdf",
        plot_rds="results/plugNplay_plots_TE/{samples}/{samples}.methylation_TEs.rds"
    params:
        TE_pp_scripts=config["plug_play_plotsTE_R"],
        bedfiles="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/chr8_TE.bed"
    log:
       "logs/plugNplay_plots_TE/{samples}.log"
    shell:
       """
Rscript {params.TE_pp_scripts} --nReads 3 --bedfile {params.bedfiles} --statsFile {input.infile} --plots_pdf {output.plot_pdf} --plots_rds {output.plot_rds} --extend_by 4000 &> {log}
        """


# this runs ok
#conda activate snakemake
# cd qc_DMR
#sh run_snakefile.sh 
# snakemake -np #test run with
# snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores
# nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &





# rule two
# rule qc_and_DMR:
#     input:
#         tum="results/splitChrom/{samples}/{samples}.chr{chr}.txt"
#         #norm="results/splitChrom/{samples}/{samples}.chr{chr}.txt",
#     output:
#         data="results/qc_and_DMR/{samples}/data/{samples}.chr{chr}.rds",

#         # data="results/qc_and_DMR/{samples}/plots/{samples}.chr{chr}.pdf",
#         # data="results/qc_and_DMR/{samples}/{samples}.chr{chr}.txt",
#         # data="results/qc_and_DMR/{samples}/{samples}.chr{chr}.txt",
#         # data="results/qc_and_DMR/{samples}/{samples}.chr{chr}.txt",
#     params:
#         qc_DMRScript=config["qc_DMRScript"],
#         # chrs=config["chrs"],
#         dir="results/splitChrom/{samples}"
#     log:
#         "logs/splitChrom/{samples}/{samples}.chr{chr}.log"
#     shell:
#         "bash {params.split_chr_script} {wildcards.samples} {input.input_file_from_config} {params.dir} &> {log}"

