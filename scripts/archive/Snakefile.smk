configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input: 
        expand("results/splitChrom/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/addHeader/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.header.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/data/{samples}.chr{chr}.methyl_UMH_Stats.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/data/{samples}.chr{chr}.data_aggregate_stats.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.density_plot_mod_prob.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.histogram_mod_prob.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.ecdf_plot_mod_prob.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.density_plot_nReads.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.histogram_nReads.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.ecdf_plot_nReads.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.dot_plot_nReadsVrsMedianProb.pdf", samples=config["samples"], chr=config["chrs"]),
        expand("results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.GenomeWide_methylation.pdf", samples=config["samples"], chr=config["chrs"])


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
       #rm {output}
       """

rule per_read_aggregate:
    input:
        infile="results/addHeader/{samples}/data/{samples}.chr{chr}.per_read_modified_base_calls.header.txt"
    output:
        MStats_UMH="results/per_read_aggregate/{samples}/data/{samples}.chr{chr}.methyl_UMH_Stats.txt",
        MStats="results/per_read_aggregate/{samples}/data/{samples}.chr{chr}.data_aggregate_stats.txt",
        densPlot_meth="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.density_plot_mod_prob.pdf",
        hist_meth="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.histogram_mod_prob.pdf",
        ecdf_meth="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.ecdf_plot_mod_prob.pdf",
        densPlot_nReads="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.density_plot_nReads.pdf",
        hist_nReads="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.histogram_nReads.pdf",
        ecdf_nReads="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.ecdf_plot_nReads.pdf",
        plot_nreads_meth="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.dot_plot_nReadsVrsMedianProb.pdf",
        genomewide_meth="results/per_read_aggregate/{samples}/figures/{samples}.chr{chr}.GenomeWide_methylation.pdf"

    params:
        R_script=config["qc_Script"]
    log:
        "logs/per_read_aggregate/{samples}/{samples}.chr{chr}.log"
    shell:
       """
        (Rscript {params.R_script} --input_file {input.infile} \
        --methyl_UMH_Stats {output.MStats_UMH} \
        --data_aggregate_stats {output.MStats} \
        --plot_density_mod_prob {output.densPlot_meth} \
        --plot_hist_mod_prob {output.hist_meth} \
        --plot_ecdf_mod_prob {output.ecdf_meth} \
        --plot_density_nReads {output.densPlot_nReads} \
        --plot_hist_nReads {output.hist_nReads} \
        --plot_ecdf_nReads {output.ecdf_nReads} \
        --plot_dot_nReadsVrsMedianProb {output.plot_nreads_meth} \
        --plot_GenomeWide {output.genomewide_meth} 2> {log}
        """


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



# this runs ok
#conda activate snakemake
# cd qc_DMR
#sh run_snakefile.sh 
# snakemake -np #test run with
# snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores
# nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &


