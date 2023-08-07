configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input: 
        expand("results/splitChrom/{samples}/{samples}.chr{chrs}.txt", samples=config["samples"], chrs=str(config["chrs"]))

rule splitChrom:
    input:
        input_file_from_config=lambda wildcards: config["samples"][wildcards.samples]
    output:
        output_file_user_sets="results/splitChrom/{samples}/{samples}.chr{chrs}.txt"
    params:
        split_chr_script=config["splitChr_script"],
        chrs=config["chrs"],
        dir="results/splitChrom/{samples}"
    log:
        "logs/splitChrom/{samples}/{samples}.chr{chrs}.log"
    shell:
        "bash {params.split_chr_script} {wildcards.samples} {input.input_file_from_config} {params.dir} 2> {log}"




#rule two
# rule qc_and_DMR:
#     input:
#         tum="results/splitChrom/{samples}/{samples}.chr"+ str(config["chrs"]) + ".txt"
#         norm="results/splitChrom/{samples}/{samples}.chr"+ str(config["chrs"]) + ".txt"
#     output:
#         data="results/qc_and_DMR/{samples}/{samples}.chr{chrs}.txt",
#         data="results/qc_and_DMR/{samples}/{samples}.chr{chrs}.txt",
#         data="results/qc_and_DMR/{samples}/{samples}.chr{chrs}.txt",
#         data="results/qc_and_DMR/{samples}/{samples}.chr{chrs}.txt",
#         data="results/qc_and_DMR/{samples}/{samples}.chr{chrs}.txt",
        
#     params:
#         split_chr_script=config["splitChr_script"],
#         chrs=config["chrs"],
#         dir="results/splitChrom/{samples}"
#     log:
#         "logs/splitChrom/{samples}/{samples}.chr{chrs}.log"
#     shell:
#         "bash {params.split_chr_script} {wildcards.samples} {input.input_file_from_config} {params.dir} 2> {log}"



# this runs ok
# cd qc_DMR
# snakemake -np #test run with
# snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores
# nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &


