configfile: "config/config.yaml"
configfile: "config/samples_automatic.yaml"

import yaml

# Load the samples from the provided YAML configuration file
with open("config/samples_automatic.yaml", 'r') as stream:
    try:
        data = yaml.safe_load(stream)
        SAMPLES = list(data['samples'].values())
    except yaml.YAMLError as exc:
        print(exc)



# rule all:
#     input:
#         "results/build_consensus"     
rule build_consensus:
    input:
        SAMPLES
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-1/D-0-1_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-3/D-A-3_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-2/D-A-2_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed",
        # "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-2/D-0-2_modpileup_5mC.bed"
        # aggregate_stat_data=expand("results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt", chr=config["samples"], allow_missing=True),
    output:
        directory("results/build_consensus")
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
        {output} \
        {input:q} \
        &> {log}
        """

#        echo "Building consensus for {input:q}..."
# sh run_snakefile.sh Snakemake_methyl_analysis_from_dorado_modkit.smk
#snakemake -s Snakemake_methyl_analysis_from_dorado_modkit.smk -np