#get consensus high qulaity cpgs for proportion of disconcordant reads analysis,
#motif reads with at least 5 reads in 75% of samples; this can be changed and we see using read distribution plots
#name: samuel ahuno
#date: sept 1st 2023 

#load library
library(tidyverse)
library(data.table)
library(matrixStats)
library(optparse)
library(dtplyr)

#get inputs from commad line
args <- commandArgs(trailingOnly = TRUE)

# Validate input arguments
if (length(args) != 4) {
  stop("Please provide all four arguments: 1) list of files (comma-separated without spaces), 2) nReads_consensus, 3) percent_Samples_consensus, 4) output directory")
}

# Parse the arguments
file_list <- unlist(strsplit(args[1], split = ","))
nReads_consensus <- as.integer(args[2])
percent_Samples_consensus <- as.double(args[3])
output_dir <- args[4]

print(file_list)
# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Rscript /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/scripts/R/test_r_commandArgs.r "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-1/D-0-1_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-3/D-A-3_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-2/D-A-2_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-2/D-0-2_modpileup_5mC.bed"          5         0.75         results/build_consensus
# Rscript /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/scripts/R/consensus_gp3.R "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-1/D-0-1_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-3/D-A-3_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-2/D-A-2_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed,/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-2/D-0-2_modpileup_5mC.bed"          5         0.75         results/build_consensus
