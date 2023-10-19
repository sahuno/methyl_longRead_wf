# samuel ahuno


# load libraries
library(data.table)
library(DESeq2)

data_path <- "/juno/work/greenbaum/projects/TRI_EPIGENETIC/RNASeq_DE_TriEpi/data/condition_AZCT_vs_DMSO/"
files_list <- list.files(data_path, full.names = TRUE)

path_de <- "/juno/work/greenbaum/projects/TRI_EPIGENETIC/RNASeq_DE_TriEpi/data/condition_AZCT_vs_DMSO/Dseq2ResultsObject_condition_AZCT_vs_DMSO_padjust.rds"


# load methylation data

dt_methyl <- fread("/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/results/gather_files/D-0-1_modpileup_5mC/D-0-1_modpileup_5mC.methyl_metrics_promoter_data_all_chroms.tsv")
