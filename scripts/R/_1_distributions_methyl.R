#samuel ahuno
#date: sept 21st 2023


library(data.table)
library(optparse)
library(janitor)
#install.packages("janitor")
#install.packages("tidyverse")
# install.packages("tidyverse")

#get inputs from commad line
option_list_val <- list(opt_str = make_option(c("-i","--input_file"), type="character", default="/lila/data/greenbaum/users/ahunos/apps/dorado_ont_wf/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed", help = "mod base pileup file"),
                        #opt_str = make_option(c("-i","--input_file"), type="character", default="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/addHeader/Spectrum-OV-009_N/data/Spectrum-OV-009_N.chr1.per_read_modified_base_calls.header.txt", help = "per read modified base file"),
                        opt_str = make_option(c("-t","--methyl_threshold"), type="double", default=0.7, help = "threshold to be called methylated between 0-1"),
                        opt_str = make_option(c("--cuttOff_nReads"), type="integer", default=5, help = "# reads needed per site"),
                        opt_str = make_option(c("--data_aggregate_stats"), type="character", default=NULL, help = "5mC and 5hmC output file name .txt"),
                        opt_str = make_option(c("--plot_density_mod_prob"), type="character", default=NULL, help = "density plot _MeanMedian_exp_prob_meth_ONT .pdf"),
                        opt_str = make_option(c("--plot_hist_mod_prob"), type="character", default=NULL, help = "histogram plot MeanMedian_exp_prob_meth_ONT .pdf"),
                        opt_str = make_option(c("--plot_ecdf_mod_prob"), type="character", default=NULL, help = "ecdf plot MeanMedian_exp_prob_meth_ONT.pdf"),
                        opt_str = make_option(c("--plot_density_nReads"), type="character", default=NULL, help = "desnity plot number of reads_meth_ONT.pdf"),
                        opt_str = make_option(c("--plot_hist_nReads"), type="character", default=NULL, help = "histogram plot _number of reads_meth_ONT.pdf"),
                        opt_str = make_option(c("--plot_ecdf_nReads"), type="character", default=NULL, help = "ecdf plot _number of reads _meth_ONT.pdf"),
                        opt_str = make_option(c("--plot_dot_nReadsVrsMedianProb"), type="character", default=NULL, help = "dot plot of exp prob and reads ONT.pdf")#,
#                        opt_str = make_option(c("--plot_GenomeWide"), type="character", default=NULL, help = "Genome wide methylation profile ONT.pdf")
                        )
# opt_str = make_option(c("--methyl_UMH_Stats"), type="character", default=NULL, help = "per_reads_aggregate methylation Stats and States.txt"),

parseObject <- OptionParser(option_list = option_list_val)
opt <- parse_args(parseObject)

########################################################
##check to see if there is a file work with
if (is.null(opt$input_file)){
  print_help(parseObject)
  stop("Please supply sample to work on", call.=FALSE)
}



#load files 
dt <- fread(opt$input_file, sep = "\t")
names(dt) <- c("chrom", "start", "end", "mod_base", "score", "strand", "start2", "end2", "color", "Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")


#1. keep only standard chroms; 2. create a minimal data table fro downstream and ploting
cols_filter <- c("chrom", "start", "end", "strand","mod_base","Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical")
dt_minimal <- dt[!chrom %like% "_"][!chrom %like% "m"][,colnames(dt) %in% cols_filter, with=FALSE]