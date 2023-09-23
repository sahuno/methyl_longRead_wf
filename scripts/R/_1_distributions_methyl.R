#samuel ahuno
#date: sept 21st 2023


#install.packages("plyranges")
library(data.table)
library(optparse)
library(janitor)
library(ggplot2)
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
dt_minimal[,`:=`(pos = paste0(chrom,":",end,strand))] #create id cols
dt_minimal[,`:=`(fraction_modified = fraction_modified/100)] #use fractions instead of percentages
#sample few data points
#dt_s <- dt_minimal[sample(.N,100)]

# dput(head(dt_s, 25))
# dt_s[,`:=`(code = fcase(fraction_modified > opt$methyl_threshold, mod_base, rep_len(TRUE, length(mod_base)), fraction_modified))]

#make methyl cutt-off vals
message("using ", opt$methyl_threshold, " as lower bound")
methyl_thresh <- seq(from = opt$methyl_threshold, to = 1, by=0.05)
#create a for loop for different methylation c(0.6, 0.07, 0.8,0.) values

#function to assign methyl state
assign_methyl_state <- function(dt_func, threshold_value) {
  dt_copy <- copy(dt_func)[, `:=`(methyl_state = fcase(fraction_modified >= threshold_value, mod_base,
                             rep(TRUE, .N), "U"))]
  return(dt_copy)
}

#find m and U with different trhesholds
dt_ls <- lapply(methyl_thresh, assign_methyl_state, dt_func=dt_minimal)
ls_names_to_append <- paste0("methyl_threshold_", methyl_thresh)
names(dt_ls) <- paste0("methyl_threshold_", methyl_thresh) #add headers to the list
# dt_min_Mstate <- assign_methyl_state(dt_minimal, opt$methyl_threshold) #test case for single threshold
saveRDS(dt_ls, "sample_methyl_bed.rds")


#############
############################
### plots
############################
#plots mean and median methylation
# str(head(dt_ls[[1]]))
plt_density_stats_meth <-  ggplot(dt_ls[[1]], 
                                    aes(x=fraction_modified, color = methyl_state)) + 
                                  geom_density() + 
                                  theme(legend.position="bottom") + 
                                  # facet_wrap(~methyl_state, scales = "free_x")
                                  labs(title = paste0("density of proportions motif reads at genomic site", ls_names_to_append[1]))
ggsave(plt_density_stats_meth, file="plot_density_mod_prop.pdf")

plt_histo_stats_meth <- ggplot(dt_ls[[1]], aes(x=fraction_modified, color = methyl_state)) + 
                  geom_histogram(alpha = 0.5, position="identity") +
                  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) + 
                  labs(title = paste0("hist of proportions motif reads at genomic site"))+
                  theme(legend.position="bottom")
ggsave(plt_histo_stats_meth, file="plot_hist_mod_prop.pdf")


plt_ecdf_stats_meth <- ggplot(dt_ls[[1]], aes(x=fraction_modified, color = methyl_state)) + 
                                stat_ecdf() + 
                                    labs(title = paste0("ecdf")) +
                                        theme(legend.position="bottom")
ggsave(plt_ecdf_stats_meth, file="plot_ecdf_mod_prob.pdf")



#dot plots of methyaltion
plot_dot_nReads_frac_mod <- ggplot(dt_ls[[1]], aes(x=fraction_modified, y=Nmod)) + 
                             geom_point(size=1, alpha=0.3) +
scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) +
                  # scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999"))  
                  facet_wrap(~methyl_state)
            #scale_x_continuous(breaks = seq(0, 1, 0.1))
ggsave(file="plot_dot_nReads_frac_mod.pdf", plot = plot_dot_nReads_frac_mod, width=9, height=7)
