#date Jun 2nd 2023
#purpose; agrregate per read methylation statistics 

library(data.table)
library(magrittr)
library(optparse)
library(tidyverse)
library(dtplyr)
library(janitor)
library(pryr)

########################################################
######accept arguments from environment
########################################################
option_list_val <- list(#opt_str = make_option(c("-i","--input_file"), type="character", default="/juno/work/greenbaum/users/ahunos/rotation/scripts/workflows/qc_DMR/results/addHeader/BRCA_13135_P_2/data/BRCA_13135_P_2.chr18.per_read_modified_base_calls.header.txt", help = "per read modified base file"),
                        opt_str = make_option(c("-i","--input_file"), type="character", default="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/addHeader/Spectrum-OV-009_N/data/Spectrum-OV-009_N.chr1.per_read_modified_base_calls.header.txt", help = "per read modified base file"),
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


#####load per read modified bed befile
#opt$input_file <- "/work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon/Spectrum-OV-044_N/per_read_modified_base_calls.txt"

message("reading files")
df_perReadModifiedBases <- data.table::fread(opt$input_file, nThread=16)
#df_perReadModifiedBases <- data.table::fread("/work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon/Spectrum-OV-009_T/per_read_modified_base_calls.txt", nThread=16)



###################################
# aggregate reads + summary stats
###################################
#df_perReadModifiedBases <- fread("/juno/work/greenbaum/users/ahunos/rotation/scripts/workflows/qc_DMR/results/addHeader/BRCA_13135_P_2/data/BRCA_13135_P_2.chr18.per_read_modified_base_calls.header.txt")
message("aggregating reads for 5mc only")
df_perReadModifiedBases <- df_perReadModifiedBases[!chrm %like% "_"][mod_base %like% "m"][,
                        `:=`(exp_mod_log_prob = exp(mod_log_prob))][,
                        list(PrM = sum(exp_mod_log_prob)),
                        by=list(chrm, pos, strand, read_id)][,
                        list(number_reads=.N, mean_Prob_Meth_perSite = mean(PrM), median_Prob_Meth_perSite = median(PrM)),by=list(chrm, pos, strand)]


message("size aggregated reads summary stats")
dim(df_perReadModifiedBases)

df_perReadModifiedBases[,.N,by=number_reads]
df_perReadModifiedBases[,.(medReads=median(number_reads),meanReads=mean(number_reads))]
#df_perReadModifiedBases
#lobstr::obj_size(df_perReadModifiedBases)
message("save long table to disk")
#write to disk as tsv 
fwrite(df_perReadModifiedBases, file = opt$data_aggregate_stats, sep="\t")
# fwrite(df_perReadModifiedBases, file = "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/addHeader/Spectrum-OV-009_N/data/Spectrum-OV-009_N.chr1.per_read_modified_base_calls.header.txt", 
# sep="\t")

## convert file into long to plot
df_perReadModifiedBases <- melt(`df_perReadModifiedBases`, measure.vars = c("mean_Prob_Meth_perSite", "median_Prob_Meth_perSite"), variable.name = "statistic", variable.factor = FALSE)


############################
### plots
############################
#plots mean and median methylation
plt_density_stats_meth <- df_perReadModifiedBases %>% ggplot(aes(x=value, fill = statistic)) + 
                                                geom_density(alpha = 0.5) + 
                                                theme(legend.position="bottom") + 
                                                labs(title = "density of mean and median of prob of methylated aggregated per site")
ggsave(plt_density_stats_meth, file=opt$plot_density_mod_prob)

plt_histo_stats_meth <- df_perReadModifiedBases %>% ggplot(aes(x=value, fill = statistic)) + 
                  geom_histogram(alpha = 0.5, position="identity") +
                  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) + 
                  labs(title = "histogram of mean and median of prob of methylated aggregated per site")+
                  theme(legend.position="bottom")

ggsave(plt_histo_stats_meth, file=opt$plot_hist_mod_prob)


plt_ecdf_stats_meth <- df_perReadModifiedBases %>% 
                            ggplot(aes(x=value, color = statistic)) + 
                                stat_ecdf() + 
                                    labs(title = "ecdf of mean and median of prob of methylated aggregated per site") +
                                        theme(legend.position="bottom")
ggsave(plt_ecdf_stats_meth, file=opt$plot_ecdf_mod_prob)

#plots counts of number of reads methylation
plt_density_stats_reads <- df_perReadModifiedBases %>% ggplot(aes(x=number_reads)) + 
                                                geom_density(alpha = 0.5) + 
                                                theme(legend.position="bottom") +
                                                scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                                labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) + 
                                                labs(title = "Density of number of aggregated reads per site") 
ggsave(plt_density_stats_reads, file=opt$plot_density_nReads)


plt_histo_stats_reads <- df_perReadModifiedBases %>% ggplot(aes(x=number_reads)) + 
                  geom_histogram(alpha = 0.5, position="identity") +
                  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) + 
                  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) +
                  theme(legend.position="bottom") + 
                  labs(title = "histogram of number of aggregated reads per site") 
ggsave(plt_histo_stats_reads, file=opt$plot_hist_nReads)


plt_ecdf_stats_reads <- df_perReadModifiedBases %>% 
                ggplot(aes(x=number_reads)) + 
                    stat_ecdf() + 
                    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) +
                    labs(title = "ecdf of number of aggregated reads per site") +
                    theme(legend.position="bottom")  
ggsave(plt_ecdf_stats_reads, file=opt$plot_ecdf_nReads)


#############################
# message("save copy of data")
# fwrite(df_perReadModifiedBases, file ="debug_methyl.tsv", sep="\t")







#dot plots of methyaltion
plot_point_nReadsMeanProb <- df_perReadModifiedBases %>% 
                            ggplot(aes(x=value, y=number_reads)) + 
                             geom_point(size=1, alpha=0.3) +
scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) +
                  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999"))  + facet_wrap(~statistic)
            #scale_x_continuous(breaks = seq(0, 1, 0.1))
ggsave(file=opt$plot_dot_nReadsVrsMedianProb, plot = plot_point_nReadsMeanProb, width=9, height=7)



#TODO
#save each plot as ggplot object, so we can gather all individual chrom files as one gaint file for visualization
#make a dir to dump figures



############################################
#plot per genome methylation
df_perReadModifiedBases[,`:=`(chr_start=paste0(chrm,"_",pos))]
setorder(df_perReadModifiedBases, cols = "chr_start")     # Reorder data.table

message("plot number of reads and methylation probabilities")
# plt_gw <- df_perReadModifiedBases[number_reads >= opt$cuttOff_nReads, ] %>% 
#           ggplot(aes(y=value, x=chr_start)) + 
#               geom_point(size=1, alpha=0.3) + 
#              scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999"))  + 
#               geom_smooth(method="gam", formula = y~s(x)) +
#               labs(title = paste0("Chromosome wide methylation, #reads cutOff = ", opt$cuttOff_nReads)) + 
#               scale_x_discrete(name ="chr_start") + 
#               theme(axis.ticks.x=element_blank(),axis.text.x=element_blank()) + 
#               facet_wrap(~statistic)

#ggsave(opt$plot_GenomeWide, plot = plt_gw, width=9, height=7)


message("successfully completed script")