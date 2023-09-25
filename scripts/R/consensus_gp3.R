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
library(cowplot)

#set
numb_nonSamples_cols <- 3 #number of non sample columns

#get inputs from commad line
args <- commandArgs(trailingOnly = TRUE)

# Check the number of arguments passed
if (length(commandArgs(trailingOnly = TRUE)) < 4) {
  cat("Usage: script.R <nReads_consensus> <percent_Samples_consensus> <output_dir> <file1> [<file2> ...]\n")
  quit("no")
}





# Parse the arguments
# file_list <- unlist(strsplit(args[1], split = ","))
# nReads_consensus <- as.integer(args[2])
# percent_Samples_consensus <- as.double(args[3])
# output_dir <- args[4]
nReads_consensus <- as.integer(args[1])
percent_Samples_consensus <- as.double(args[2])
output_dir <- args[3]
file_list <- args[4:length(args)]

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

##########################################################################################################################################
#use this to test
file_list <- c("/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-1/D-0-1_modpileup_5mC.bed",
        "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-3/D-A-3_modpileup_5mC.bed",
        "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-2/D-A-2_modpileup_5mC.bed",
        "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed",
        "/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-2/D-0-2_modpileup_5mC.bed")

percent_Samples_consensus=0.75
nReads_consensus=5
output_dir="results/build_consensus/"
##########################################################################################################################################

#to do: delete files
# test_files <- list.files("/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit", full.names = TRUE, pattern = "*_modpileup_5mC.bed", recursive = TRUE)
# file_list <- test_files
#get file names
fileNames <- gsub(".bed","",basename(file_list))

#read all files
lapply(file_list, function(x) fread(x, sep = "\t", select=c("V1", "V3", "V6", "V10"))) -> dt_list
message("done reading in bed files as data.table list")


#rename files and remove non standard chroms
names(dt_list) <- fileNames
dt_list_named <- lapply(dt_list, function(x) {
      setnames(x, old = names(x), new = c("chrom", "pos", "strand", "Nvalid_cov"))
      x <- x[!chrom %like% "_"][!chrom %like% "m"]
  return(x)  # Return the modified data.table
    })
message("renamed columns of data.table")


#function to compute stats, over columns
compute_stats_dt <- function(dt_in, x_var){
  dt_stats <- copy(dt_in)[,list(Counts=.N, Mean=mean(get(x_var)), 
                                Max=max(get(x_var)), 
                                Min=min(get(x_var)), 
                                Median=as.numeric(median(get(x_var))), 
                                Std=sd(get(x_var)))]
  return(dt_stats)
}
#stats_D01 <- compute_stats_dt(dt_list_named[[1]], x_var="Nvalid_cov")

##plot distribution of reads per sample
# plt_h <- hist_plt(data_in=dt_list_named[[1]], x_var=Nvalid_cov) + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                   labels = scales::trans_format("log10", scales::math_format(expr = 10^.x)))
# ggsave(plt_h, filename="hist_D-0-1_modpileup_5mC.pdf")

#compute stats over all samples of 
# Nvalid_cov_stats_ls_dt <- lapply(names(dt_list_named),function(x){compute_stats_dt(dt_list_named[[x]], x_var="Nvalid_cov")})

plt_hist_ls <- lapply(names(dt_list_named), function(x) {
  dt_stat <- compute_stats_dt(dt_list_named[[x]], x_var="Nvalid_cov")

  p <- hist_plt(data_in=dt_list_named[[x]], x_var=Nvalid_cov) + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) +
   labs(title = paste0(x)) +
  #  print(Nvalid_cov_stats_ls_dt[[x]])
  # p +
  geom_text(x = 2, y = 4.0e+06, 
            aes(label = paste0("\nCounts: ", Counts,
                              "\nRange: ", Min, "-", Max,
                              "\nmean: ", Mean, 
                              "\nmedian: ", Median,
                              "\nSD: ", Std)
                ),
            data = dt_stat, size=3)
  })

#ggsave(plt_hist_ls[[1]], file="hist_D-0-1_modpileup_5mC_cowPlot.pdf")


#combine plots and save
plt_hist_cow <- plot_grid(plotlist=plt_hist_ls) 
ggsave(plt_hist_cow,filename=file.path(output_dir, "hist_prefiltering_consensus_HighQual_reads_at_motif_sites.pdf"))




#merge all files into 1 data table
# do.call(cbind, dt_list_named) -> dt_cbind

dt_rbindlist <- rbindlist(dt_list_named, idcol="sample") #merge list of data.tables into one data.table
dt_wider <- dcast(copy(dt_rbindlist), formula = chrom + pos + strand ~ sample, value.var = "Nvalid_cov")
message("created wide data.table")
#get dt commands 
# # dt_rbindlist[, sample := gsub(".data.*","",file)] #remove .data.* from filenames
# dt_small <- dt_rbindlist[sample(.N,100)]
# dt_small
# lazy_dt1 <- lazy_dt(dt_small)
# lazy_dt1
# lazy_dt1_wide <- pivot_wider(lazy_dt1, id_cols = c("chrom", "pos", "strand"), names_from = "sample", values_from = "Nvalid_cov")
# names(as.data.frame(lazy_dt1_wide))

##############
#set values for filtering
max_val=15
max_val_low=5
## percentage of samples needed that to support specific cpg sites as high quality
#percThreshold_hiQ_samples <- .75 #
# how many samples needed to support cpg sites as high quality
# minReads_threshold <- 5 #set minimum reads supporting each cpg site
#create toy data
# dt_s <- data.table(chrm = paste0("chrm",1), 
#                    pos = floor(runif(10, min=10000, max=50000)),
#                    sample1 = floor(runif(10, min=0, max=max_val_low)),
#                    sample2 = floor(runif(10, min=0, max=max_val)), 
#                    sample3 = floor(runif(10, min=0, max=max_val)), 
#                    sample4 = floor(runif(10, min=0, max=max_val)))
# dt_s
# dput(dt_s)

#for real data 
# fread(cmd=paste("grep",word,"filename"))
# dtF <- fread(cmd=paste("awk -F'\t' '{print $1, $3, $14}' ",input_file)) #not efficient
# dtF <- fread(cmd=paste("head ",input_file)) #memory efficeint 
# dtF2 <- fread(input_file, select=c("V1", "V3", "V6", "V10", "V11")) #more efficient

# fread("awk -F'\t' '{print $1, $3, $14}' your_file.txt")

# how many samples needed to support cpg sites as high quality
nThreshold_hiQ_samples <- round(percent_Samples_consensus*(ncol(dt_wider)-numb_nonSamples_cols))
message("will filter for motif sites with at least ", nReads_consensus, " reads in ",  nThreshold_hiQ_samples, " of samples")

dt_wider[, num_samples_with_minReads := rowSums(.SD >= nReads_consensus, na.rm=T), .SDcols = !c("chrom", "pos", "strand")]
# Filter CpGs (rows) where 75% of samples have at least 10 (min_reads) reads
consensus_hiQ_CpGs <- dt_wider[num_samples_with_minReads >= nThreshold_hiQ_samples,]
message("done filtering ", nReads_consensus, " reads in ",  nThreshold_hiQ_samples, " of samples")

#find excluded regions
dt_excluded_regions <- dt_wider[!consensus_hiQ_CpGs, on = c("chrom"  ,    "pos" ,"strand")]


# Remove the num_samples_with_minReads column, no longer need it
consensus_hiQ_CpGs[, num_samples_with_minReads := NULL]


#generate stats for consensus high_quality reads
message("computing stats on remaining motif sites after filtering ")
# consensus_hiQ_CpGs[,`:=`(MIN = rowMins(as.matrix(.SD), na.rm=T),SUM=rowSums(.SD, na.rm=T), AVG = rowMeans(.SD, na.rm = T), MED = rowMedians(.SD, na.rm = T)), .SDcols = !c("chrom", "pos", "strand")]
# Get the names of columns that are not in the vector c("chrom", "pos", "strand")
cols_to_include <- setdiff(names(consensus_hiQ_CpGs), c("chrom", "pos", "strand"))

# Use these columns in your data.table operation
consensus_hiQ_CpGs[, `:=`(
  MIN = rowMins(as.matrix(.SD), na.rm=T),
  SUM = rowSums(.SD, na.rm=T),
  AVG = rowMeans(.SD, na.rm = T),
  MED = rowMedians(as.matrix(.SD), na.rm = T)
), .SDcols = cols_to_include]


#save as tsv file
message("writing consensus high quality motif sites to .tsv file")
fwrite(consensus_hiQ_CpGs, file = file.path(output_dir,"consensus_HighQual_motif_sites.tsv"), sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)

#save as bed file
message("generating .bed for consensus high quality motif sites")
fwrite(consensus_hiQ_CpGs[,.(chrom, start=pos-1, pos, strand)], file =file.path(output_dir, "consensus_HighQual_motif_sites.bed") , sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = FALSE)

#save as tsv file
message("writing excluded motif sites as .tsv file")
fwrite(dt_excluded_regions, file = file.path(output_dir,"dt_excluded_motif_sites.tsv"), sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)




#dot plots of methyaltion
# plot_dot_nReads_frac_mod <- ggplot(consensus_hiQ_CpGs, aes(x=fraction_modified, y=Nmod)) + 
#                              geom_point(size=1, alpha=0.3) + labs("Average #Consensus reads per motif site")
#                              # ggsave("plot_dot_nReads_frac_mod.tiff", plot = plot_dot_nReads_frac_mod, width=300, height=225, units="mm", dpi=300, compression = "lzw")


# test_dt2 <- melt(copy(consensus_hiQ_CpGs), id.vars = c("chrom","pos","strand"),
#                 measure.vars = c("AVG", "MED"),
#                 variable.name = "samples", value.name = "Nvalid_cov")



# library(ggplot2)
# library(data.table)
# test_dt <- fread("/juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/results/build_consensus/consensus_HighQual_motif_sites.tsv")
dt_AvgReads <- copy(consensus_hiQ_CpGs)[,list(Counts=.N, Mean=mean(AVG), Max=max(AVG), Min=min(AVG), 
                      Median=as.numeric(median(AVG)), Std=sd(AVG))]

geom_text_size = 3
plt_h <- ggplot(consensus_hiQ_CpGs, aes(x=AVG)) + 
  geom_histogram(bins=50) + 
  labs(title = "Average #Consensus reads per motif site")+
                                        geom_text(x = 2, y = 6.0e+06, 
                                        aes(label = paste0("\nCounts: ", Counts,
                                        "\nRange: ", Min, "-", Max,
                                        "\nmean: ", Mean, 
                                      "\nmedian: ", Median,
                                      "\nSD: ", Std)), 
                              data = dt_AvgReads, size=geom_text_size) +
                              scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(expr = 10^.x)))
ggsave(filename  = file.path(output_dir, "hist_Avg_Consensus_Reads_per_site.tiff"), plot = plt_h, width=300, height=225, units="mm", dpi=300, compression = "lzw")



# #cols to disregard
# omit_col <- paste0(c("chrom","pos","strand", "SUM","AVG", "MIN", "MED"), collapse = "|")
# #melt data.table
# test_dt2 <- melt(consensus_hiQ_CpGs, id.vars = c("chrom","pos","strand"),
#                 measure.vars = names(test_dt)[!str_detect(names(consensus_hiQ_CpGs), omit_col)],
#                 variable.name = "samples", value.name = "Nvalid_cov")

# plt_gw_line <- hDT2 %>% 
#                 ggplot(aes(y=median_exp_mod_log_prob, x=chr_start, color =methylationState )) + 
#                   geom_line(aes(group = 1)) + labs(title = "genomewide methylation")
# ggsave("plot_dot_nReads_frac_mod.tiff", plot = plot_dot_nReads_frac_mod, width=300, height=225, units="mm", dpi=300, compression = "lzw")


        # Rscript {params.R_script_cons} \
        # {params.nReads_build_consen} \
        # {params.Perc_samples_conse} \
        # {output} \
        # {input:q}