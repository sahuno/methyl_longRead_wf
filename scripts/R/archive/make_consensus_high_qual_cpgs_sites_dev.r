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
option_list_val <- list(#opt_str = make_option(c("-i","--input_file"), type="character", default="/lila/data/greenbaum/users/ahunos/apps/dorado_ont_wf/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed", help = "mod base pileup file"),
                        opt_str = make_option(c("-i","--input_files"), type="character", default="/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-1/D-A-1_modpileup_5mC.bed", help = "juno file"),
                        opt_str = make_option(c("-o","--output_dir"), type="character", default="consensus_motif_sites_modpileup_5mC.bed", help = "juno file"),
                        opt_str = make_option(c("--percent_Samples_consensus"), type="double", default=0.75, help = "percentage of sample needed to support a cpg site"),
                        opt_str = make_option(c("--nReads_consensus"), type="integer", default=5, help = "set a number based on read distribution plots")
                        )

parseObject <- OptionParser(option_list = option_list_val)
opt <- parse_args(parseObject)


########################################################
##check to see if there is a file work with
if (is.null(opt$input_files)){
  print_help(parseObject)
  stop("Please supply sample to work on", call.=FALSE)
}

#to do: delete files
# test_files <- list.files("/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit", full.names = TRUE, pattern = "*_modpileup_5mC.bed", recursive = TRUE)
# opt$input_files <- test_files
#get file names
fileNames <- gsub(".bed","",basename(opt$input_files))

#read all files
lapply(opt$input_files, function(x) fread(x, sep = "\t", select=c("V1", "V3", "V6", "V10"))) -> dt_list

#rename files and remove non standard chroms
names(dt_list) <- fileNames
dt_list_named <- lapply(dt_list, function(x) {
      setnames(x, old = names(x), new = c("chrom", "end", "strand", "Nvalid_cov"))
      x <- x[!chrom %like% "_"][!chrom %like% "m"]
  return(x)  # Return the modified data.table
    })

#merge all files into 1 data table
dt_rbindlist <- rbindlist(dt_list_named, idcol="sample") #merge list of data.tables into one data.table
dt_wider <- dcast(copy(dt_rbindlist), formula = chrom + end + strand ~ sample, value.var = "Nvalid_cov")

#get dt commands 
# # dt_rbindlist[, sample := gsub(".data.*","",file)] #remove .data.* from filenames
# dt_small <- dt_rbindlist[sample(.N,100)]
# dt_small
# lazy_dt1 <- lazy_dt(dt_small)
# lazy_dt1
# lazy_dt1_wide <- pivot_wider(lazy_dt1, id_cols = c("chrom", "end", "strand"), names_from = "sample", values_from = "Nvalid_cov")
# names(as.data.frame(lazy_dt1_wide))

##############
#set values for filtering
max_val=15
max_val_low=5
## percentage of samples needed that to support specific cpg sites as high quality
percThreshold_hiQ_samples <- .75 #
# how many samples needed to support cpg sites as high quality
# minReads_threshold <- 5 #set minimum reads supporting each cpg site
numb_nonSamples_cols <- 3 #number of non sample columns
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
# dtF <- fread(cmd=paste("awk -F'\t' '{print $1, $3, $14}' ",opt$input_file)) #not efficient
# dtF <- fread(cmd=paste("head ",opt$input_file)) #memory efficeint 
# dtF2 <- fread(opt$input_file, select=c("V1", "V3", "V6", "V10", "V11")) #more efficient

# fread("awk -F'\t' '{print $1, $3, $14}' your_file.txt")

# how many samples needed to support cpg sites as high quality
nThreshold_hiQ_samples <- round(opt$percent_Samples_consensus*(ncol(dt_wider)-numb_nonSamples_cols))
dt_wider[, num_samples_with_minReads := rowSums(.SD >= opt$nReads_consensus, na.rm=T), .SDcols = !c("chrom", "end", "strand")]
# Filter CpGs (rows) where 75% of samples have at least 10 (min_reads) reads
consensus_hiQ_CpGs <- dt_wider[num_samples_with_minReads >= nThreshold_hiQ_samples,]
# Remove the num_samples_with_minReads column, no longer need it
consensus_hiQ_CpGs[, num_samples_with_minReads := NULL]


#generate stats for consensus high_quality reads
consensus_hiQ_CpGs[,`:=`(MIN = rowMins(as.matrix(.SD), na.rm=T),SUM=rowSums(.SD, na.rm=T) ,AVG = rowMeans(.SD, na.rm = T)), .SDcols = !c("chrom", "end", "strand")]


#save as bed file
fwrite(consensus_hiQ_CpGs, file = file.path(opt$output_dir,"consensus_HighQual_motif_sites.tsv"), sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)

#save as bed file
fwrite(consensus_hiQ_CpGs[,.(chrom, start=end-1, end, strand)], file =file.path(opt$output_dir, "consensus_HighQual_motif_sites.bed") , sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = FALSE)