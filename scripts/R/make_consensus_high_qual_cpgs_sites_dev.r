#get consensus high qulaity cpgs for proportion of disconcordant reads analysis
#name: samuel ahuno
#date: sept 1st 2023 

#load library
library(tidyverse)
library(data.table)
library(matrixStats)

#set vals
max_val=15
max_val_low=5
## percentage of samples needed that to support specific cpg sites as high quality
percThreshold_hiQ_samples <- .75 #
# how many samples needed to support cpg sites as high quality
minReads_threshold <- 5 #set minimum reads supporting each cpg site

#create toy data
dt_s <- data.table(chrm = paste0("chrm",1), 
                   pos = floor(runif(10, min=10000, max=50000)),
                   sample1 = floor(runif(10, min=0, max=max_val_low)),
                   sample2 = floor(runif(10, min=0, max=max_val)), 
                   sample3 = floor(runif(10, min=0, max=max_val)), 
                   sample4 = floor(runif(10, min=0, max=max_val)))
dt_s
# dput(dt_s)


# how many samples needed to support cpg sites as high quality
nThreshold_hiQ_samples <- round(percThreshold_hiQ_samples*(ncol(dt_s)-2))
dt_s[, num_samples_with_minReads := rowSums(.SD >= minReads_threshold, na.rm=T),.SDcols = !c("chrm", "pos")]
# Filter CpGs (rows) where 75% of samples have at least 10 (min_reads) reads
consensus_hiQ_CpGs <- dt_s[num_samples_with_minReads >= nThreshold_hiQ_samples,]
# Remove the num_samples_with_minReads column, no longer need it
consensus_hiQ_CpGs[, num_samples_with_minReads := NULL]


#generate stats for consensus high_quality reads
consensus_hiQ_CpGs[,`:=`(MIN = rowMins(as.matrix(.SD), na.rm=T),SUM=rowSums(.SD, na.rm=T) ,AVG = rowMeans(.SD, na.rm = T)), .SDcols = !c("chrm", "pos")]

consensus_hiQ_CpGs