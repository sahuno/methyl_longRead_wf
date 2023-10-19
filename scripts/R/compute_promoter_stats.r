#get consensus high qulaity cpgs for proportion of disconcordant reads analysis,
#computes stats for promoter regions and plots the results
#name: samuel ahuno
#date: sept 1st 2023 

#load library
library(tidyverse)
library(data.table)
library(matrixStats)
library(optparse)
library(dtplyr)
library(cowplot)
library(parallel)
#set
# Detect the number of cores
num_cores <- detectCores()



option_list <- list(
   make_option(("--input_file"), type = "character", default="/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/results/promoters/D-A-3/data/D-A-3.chr1.filtered_Consensus.overlaps.bed"),
    make_option(("--nReads"), type="integer", default=5,help ="threshold for number minimum numvber reads per site"),
    # make_option(("--stats2Use"), type="character", default="median",help ="use mean or median aggreagate prob per site"),
    make_option(("--methyl_percent"), type="numeric", default=70,help ="threshold for site to be called methylated, used in combination with nReads"),
    make_option(("--methyl_metrics_promoter_data"), type="character", default="methyl_metrics_promoter_data.tsv",help ="methylation activity per gene promoter data in .txt")
    # make_option(("--overlaps_rds"), type="character", default="overlaps.rds",help ="data for overlaps"),
  #  make_option("--plots_pdf", type="character", default="plots.pdf", help="plot name"),
  #  make_option("--plots_rds", type="character", default="plot_object.rds", help="name of plot file")
   )

opt <- parse_args(OptionParser(option_list = option_list))


message("done reading in bed files as data.table")
dt_ov <- fread(opt$input_file)
# dt_ov[, .N, by = "V19"] #check number of genes we found overlaps for
# dt_ov[, .(MIN=min(V6)), by = "V19"][,.(min(MIN))] #sanity check min number of reads shouldn't be less than number set for consensus

#c("chrom", "start", "end", "mod_base", "score", "strand", "start2", "end2", "color", "Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")
setnames(dt_ov, c("V1",  "V2",  "V3",  "V4",  "V5",  "V6",  "V7",  "V8",  "V9",  "V10", "V11", "V12","V13", "V14" ,"V15", "V16", "V17", "V18", 
"V19", "V20", "V21", "V22", "V23", 
"V24", "V25", "V26", "V27", "V28", "V29"
), 
c("chrom", "start", "end", "mod_base", "score", "strand", "start2", "end2", "color", "Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall",
"prom_chrom", "prom_start", "prom_end", "prom_strand", "key",
"ensgene.version", "symbol", "biotype","transcription_start_site", "numGCs", "numCGs"))


dt_ov_ls  <- split(dt_ov, dt_ov$key) #split by keys

# dt_ov_ls[[1]]

##############################################################
message("computing methylation statistics per gene promoter in chr", opt$chrom)



#Revised entropy
compute_entropy_Avg <- function(Px){
  OneMinusPx =  1-Px #cal prob of not methylated per spg site
  entropy_Px_1minPx <- -(Px * log(Px)) - (OneMinusPx*log(OneMinusPx))
  sum_methyl_entropy <- sum(entropy_Px_1minPx, na.rm = TRUE)
  avg_methyl_entropy <- mean(entropy_Px_1minPx, na.rm = TRUE)
  
  return(avg_methyl_entropy)
}

#function to compute entropy
compute_entropy_shannon <- function(Px){
  shann_methyl_entropy <- -sum((Px * log2(Px)), na.rm = TRUE)  
  return(shann_methyl_entropy)
}


compute_geom_mean <- function(x){
  #compute geometric mean of methylation values
  gm <- exp(mean(log(x[x > 0]),na.rm=TRUE))
  return(gm)
}



# ##############################################################
#Function to compute average methylation per gene promoter
methylation_prom <- function(dat_in, statistic = NULL, cov_col=NULL){
    #compute average methylation per gene promoter; no filtering just compute avaerge of methylation values
    #compute proportion of methylated cps per gene promoter; minimun number of reads per site and minimum methylation percentage
    promoter_methyl_rate <- dat_in[,Meth_or_UnMeth := 
                                    data.table::fcase(get(statistic) >= opt$methyl_percent & get(cov_col) >= opt$nReads,"M" , default = "U")][
                                        ,.(
                                          gene.id=unique(ensgene.version),
                                        median_promoter_methyl = median(get(statistic), na.rm=TRUE), 
                                        geom_mean_promoter_methyl = compute_geom_mean(get(statistic)), 
                                            methyl_promoter_entropy_Avg = compute_entropy_Avg(get(statistic)), 
                                            methyl_promoter_entropy_shann = compute_entropy_shannon(get(statistic)), 
                                            nCpGs_promoter_observed_data = length(get(statistic)),
                                            nGCs_promoter_expected_ref = unique(numGCs) * 2,
                                            nCpGs_promoter_expected_ref = unique(numCGs) * 2,
                                            prop_methyl_in_promoter_data = sum(Meth_or_UnMeth == "M")/length(Meth_or_UnMeth)), by = key]
return(promoter_methyl_rate)
}

# #compute methylation rate per gene promoter
methyl_rate_promoters_ls <- lapply(dt_ov_ls, methylation_prom, 
                              statistic = "fraction_modified", 
                              cov_col="Nvalid_cov")




methyl_rate_promoters_dt <- rbindlist(methyl_rate_promoters_ls, fill=FALSE, idcol=NULL) #bind all data.tables together
setkey(methyl_rate_promoters_dt, key) #set key
# duplicated(methyl_rate_promoters_dt, by="key") #check for duplicates
# methyl_rate_promoters_dt[, .N, by = "key"] #check number of genes we found overlaps for
# methyl_rate_promoters_dt[, `:=`(seqnames = paste0("chr",opt$chrom))]

methyl_rate_promoters_dt[,min(nCpGs_promoter_observed_data)] #check number of genes we found overlaps for
# # methyl_rate_promoters_dt[,fD := .N > 1, by = "key"]

print(head(methyl_rate_promoters_dt))
message("saving methylation rate data as tsv")
fwrite(methyl_rate_promoters_dt, file = paste0(opt$methyl_metrics_promoter_data), sep="\t")
message("done saving methylation rate data as tsv")
 