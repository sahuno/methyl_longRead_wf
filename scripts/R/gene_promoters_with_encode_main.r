#name - samuel ahuno
#purpose - find methylation activity per gene promoter

library(optparse)
library(data.table)
library(plyranges)

option_list <- list(
    make_option(("--chrom"), type="character", default="6",help ="chromosome to work on"),
   make_option(("--input_file"), type = "character", default="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/per_read_aggregate/BRCA_13135_P_1/chr6/data/BRCA_13135_P_1.chr6.data_aggregate_stats.txt"),
   make_option(("--input_promoter"), type = "character", default="/juno/work/greenbaum/users/ahunos/apps/methylONT/data/gene_promoters_encode1kb_proteinCoding.rda"),
    make_option(("--nReads"), type="integer", default=3,help ="threshold for number minimum numvber reads per site"),
    make_option(("--stats2Use"), type="character", default="median",help ="use mean or median aggreagate prob per site"),
    make_option(("--methyl_percent"), type="numeric", default=0.9,help ="threshold for site to be called methylated, used in combination with nReads"),
    make_option(("--methyl_metrics_promo_data"), type="character", default=NULL,help ="methylation activity per gene promoter data in .txt"),
    make_option(("--overlaps_rds"), type="character", default=NULL,help ="data for overlaps"),
   make_option("--plots_pdf", type="character", default=NULL, help="plot name"),
   make_option("--plots_rds", type="character", default=NULL, help="name of plot file")
   )

opt <- parse_args(OptionParser(option_list = option_list))



##################################
#################################
## to delete after testing
# DT_in <- fread(opt$input_file)
# DT_cp[,`:=`(chr_pos=paste0(chrm,"_",pos), Position=pos)] # make chrom_position index for plotting
# DT_cp <- setnames(DT_cp, "chrm", "seqnames")

# #select mean or median methylation as stats to use
# if(opt$stats2Use == "median"){
#     message("using ",opt$stats2Use, " methylation percentage")
#     setnames(DT_cp, paste0(opt$stats2Use,"_Prob_Meth_perSite"), "statistic")
# }else{
#     stop("please specify mean or median for stats2Use")
# }

promoters_genes <- function(path_dt="/Users/samuelahuno/PhD/research/GreenbaumLab/data/D-A-1_modpileup_5mC.bed", extra, path_promoter_file=opt$input_promoter){
#function takes in paths to dt
DT_in <- fread(path_dt)
#names(DT_in) <- c("chrom", "start", "end", "mod_base", "score", "strand", "start2", "end2", "color", "Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")



###condition;; 
# filter; keep only consensus regions from sample
DT_in_filtered <- DT_in[, DT_in_filter, on = c("chrom", "start", "end")]

#specify if zero-based or 1-based bed file
#make GRanges object of aggregate methylation data
DT_subj <- makeGRangesFromDataFrame(DT_in_filtered, start.field = "start", 
                                    end.field = "end", 
                                    keep.extra.columns=TRUE,
                                    starts.in.df.are.0based = TRUE)
####################################


######################main promoter scripts ################
##############################################################
#read in gene promoters and split
load(path_promoter_file)
gene_promoters_encode1kb <- gene_promoters_encode1kb_proteinCoding


#filter only for chr of interest
if(!is.null(opt$chrom)){
    message("limiting further anlysis to chromosome ", opt$chrom)
gene_promoters_encode1kb_chr_subset <- gene_promoters_encode1kb %>% plyranges::filter(seqnames == paste0("chr", opt$chrom))
} else {
    gene_promoters_encode1kb_chr_subset <- gene_promoters_encode1kb
}

#convert to list where each element in the list is genes promoter to make it easier to work with
gene_promoters_encode1kb_ls <- split(gene_promoters_encode1kb_chr_subset, gene_promoters_encode1kb_chr_subset$key) #split by seqnames

sanity_check_length <- gene_promoters_encode1kb_ls[lapply(gene_promoters_encode1kb_ls, function(x) length(x)) > 1]
if(length(sanity_check_length)==0){
    message("all promoters are unique")
} else {
    message("possibible duplicate promoters persent in query file")
    message("\n",sanity_check_length,"\n")
}



message("finding overlaps between promoters and methylation data")
#join overlaps bwteen transcripts and methylation data
ov_methyl_prom <- lapply(gene_promoters_encode1kb_ls, function(x) plyranges::join_overlap_left(x, DT_subj))


message("convert overlapped granges to data.tables")
ov_methyl_prom_dt <- lapply(ov_methyl_prom, function(x) x %>% as.data.table())
message("removing empty data.tables")
ov_methyl_prom_dt <- ov_methyl_prom_dt[sapply(ov_methyl_prom_dt, function(x) dim(x)[1]) > 0]


#save plot data- prefiltering
message("saving overlaps Data as RDS")
saveRDS(ov_methyl_prom_dt, file = paste0(opt$overlaps_rds))


###################
#plot methylations
#function to check if data.table is empty
    uniq_keys <- function(x){
        length(unique(x[, key]))
    }

#functiom to check if there are valid stats to plot
    valid_stats <- function(x){
        yy <- any(x[, !is.na(statistic)])
        return(yy)
    }


#convert 1 col as stats to use
# DT <- copy(ov_methyl_prom_dt[[1]]) #uncomment out for testing

plot_promoters <- function(DT, featureName){
#plot if gene has multiple transcripts and there is a valid number to plot
   if(uniq_keys(DT) > 1 & valid_stats(DT)){
#option 1: plot all isoforms of same gene in single plot
      DT <- DT[order(seqnames, Position)] #order by position
      tss_cord <- copy(DT)[,.(tss = unique(transcription_start_site),start=unique(start),end=unique(end)),by=key]

      plot_addOn_data <- tss_cord[,-"tss"] %>% pivot_longer(!c(key), values_to="Position") %>% mutate(statistic = case_when(name == "start" ~ 0, name == "end" ~ 1)) %>% dplyr::select(-key)

    plt_promoter_methyl_facet <- ggplot(DT, aes(x=Position, y={{statistic}})) + 
                            geom_point() + 
                            geom_smooth(color="black", span = 0.5) +
                            theme(axis.ticks.x=element_blank(), 
                            axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)
                            ) + 
                            facet_wrap(~key, scales = "free_x") + 
                            geom_blank(data=plot_addOn_data) + 
                            coord_cartesian(ylim = c(0,1)) + 
                            geom_vline(aes(xintercept = tss), colour="blue", tss_cord)
#ggsave(paste0(featureName,"_promoter_methyl_facet.png"), plt_promoter_methyl_facet, width = 12, height = 9)
   return(plt_promoter_methyl_facet)
   } else if (uniq_keys(DT) == 1 & valid_stats(DT)) {
    #option 2: assume there's only 1 isoform per gene. 

      DT <- DT[order(seqnames, Position)] #order by position
      tss_cord2 <- copy(DT)[,.(tss = unique(transcription_start_site), start=unique(start), end=unique(end)), by=key]
      plt_promoter_methyl <- ggplot(DT, aes(x=Position, y={{statistic}})) + 
                            geom_point() + 
                            geom_smooth(color="black", span = 0.5) +
                            theme(axis.ticks.x=element_blank(), 
                            axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)
                            ) +
                            coord_cartesian(ylim = c(0,1)) + 
                           labs(title = paste0(featureName))+
                            geom_vline(aes(xintercept = tss), colour="blue", tss_cord2)
    #ggsave(paste0(featureName,"_promoter_methyl.png"), plt_promoter_methyl, width = 12, height = 9)
      return(plt_promoter_methyl)
   } 
}

#plot and store as list
plots_list <- imap(ov_methyl_prom_dt, ~plot_promoters(DT=.x , featureName=.y))
#plots_list <- iwalk(copy(ov_methyl_prom_dt[1]), ~plot_promoters(DT=.x , featureName=.y))

#save as rds on disk
saveRDS(plots_list, file = paste0(opt$plots_rds))
message("saving list of plots as RDS")

#save plots as pdf
pdf(file=paste0(opt$plots_pdf),height = 12, width = 30)
print(plots_list)
dev.off()
######################





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



##############################################################
#Function to compute average methylation per gene promoter
methylation_prom <- function(dat_in){
    #compute average methylation per gene promoter; no filtering just compute avaerge of methylation values
    #compute proportion of methylated cps per gene promoter; minimun number of reads per site and minimum methylation percentage
    promoter_methyl_rate <- dat_in[,Meth_or_UnMeth := 
                                    data.table::fcase(statistic >= opt$methyl_percent & number_reads >= opt$nReads,"M" , default = "U")][
                                        ,.(median_promoter_methyl = median(statistic, na.rm=TRUE), 
                                        geom_mean_promoter_methyl = compute_geom_mean(statistic), 
                                            methyl_promoter_entropy_Avg = compute_entropy_Avg(statistic), 
                                            methyl_promoter_entropy_shann = compute_entropy_shannon(statistic), 
                                            nCpGs_promoter_observed = length(statistic),
                                            nGCs_promoter_expected = unique(numGCs) * 2,
                                            nCpGs_promoter_expected = unique(numCGs) * 2,
                                            proportion_methyl_promoter = sum(Meth_or_UnMeth == "M")/length(Meth_or_UnMeth)), by = key]
return(promoter_methyl_rate)
}

#compute methylation rate per gene promoter
methyl_rate_promoters_ls <- lapply(ov_methyl_prom_dt, methylation_prom)


methyl_rate_promoters_dt <- rbindlist(methyl_rate_promoters_ls, fill=FALSE, idcol=NULL) #bind all data.tables together
setkey(methyl_rate_promoters_dt, key) #set key
# duplicated(methyl_rate_promoters_dt, by="key") #check for duplicates

methyl_rate_promoters_dt[, `:=`(seqnames = paste0("chr",opt$chrom))]

# methyl_rate_promoters_dt[,fD := .N > 1, by = "key"]
message("saving methylation rate data as txt")
fwrite(methyl_rate_promoters_dt, file = paste0(opt$methyl_metrics_promo_data), sep="\t")

}
