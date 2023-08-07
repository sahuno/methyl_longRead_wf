#name - samuel ahuno
#purpose - find methylation activity per gene promoter

library(optparse)
library(data.table)
library(plyranges)

option_list <- list(
    make_option(("--chrom"), type="character", default="6",help ="chromosome to work on"),
   make_option(("--input_file"), type = "character", default="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/per_read_aggregate/BRCA_13135_P_1/chr6/data/BRCA_13135_P_1.chr6.data_aggregate_stats.txt"),
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
DT_cp <- fread(opt$input_file)
DT_cp[,`:=`(chr_pos=paste0(chrm,"_",pos), Position=pos)] # make chrom_position index for plotting
DT_cp <- setnames(DT_cp, "chrm", "seqnames")

#select mean or median methylation as stats to use
if(opt$stats2Use == "median"){
    message("using ",opt$stats2Use, " methylation percentage")
    setnames(DT_cp, paste0(opt$stats2Use,"_Prob_Meth_perSite"), "statistic")
}else{
    stop("please specify mean or median for stats2Use")
}

#make GRanges object of aggregate methylation data
DT_subj <- makeGRangesFromDataFrame(DT_cp, start.field = "pos", end.field = "pos", keep.extra.columns=TRUE)
####################################


######################main promoter scripts ################
##############################################################
#read in gene promoters and split
load("/juno/work/greenbaum/users/ahunos/apps/methylONT/data/gene_promoters_encode1kb.rda")


#filter only for chr of interest
if(!is.null(opt$chrom)){
    message("limiting further anlysis to chromosome ", opt$chrom)
gene_promoters_encode1kb_chr_subset <- gene_promoters_encode1kb %>% plyranges::filter(seqnames == paste0("chr", opt$chrom))
} else {
    gene_promoters_encode1kb_chr_subset <- gene_promoters_encode1kb
}

#GRangesList(gene_promoters_encode1kb_chr_subset)
#convert to list where each element in the list is genes promoter to make it easier to work with
gene_promoters_encode1kb_ls <- split(gene_promoters_encode1kb_chr_subset, gene_promoters_encode1kb_chr_subset$key) #split by seqnames
#keep = (!names(ov_methyl_prom_dt) == "")
sanity_check_length <- gene_promoters_encode1kb_ls[lapply(gene_promoters_encode1kb_ls, function(x) length(x)) > 1]

if(length(sanity_check_length)==0){
    message("all promoters are unique")
} else {
    message("possibible duplicate promoters persent in query file")
    message("\n",sanity_check_length,"\n")
}

#find length of promoters
# gene_promoters_encode1kb_ls[[1]]
# length(gene_promoters_encode1kb_ls[[1000]])

message("finding overlaps between promoters and methylation data")
#join overlaps bwteen transcripts and methylation data
ov_methyl_prom <- lapply(gene_promoters_encode1kb_ls, function(x) plyranges::join_overlap_left(x, DT_subj))
#ov_methyl_prom <- lapply(gene_promoters_encode1kb_ls[grep(c("KRAS|BRCA"),gene_promoters_encode1kb_ls)], function(x) plyranges::join_overlap_left(x, gr_test))
#ov_methyl_prom[["OR2J1_ENSG00000204702.5"]]


#function to check if data.table is empty
# check_non_empty_dt <- function(x){
#     d <- dim(x)
#     if(!all(dim(ov_methyl_prom_dt[[1]]) == 0)){
#         return(d)
#         }
# }
###sanity checks; make sure all promoters are present
# lapply(ov_methyl_prom, )

message("convert overlapped granges to data.tables")
ov_methyl_prom_dt <- lapply(ov_methyl_prom, function(x) x %>% as.data.table())
message("removing empty data.tables")
ov_methyl_prom_dt <- ov_methyl_prom_dt[sapply(ov_methyl_prom_dt, function(x) dim(x)[1]) > 0]


#save plot data- prefiltering
message("saving overlaps Data as RDS")
saveRDS(ov_methyl_prom_dt, file = paste0(opt$overlaps_rds))
# lapply(ov_methyl_prom_dt, check_non_empty_dt)

# all(dim(ov_methyl_prom_dt[[1]]))
#dt_test2[,Meth_or_UnMeth := data.table::fcase(methylation_percentage_5mC_5hmC >= methyl_percentage_threshold & number_reads >= methyl_numb_reads_threshold,"M" , default = "U")]

#test data
# dtest <- copy(ov_methyl_prom_dt[["OR1I1_ENSG00000094661.3"]])
# dtest2 <- dtest[,Meth_or_UnMeth := data.table::fcase(median_PrUnM >= opt$methyl_percent & number_reads >= opt$nReads,"M" , default = "U")]

#[,`:=`(mean_prom_methyl = median(median_PrUnM), prop_methyl = sum(Meth_or_UnMeth == "M")/length(Meth_or_UnMeth))]






###################
#plot methylations
#function to check if data.table is empty
    uniq_keys <- function(x){
        length(unique(x[, key]))
    }
#uniq_keys(DT, "key")
#uniq_keys(DT[[1]])
#functiom to check if there are valid stats to plot
    valid_stats <- function(x){
        #any(x[, !is.na(..target_cols)])
        yy <- any(x[, !is.na(statistic)])
        return(yy)
        #x[, ..target_cols]
    }
#valid_stats(DT[[1]], target_cols="median_PrUnM")
# any(!is.na(DT[[1]]$median_PrUnM))
# DT[[1]][, (!is.na(median_PrUnM))]

#convert 1 col as stats to use
# DT <- copy(ov_methyl_prom_dt[[1]]) #uncomment out for testing

plot_promoters <- function(DT, featureName){
#plot if gene has multiple transcripts and there is a valid number to plot
   if(uniq_keys(DT) > 1 & valid_stats(DT)){
#option 1: plot all isoforms of same gene in single plot
      DT <- DT[order(seqnames, Position)] #order by position
      tss_cord <- copy(DT)[,.(tss = unique(transcription_start_site),start=unique(start),end=unique(end)),by=key]

      plot_addOn_data <- tss_cord[,-"tss"] %>% pivot_longer(!c(key), values_to="Position") %>% mutate(statistic = case_when(name == "start" ~ 0, name == "end" ~ 1)) %>% dplyr::select(-key)

    plt_promoter_methyl_facet <- ggplot(DT, aes(x=Position, y=statistic)) + 
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
      plt_promoter_methyl <- ggplot(DT, aes(x=Position, y=statistic)) + 
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





#ov_methyl_prom_dt[["OR2J1_ENSG00000204702.5"]]
##############################################################
##############################################################
message("computing methylation statistics per gene promoter in chr", opt$chrom)
#function to compute entropy
compute_entropy <- function(Px){
methyl_entropy <- -sum(Px * log2(Px))
return(methyl_entropy)
}

##############################################################
#Function to compute average methylation per gene promoter
methylation_prom <- function(dat_in){
    #compute average methylation per gene promoter; no filtering just compute avaerge of methylation values
    #compute proportion of methylated cps per gene promoter; minimun number of reads per site and minimum methylation percentage
    promoter_methyl_rate <- dat_in[,Meth_or_UnMeth := 
                                    data.table::fcase(statistic >= opt$methyl_percent & number_reads >= opt$nReads,"M" , default = "U")][
                                        ,.(median_promoter_methyl = median(statistic, na.rm=TRUE), 
                                            methyl_promoter_entropy = compute_entropy(statistic), 
                                            nCpGs_promoter_observed = length(statistic),
                                            nGCs_promoter_expected = unique(numGCs) * 2,
                                            nCpGs_promoter_expected = unique(numCGs) * 2,
                                            proportion_methyl_promoter = sum(Meth_or_UnMeth == "M")/length(Meth_or_UnMeth)), by = key]
return(promoter_methyl_rate)
}

# tb <- table(names(ov_methyl_prom_dt))
# which(tb == 1)

#compute methylation rate per gene promoter
methyl_rate_promoters_ls <- lapply(ov_methyl_prom_dt, methylation_prom)
#methyl_rate_promoters_ls[["OR2J1_ENSG00000204702.5"]]
##to do: add gc content to promoter data.table

#gsub(".*_","","OR2J1_ENSG00000204702.5")
#gsub("\\..*","","OR2J1_ENSG00000204702.5")

methyl_rate_promoters_dt <- rbindlist(methyl_rate_promoters_ls, fill=FALSE, idcol=NULL) #bind all data.tables together
setkey(methyl_rate_promoters_dt, key) #set key
# duplicated(methyl_rate_promoters_dt, by="key") #check for duplicates

methyl_rate_promoters_dt[, `:=`(seqnames = paste0("chr",opt$chrom))]

# methyl_rate_promoters_dt[,fD := .N > 1, by = "key"]
message("saving methylation rate data as txt")
fwrite(methyl_rate_promoters_dt, file = paste0(opt$methyl_metrics_promo_data), sep="\t")



##############################################################
# ##############################################################
# ### To do: merge with expression data in different script
# expr_data <- "/juno/work/greenbaum/users/ahunos/rotation/data/dds_normalized_cnts_parp_inhibitor.tsv"
# expr_dt <- fread(expr_data)

# # expr_dt[ensgene %like% "ENSG00000178605"] #this supports the fact that it makes sense to merge on `ensgene.version` since dseq object had `ENSG00000178605.13` & `ENSG00000178605.13_PAR_Y` in it's count matrix 
# # expr_dt[,`:=`(ensgene2 = gsub("\\..*","",ensgene))][,`:=`(num2 = .N),by=ensgene2][num2>1,]

# #methyl_rate_promoters_dt[key =="DUSP22_ENSG00000112679.14",]
# methyl_rate_promoters_exprs_dt <- methyl_rate_promoters_dt[,`:=`(ensgene = gsub(".*_","",key))][expr_dt, on = "ensgene", nomatch = NULL]
# # methyl_rate_promoters_exprs_dt <- methyl_rate_promoters_dt[,`:=`(ensgene = gsub(".*_","",gsub("\\..*","",key)))][expr_dt[,`:=`(ensgene = gsub("\\..*","",ensgene))], on = "ensgene", nomatch = NULL] #don't use, cos no specific mapping of exp and methylation data
# # sum(!is.na(methyl_rate_promoters_exprs_dt$key))
# setkey(methyl_rate_promoters_exprs_dt, key)
# # install.packages("GGally")
# #duplicated(methyl_rate_promoters_exprs_dt, by="key")

# library(GGally)
# plt_ggpairs <- ggpairs(methyl_rate_promoters_exprs_dt[,.(RNA_expres = P_1, med_methyl=median_promoter_methyl, entropy_methyl = methyl_promoter_entropy, nCpGs_obs = nCpGs_promoter_observed, nGCs_exp=nGCs_promoter_expected, nCpGs_exp=nCpGs_promoter_expected, met_prop=proportion_methyl_promoter)], 
#                        lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)), 
#                        diag = list(continuous = wrap("barDiag", binwidth = 0.1, alpha = 0.3, size = 0.5)), 
#                        title = paste0("Gene Promoter Methylation activity per Chrm", opt$chrom)) + 
#                        theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#                        theme(panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# ggsave(plt_ggpairs, file="test_methylation_rate_promoters.pdf", width = 12, height=7)









#test
# gr <- GRanges(seqnames = "chr6", strand = c("*", "*", "*"),
#               ranges = IRanges(start = c(1,3,5), width = 1000))
# values(gr) <- DataFrame(score = c(0.1, 0.5, 0.3))
# gr
#grt <- start(Biostrings::matchPattern("CG", BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]]))
#Biostrings::countPattern("CG", BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]])
# BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]]

##questions: do number of cpgs correlate with promoter methylation?
##questions: do number of cpgs correlate with promoter methylation entropy?
##questions: do number of cpgs correlate with proportion of methylated cpgs?
##questions: do proportion of methylated cpgs correlate with gene expression?




# library(BSgenome.Hsapiens.UCSC.hg38)
# lapply(gr, function(x) Biostrings::countPattern("CG", BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]]))
# BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[gr[1]]]
# library(BSgenome.Hsapiens.UCSC.hg38)
# gr_pro <- dput(head(gene_promoters_encode1kb_chr_subset[,c("symbol", "key")], 3))
# lapply(gr_pro, function(x) Biostrings::countPattern("CG", BSgenome.Hsapiens.UCSC.hg38::Hsapiens[x]))
# #this is the same as above
# ov_gr_prom <- getSeq(Hsapiens, gr_pro)
# Biostrings::vcountPattern("GC", ov_gr_prom)

# BSgenome.Hsapiens.UCSC.hg38::Hsapiens[gr_pro]
