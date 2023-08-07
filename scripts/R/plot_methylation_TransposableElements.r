
#ahuno samuel
##purpose  - line1 methylation profile
library(data.table)
library(optparse)
library(plyranges)
library(dtplyr)


#make option list
option_list <- list(
make_option(("--chrom"), default="18"),
make_option(("--output_rds"), type="character", default = "/juno/work/greenbaum/users/ahunos/rotation/scripts/workflows/qc_DMR/results/per_read_aggregate/BRCA_13135_P_2/chr18/data/BRCA_13135_P_2.chr18.LINE_methylation.rds"),
make_option(("--nReads"), type="integer", default=1,help ="threshold for number minimum numvber reads per site"),
make_option(("--stats2Use"), type="character", default="median",help ="use mean or median aggreagate prob per site"),
make_option(("--input_file"), type="character", default = "/juno/work/greenbaum/users/ahunos/rotation/scripts/workflows/qc_DMR/results/per_read_aggregate/BRCA_13135_P_2/chr18/data/BRCA_13135_P_2.chr18.data_aggregate_stats.txt"),
make_option(("--TE_ref"), type="character", default = "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/hsflil1_8438.bed"),
make_option("--plots_pdf", type="character"),
make_option("--plots_rds", type="character", default=NULL, help="path to raw plot")
)

opt <- parse_args(OptionParser(option_list = option_list))


#read TE reference
#path <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/hsflil1_8438.bed"
TE_dt <- fread(opt$TE_ref, skip=1)
names(TE_dt) <- c("seqnames", "start", "end", "id", "score", "strand", "pos1", "pos2", "colorID")

#filter only chrom of interest
TE_dt <- TE_dt[seqnames == paste0("chr",opt$chrom), .SD, .SDcols=c("seqnames", "start", "end","strand", "id")]
#TE_dt <- TE_dt[seqnames == paste0("chr",opt$chrom), .SD(), .SDcols=c("seqnames", "start", "end","c", "id")]
TE_dt[,`:=`(id = paste0(seqnames,":" ,start,"-",end,strand,id), Start_Position= start)]
TE_dt[,`:=`(Start_Position = fcase(strand == "+", end, strand == "-", start, rep(TRUE, .N), start))] #get values for start sites
# ldt <- lazy_dt(TE_dt)
# ldt %>% mutate(Start_Position = case_when(strand == "+" ~ end, 
#                                         strand == "-" ~ start,
#                                         TRUE ~ start))


TE_ls <- split(TE_dt, TE_dt$id) #split TE by ids
TE_grlist <- lapply(TE_ls, function(x) makeGRangesFromDataFrame(x, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE))
TE_grlist <- lapply(TE_grlist, function(x) x %>% mutate(Start_Position = start)) #attach a start position for ploting



#TE_grlist[[1]] %>% mutate(Start_Position = start)
#promoters(TE_grlist[[1]], upstream = 1000, downstream=1000)
##
#set promoter regions of Transposable elements
extend_by <- 1000 #to get promoter regions
TE_promoters_grlist <- lapply(TE_grlist, function(x) stretch(anchor_5p(x), extend=extend_by))
#TE_promoters_grlist <- lapply(TE_promoters_grlist, function(x) tile_ranges(x, width=1)) #tile granges

#length(TE_promoters_grlist)

#stretch(anchor_5p(gr_test), extend=10) #this is what we need for line promoter
#convert to granges
#TE_gr <- makeGRangesFromDataFrame(TE_dt, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)



###get methylation data
##load data
#DT_cp <- fread("/juno/work/greenbaum/users/ahunos/rotation/scripts/workflows/qc_DMR/results/per_read_aggregate/BRCA_13135_P_2/chr18/data/BRCA_13135_P_2.chr18.data_aggregate_stats.txt")
DT_cp <- fread(opt$input_file)
DT_cp[,`:=`(chr_pos = paste0(chrm,"_",pos), Position=pos)] # make chrom_position index for plotting
DT_cp <- setnames(DT_cp, "chrm", "seqnames")

#select mean or median methylation as stats to use
if(opt$stats2Use == "median"){
    message("using ",opt$stats2Use, " methylation percentage")
    setnames(DT_cp, paste0(opt$stats2Use,"_Prob_Meth_perSite"), "statistic")
}else{
    stop("please specify mean or median for stats2Use")
}

#if assume methylation of low coverage regions, if less than minimum number of reads assign meethylated
#DT_cp[,`:=`(statistic = fcase(number_reads <= 3, 0, default = x))] #assume methylated in regions with low read counts
DT_subj <- makeGRangesFromDataFrame(DT_cp, start.field = "pos", end.field = "pos", keep.extra.columns=TRUE)


##############
##
#join overlaps bwteen TE and methylation data
overlaps_TE_methyl_grList <- lapply(TE_promoters_grlist, function(x) join_overlap_left(x, DT_subj))
overlaps_TE_methyl_DTList <- lapply(overlaps_TE_methyl_grList, function(x) x %>% as.data.table()) #convert to data table
# dput(head(overlaps_TE_methyl_DTList[[1]]))
# names(overlaps_TE_methyl_DTList)

#lapply(overlaps_TE_methyl_DTList, dim)
############
#remove nameless TE
keep = (!names(overlaps_TE_methyl_DTList) == "")
overlaps_TE_methyl_DTList <- overlaps_TE_methyl_DTList[keep]

#save plot data- prefiltering
saveRDS(overlaps_TE_methyl_DTList, file = paste0(opt$output_rds))

overlaps_TE_methyl_DTList <- lapply(overlaps_TE_methyl_DTList, function(x) x[(number_reads >= opt$nReads),])



##############################
###plots
plot_TE_plus_promoters <- function(DT, TE_id){

#plot if a valid number to statistic exists
   if(any(!is.na(DT$statistic))){

        DT <- DT[order(seqnames, Position)]
        start_cord <- copy(DT)[,.(start = unique(Start_Position))]
# print(start_cord)
        plt_TE <- ggplot(DT, aes(x=Position, y=statistic)) + 
                            geom_point() + 
                            geom_smooth(color="black", span = 0.5) +
                            coord_cartesian(ylim = c(0,1)) + 
                            labs(title = paste0(TE_id))+
                            theme(axis.ticks.x=element_blank(), 
                            axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
                            geom_vline(aes(xintercept = start), colour="blue",start_cord)

   return(plt_TE)
   }
}

#plot and store as list
plots_list <- imap(overlaps_TE_methyl_DTList, ~plot_TE_plus_promoters(DT=.x , TE_id=.y))

#save as rds on disk
saveRDS(plots_list, file = paste0(opt$plots_rds))

#Export as pdf
pdf(file=paste0(opt$plots_pdf),height = 12, width = 30)
#pdf(file=paste0("test_TE_methylation.pdf"),height = 9, width = 24)
print(plots_list)
dev.off()
