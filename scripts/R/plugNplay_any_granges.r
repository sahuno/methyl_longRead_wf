
library(plyranges)
library(tidyverse)
library(optparse)
library(data.table)

option_list <- list(make_option(c("--bedfile"), default = "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/chr8_TE.bed"),
#                    make_option(c("--statsFile"), default = "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/per_read_aggregate/Spectrum-OV-009_N/chr8/data/Spectrum-OV-009_N.chr8.data_aggregate_stats.txt"),
make_option(("--stats2Use"), type="character", default="median",help ="use mean or median aggreagate prob per site"),
                    make_option(c("--statsFile"), default = "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/gather_files/Spectrum-OV-009_T/Spectrum-OV-009_T.data_aggregate_stats_all_chroms.txt"),
                    make_option(("--nReads"), type="integer", default=1,help ="threshold for number minimum numvber reads per site"),
                    make_option("--plots_pdf", type="character"),
                    make_option("--extend_by", type="integer", default = 4000),
make_option("--plots_rds", type="character", default=NULL, help="path to raw plot"))

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)


#make granges from input
make_gr <- function(bedFile){
     dt <- fread(bedFile)
     if(ncol(dt) == 3){
             names(dt) <- c("seqnames", "start", "end")
             dt[,`:=`(strand = "*")]
             dt[,`:=`(id = paste0(seqnames,":" ,start,"-",end), Start_Position= start)]
             dt[,`:=`(id = paste0(seqnames,":" ,start,"-",end,strand,id), Start_Position= start)]
            dt[,`:=`(Start_Position = fcase(strand == "+", end, strand == "-", start, rep(TRUE, .N), start))] #get values for start sites
             #print(dt)
             gr <- makeGRangesFromDataFrame(dt, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
     }
return(gr)
}






#read input bed file
gr <- make_gr(bedFile=opt$bedfile)
gr_ls <- split(gr, gr$id) #make list of ids in granges
# extend_by <- 4000 #to get promoter regions
gr_ls_extended <- lapply(gr_ls, function(x) stretch(anchor_center(x), extend=opt$extend_by)) #add extend_by/2 basePairs each to both ends
#gr_ls_extended <- lapply(gr_ls_extended, function(x) unanchor(x)) #not nvessary cos no more anchored




#########
#read input methylation file
DT <- fread(opt$statsFile)
names(DT) <- c("chrm" ,"pos" ,"strand","number_reads", "mean_Prob_Meth_perSite", "median_Prob_Meth_perSite") #assume there's no header on concatenated files
DT[,`:=`(chr_pos = paste0(chrm,"_",pos), Position=pos)] # make chrom_position index for plotting
DT <- setnames(DT, "chrm", "seqnames")

#select mean or median methylation as stats to use
if(opt$stats2Use == "median"){
    message("using ",opt$stats2Use, " methylation percentage")
    setnames(DT, paste0(opt$stats2Use,"_Prob_Meth_perSite"), "statistic")
}else{
    stop("please specify mean or median for stats2Use")
}

DT_gr <- makeGRangesFromDataFrame(DT, start.field = "pos", end.field = "pos", keep.extra.columns=TRUE)


##############
#join overlaps bwteen grangesList and methylation data
overlaps_methyl_grangesList <- lapply(gr_ls_extended, function(x) join_overlap_left(x, DT_gr))
overlaps_methyl_grangesList_DT <- lapply(overlaps_methyl_grangesList, function(x) x %>% as.data.table()) #convert to data table
names(overlaps_methyl_grangesList_DT)

#remove nameless TE
keep = (!names(overlaps_methyl_grangesList_DT) == "")
overlaps_methyl_grangesList_DT <- overlaps_methyl_grangesList_DT[keep]
names(overlaps_methyl_grangesList_DT)
#save plot data- prefiltering
#saveRDS(overlaps_TE_methyl_DTList, file = paste0(opt$output_rds))

#filter low numb reads
overlaps_methyl_grangesList_DT <- lapply(overlaps_methyl_grangesList_DT, function(x) x[(number_reads >= opt$nReads),])


##############################
###plots
plot_TE_plus_promoters <- function(DT, TE_id){

#plot if a valid number to median_PrMethylated exists
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
                            axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) #+
                            # geom_vline(aes(xintercept = start), colour="blue",start_cord)

   return(plt_TE)
   }
}

#plot and store as list
plots_list <- imap(overlaps_methyl_grangesList_DT, ~plot_TE_plus_promoters(DT=.x , TE_id=.y))

#save as rds on disk
saveRDS(plots_list, file = paste0(opt$plots_rds))

#Export as pdf
#pdf(file=paste0("test_methylation_plugNplay.pdf"), height = 12, width = 30)
pdf(file=paste0(opt$plots_pdf),height = 12, width = 30)
print(plots_list)
dev.off()