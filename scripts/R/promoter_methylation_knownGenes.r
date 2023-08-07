# purpose: find promoter methylation of known genes 
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(AnnotationDbi)
# library(org.Hs.eg.db)
library(data.table)
library(plyranges)
library(optparse)
#BiomaRt::biomartCacheClear()
#library(biomaRt)

# specify our desired options in a list
option_list <- list(make_option(("--chrom"), type="character", default="1",help ="chromosome to work on"),
make_option(("--nReads"), type="integer", default=3,help ="threshold for number minimum numvber reads per site"),
make_option(("--data_rds"), type="character", default=NULL,help ="intersected promoter data for plot, pre filtered"),
   make_option("--plots_pdf", type="character", default=NULL, help="plot name"),
   make_option("--cache", type="logical", default=FALSE, help="use biomart cache"),
   make_option("--version", type="character", default="latest", help="which version of human genome to use"),
   make_option("--plots_rds", type="character", default=NULL, help="directory to store"),
   make_option(("--input_file"), type = "character", default="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/per_read_aggregate/Spectrum-OV-009_N/chr1/data/Spectrum-OV-009_N.chr1.data_aggregate_stats.txt")
   )
opt <- parse_args(OptionParser(option_list = option_list))


##load data
#DT_cp <- fread("/juno/work/greenbaum/users/ahunos/rotation/scripts/workflows/qc_DMR/results/per_read_aggregate/BRCA_13135_P_2/chr18/data/BRCA_13135_P_2.chr18.data_aggregate_stats.txt")
DT_cp <- fread(opt$input_file)
DT_cp[,`:=`(chr_pos=paste0(chrm,"_",pos), Position=pos)] # make chrom_position index for plotting
DT_cp <- setnames(DT_cp, "chrm", "seqnames")
DT_subj <- makeGRangesFromDataFrame(DT_cp, start.field = "pos", end.field = "pos", keep.extra.columns=TRUE)

#class(DT_cp)
dim(DT_cp)
######
#get gene promoters according to chromo
source("/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/get_promoters_and_TSS_ensembl_genes.r")
#gene_promoters <- make_promoters(fetch_transcripts(chr = opt$chrom, useCache_option=opt$cache), up_dis = 1000, dwn_dis = 1000)
gene_promoters <- txdb_promoters(species="Hs", version = opt$version, chrom=opt$chrom, up_dis = 1000, dwn_dis = 1000)
#gene_promoters <- txdb_promoters(species="Hs", version = "latest", chrom=1, up_dis = 1000, dwn_dis = 1000)
GenomeInfoDb::seqlevelsStyle(gene_promoters) <- 'UCSC' #change seqlevelStyle to UCSC


#join overlaps bwteen transcripts and methylation data
ov_methyl_prom <- lapply(gene_promoters, function(x) join_overlap_left(x, DT_subj))
ov_methyl_prom_dt <- lapply(ov_methyl_prom, function(x) x %>% as.data.table())
#lapply(ov_methyl_prom_dt, dim)


############
#remove nameless genes
keep = (!names(ov_methyl_prom_dt) == "")
ov_methyl_prom_dt <- ov_methyl_prom_dt[keep]

#save plot data- prefiltering
saveRDS(ov_methyl_prom_dt, file = paste0(opt$data_rds))
message("saving Data as RDS")

ov_methyl_prom_filtered_dt <- lapply(ov_methyl_prom_dt, function(x) x[(number_reads >= opt$nReads),])
# ov_methyl_prom_filtered_dt_ls <- lapply(ov_methyl_prom_filtered_dt, function(x) split(x, x$key))
#lapply(ov_methyl_prom_filtered_dt, dim)



#plot transcipt methylation
##notes about plot
#plots, by positions
## doesn't discrimate agiants posuitive and negative strands
plot_promoters <- function(DT, geneName){

#plot if gene has multiple transcripts and there is a valid number to plot
   if(length(unique(DT$key)) > 1 & any(!is.na(DT$median_PrUnM))){

      DT <- DT[order(seqnames, Position)]
      tss_cord <- copy(DT)[,.(tss = unique(transcription_start_site),start=unique(start),end=unique(end)),by=key]
      dummy_data <- tss_cord[,-"tss"] %>% pivot_longer(!c(key), values_to="Position") %>% mutate(median_PrUnM = case_when(name == "start" ~ 0, name == "end" ~ 1)) %>% dplyr::select(-name)

    plt_promoter_methyl_facet <- ggplot(DT, aes(x=Position, y=median_PrUnM)) + 
                            geom_point() + 
                            geom_smooth(color="black", span = 0.5) +
                            theme(axis.ticks.x=element_blank(), 
                            axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)
                            ) + 
                            facet_wrap(~key, scales = "free_x") + 
                            geom_blank(data=dummy_data) + 
                             coord_cartesian(ylim = c(0,1)) + 
                            geom_vline(aes(xintercept = tss), colour="blue", tss_cord)

   return(plt_promoter_methyl_facet)
   } else if (length(unique(DT$key)) == 1 & any(!is.na(DT$median_PrUnM))) {
            DT <- DT[order(seqnames, Position)]
      tss_cord2 <- copy(DT)[,.(tss = unique(transcription_start_site),start=unique(start),end=unique(end)),by=key]
      plt_promoter_methyl <- ggplot(DT, aes(x=Position, y=median_PrUnM)) + 
                            geom_point() + 
                            geom_smooth(color="black", span = 0.5) +
                            theme(axis.ticks.x=element_blank(), 
                            axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)
                            ) +
                            coord_cartesian(ylim = c(0,1)) + 
                           labs(title = paste0(geneName))+
                            geom_vline(aes(xintercept = tss), colour="blue",tss_cord2)
      return(plt_promoter_methyl)
   } 
}

#plot and store as list
plots_list <- imap(ov_methyl_prom_filtered_dt, ~plot_promoters(DT=.x , geneName=.y))

#save as rds on disk
saveRDS(plots_list, file = paste0(opt$plots_rds))
message("saving list of plots as RDS")



#Export as pdf
#pdf(file="/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/gene_promoters/Spectrum-OV-009_N/chr1/Spectrum-OV-009_N.chr1.methylation_gene_promoters.pdf",height = 12, width = 30)
pdf(file=paste0(opt$plots_pdf),height = 12, width = 30)
print(plots_list)
dev.off()

message("done ploting gene promoters")

