
#remove nameless data tables
keep = (!names(dt_ov_ls) == "")
dt_ov_ls <- dt_ov_ls[keep]
dt_ov_ls_cp <- copy(dt_ov_ls)
dt_ov_ls_cp1 <- dt_ov_ls_cp[[1]]
#filter low numb reads
# dt_ov_ls_ <- lapply(dt_ov_ls, function(x) x[(Nvalid_cov >= opt$nReads),])

dt_ov_ls_cp1 <- dt_ov_ls_cp1[order(start)]

dt_ov_ls_cp1_2 <- dt_ov_ls_cp[1:2]
dt_ov_ls_cp1_2_GRlist <- lapply(dt_ov_ls_cp1_2, function(x)
makeGRangesFromDataFrame(x, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
)

dt_ov_ls_cp1_2_GRlist <- GRangesList(dt_ov_ls_cp1_2_GRlist)
overlaps <- findOverlaps(query = dt_ov_ls_cp1_2_GRlist[[1]], subject = out_cps[[1]])
findOverlaps(query = dt_ov_ls_cp1_2_GRlist, subject = out_cps)

class(dt_ov_ls_cp1_2_GRlist)
class(out_cps)
#antijoin cpgs with 
# library(plyranges)
# Null_ov <- dt_ov_ls_cp1[, .(prom_chrom, prom_start,  prom_end, mod_base = "*", score = 0, prom_strand, key, Nvalid_cov=0,fraction_modified=0.0)] %>% head(1) %>% as.data.frame()
# setnames(Null_ov, c("prom_chrom", "prom_start", "prom_end", "prom_strand"), c("seqnames", "start", "end", "strand"))
# rng <- Null_ov %>% as_granges()

# #tile the genome
# bin_1bp <- tile(rng, width=1)


##############################
###plots
    plot_promoters <- function(dat_in, key, statistic = NULL, add_tss=TRUE, cov_col=NULL){
            start_prom <- unique(dat_in$prom_start)
            end_prom <- unique(dat_in$prom_end)


            plt_gene <- ggplot(dat_in, aes(x=start, y={{statistic}})) + 
                                geom_point() + 
                                geom_smooth(color="black", span = 0.5) +
                                scale_x_continuous(breaks = seq(start_prom, end_prom, by = 50)) +
                                coord_cartesian(ylim = c(0,100), xlim = c(start_prom, end_prom)) + 
                                labs(title = paste0(key))+
                                theme(axis.ticks.x=element_blank(), 
                                axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) 
            if(add_tss == TRUE){
                plt_gene <- plt_gene + geom_vline(aes(xintercept = unique(dat_in$transcription_start_site)), colour="blue")
            }

    return(plt_gene)

    }

# run function
# plt <- plot_promoters(dat_in=dt_ov_ls_cp1, key="1700016C15Rik_ENSMUSG00000015962.5+", statistic = fraction_modified, cov_col=NULL)
# ggsave(plt, file=paste0("test_promoter2.pdf"), height = 12, width = 30)

#plot and store as list
plots_list <- imap(dt_ov_ls, ~plot_promoters(dat_in=.x , key=.y, statistic = fraction_modified, cov_col=NULL))

#save as rds on disk
saveRDS(plots_list, file = paste0(opt$plots_rds))

#Export as pdf
#pdf(file=paste0("test_methylation_plugNplay.pdf"), height = 12, width = 30)
pdf(file=paste0(opt$plots_pdf),height = 12, width = 30)
print(plots_list)
dev.off()