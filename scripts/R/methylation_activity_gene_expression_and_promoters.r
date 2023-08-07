### plot gene expression and promoters 

library(data.table)
library(GGally)

##read multiple files
# library(data.table)
path_methyl_rate <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/gather_files"
methyl_rate_files <- list.files(path_methyl_rate, full.names = TRUE, pattern = "*.data_methylation_gene_promoters_metrics_all_chroms.txt$", recursive = TRUE)
methyl_rate_dt_ls <- lapply(methyl_rate_files, function(x) fread(x))
lapply(methyl_rate_dt_ls, function(x) setnames(x, c("key", "median_promoter_methyl" , "methyl_promoter_entropy", "nCpGs_promoter_observed", "nGCs_promoter_expected", "nCpGs_promoter_expected","proportion_methyl_promoter", "seqnames")))
names(methyl_rate_dt_ls) <- gsub(".data.*", "", basename(methyl_rate_files)) #add file names

#setnames(x, paste0(opt$stats2Use,"_Prob_Meth_perSite"), "statistic")
#names(x) <- c("key", "median_promoter_methyl" , "methyl_promoter_entropy", "nCpGs_promoter_observed", "nGCs_promoter_expected", "nCpGs_promoter_expected","proportion_methyl_promoter", "seqnames")
#lapply(methyl_rate_dt_ls, function(x) dim(x))
methyl_rate_dt <- rbindlist(methyl_rate_dt_ls, idcol = "sample")
head(methyl_rate_dt)



#read gene expression data
### To do: merge with expression data in different script
expr_data <- "/juno/work/greenbaum/users/ahunos/rotation/data/dds_normalized_cnts_parp_inhibitor.tsv"
expr_dt <- fread(expr_data)

#convert expression data to long format
expr_dt_melt <- melt(expr_dt, id.vars = "ensgene", variable.name = "sample", value.name = "Gene_expr")
#optional rename sample  names
expr_dt_melt[,`:=`(sample = paste0("BRCA_13135_",sample))]
setkey(expr_dt_melt, ensgene, sample)

#merge expression and methylation data
#methyl_rate_dt[expr_dt_melt, ]
methyl_rate_dt[,`:=`(ensgene = gsub(".*ENSG","ENSG",key))][,`:=`(ensgene = gsub("\\+|\\-*$","",ensgene))]
setkey(methyl_rate_dt, ensgene, sample)
methyl_rate_and_exprs_dt <- methyl_rate_dt[expr_dt_melt, nomatch = NULL]
#methyl_rate_dt[198086,]
key(methyl_rate_and_exprs_dt)
setkey(methyl_rate_and_exprs_dt, ensgene, sample)

#find correlation between genes and methylation
# methyl_rate_and_exprs_dt[,.(exprXmedMethyl = cor(Gene_expr, median_promoter_methyl))]
# methyl_rate_and_exprs_dt[,.(exprXmedMethyl = cor.test(Gene_expr, median_promoter_methyl),
#                             exprXproportionMethyl = cor.test(Gene_expr, proportion_methyl_promoter))]
#cor.test(methyl_rate_and_exprs_dt$Gene_expr, methyl_rate_and_exprs_dt$median_promoter_methyl)
methyl_rate_and_exprs_dt[,head(.SD, 2), by = list(sample)]
methyl_rate_and_exprs_dt %>% arrange(desc(Gene_expr)) %>% group_by(sample) %>% slice(1:2) %>% ungroup() %>% as.data.table() #what is methyllation activity of top expressed genes
# expr_dt[ensgene %like% "ENSG00000178605"] #this supports the fact that it makes sense to merge on `ensgene.version` since dseq object had `ENSG00000178605.13` & `ENSG00000178605.13_PAR_Y` in it's count matrix 
# expr_dt[,`:=`(ensgene2 = gsub("\\..*","",ensgene))][,`:=`(num2 = .N),by=ensgene2][num2>1,]

# gsub(".*_|\\+|-|*$","","DUSP22_ENSG00000112679.14+")
# gsub(".*_|\\+|-|\\*$","","DUSP22_ENSG00000112679.14-")
#methyl_rate_promoters_dt[key =="DUSP22_ENSG00000112679.14",]
#methyl_rate_promoters_exprs_dt <- methyl_rate_dt[,`:=`(ensgene = gsub(".*_","",key))][expr_dt, on = "ensgene", nomatch = NULL]
######### methyl_rate_promoters_exprs_dt <- methyl_rate_promoters_dt[,`:=`(ensgene = gsub(".*_","",gsub("\\..*","",key)))][expr_dt[,`:=`(ensgene = gsub("\\..*","",ensgene))], on = "ensgene", nomatch = NULL] #don't use, cos no specific mapping of exp and methylation data
# sum(!is.na(methyl_rate_promoters_exprs_dt$key))
#setkey(methyl_rate_promoters_exprs_dt, key)
# install.packages("GGally")
#duplicated(methyl_rate_promoters_exprs_dt, by="key")

methyl_rate_and_exprs_dt_ls <- split(methyl_rate_and_exprs_dt,  methyl_rate_and_exprs_dt$sample) #split by sample

plt_ggairs_diagnostic <- function(DT, sample_name){
    plt_ggpairs <- ggpairs(DT[,.(Gene_expr = Gene_expr, 
                                                med_methyl=median_promoter_methyl, 
                                                entropy_methyl = methyl_promoter_entropy, 
                                                nCpGs_obs = nCpGs_promoter_observed, 
                                                nGCs_exp=nGCs_promoter_expected, 
                                                nCpGs_exp=nCpGs_promoter_expected, 
                                                met_prop=proportion_methyl_promoter)], 
                       lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)), 
                       diag = list(continuous = wrap("barDiag", binwidth = 0.1, alpha = 0.3, size = 0.5)), 
                       title = paste0("Gene Promoter Methylation activity - ", sample_name)) + 
                       theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                       theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    return(plt_ggpairs)
}

#plot and store as list
plots_list_ggpairs <- imap(methyl_rate_and_exprs_dt_ls, ~plt_ggairs_diagnostic(DT=.x , sample_name=.y))
#save as rds on disk
#save plots as pdf
pdf(file=paste0("methylation_rate_promoters_ggpairs_BRCA_parp_inhibitor.pdf"),height = 12, width = 30)
print(plots_list_ggpairs)
dev.off()


library(ggpubr)
func_plt_expression_methyly <- function(DT, sample_name){
    plt_expression_methyly <- DT[,.(Gene_expr = Gene_expr, 
                                                med_methyl=median_promoter_methyl, 
                                                entropy_methyl = methyl_promoter_entropy, 
                                                nCpGs_obs = nCpGs_promoter_observed, 
                                                nGCs_exp=nGCs_promoter_expected, 
                                                nCpGs_exp=nCpGs_promoter_expected, 
                                                met_prop=proportion_methyl_promoter)] %>% 
                                                ggscatter( x="Gene_expr", y="med_methyl",
   #add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson", label.x = 10000, label.y = 0.7) + 
    labs(title = paste0("Gene Promoter Methylation activity - ", sample_name))


                        #ggplot(aes(x=Gene_expr, y=med_methyl)) + geom_point() +
                    #    lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)), 
                    #    diag = list(continuous = wrap("barDiag", binwidth = 0.1, alpha = 0.3, size = 0.5)), 
                    #    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                    #    theme(panel.grid.major = element_blank(),
                    #     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    return(plt_expression_methyly)
}
plots_list_expr_meth <- imap(methyl_rate_and_exprs_dt_ls, ~func_plt_expression_methyly(DT=.x , sample_name=.y))
pdf(file=paste0("methylation_rate_promoters_and_Expression_BRCA_parp_inhibitor.pdf"),width = 12, height = 9)
print(plots_list_expr_meth)
dev.off()
#ggsave(plt_ggpairs, file="test_methylation_rate_promoters.pdf", width = 12, height=7)

