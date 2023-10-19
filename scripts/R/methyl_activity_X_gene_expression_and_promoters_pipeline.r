#author samuel ahuno
#date august 8th 2023
#purpose: plot gene expression and promoters methylation


#load libraries
library(data.table)
library(GGally)
library(optparse)
library(tidyverse)


option_list <- list(
   make_option(("--path_methyl_rate"), type = "character", default="/juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/results/gather_files"),
   make_option(("--Diff_Expr"), type = "character", default="/juno/work/greenbaum/projects/TRI_EPIGENETIC/RNASeq_DE_TriEpi/data/condition_AZCT_vs_DMSO/Dseq2Results_condition_AZCT_vs_DMSO.tsv"), #results would look this way now
make_option(("--RNA_trans_counts"), type = "character", default="/juno/work/greenbaum/projects/TRI_EPIGENETIC/RNASeq_DE_TriEpi/data/condition_AZCT_vs_DMSO/Dseq2Results_condition_AZCT_vs_DMSO.tsv"), #results would look this way now

   make_option(("--orgn"), type = "character", default="mouse")
                )

opt <- parse_args(OptionParser(option_list = option_list))



##########################################################################################################################
### step 1: read methylation data
#path_methyl_rate <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/gather_files"
methyl_rate_files <- list.files(opt$path_methyl_rate, 
                                full.names = TRUE, 
                                pattern = "*data_all_chroms.tsv$",
                                recursive = TRUE)##read multiple files
methyl_rate_dt_ls <- lapply(methyl_rate_files, function(x) fread(x)) #read merged methylation activity data

# #add header to file
# lapply(methyl_rate_dt_ls, function(x) setnames(x, c("key", "median_promoter_methyl" , "methyl_promoter_entropy", 
#                                                     "nCpGs_promoter_observed", 
#                                                     "nGCs_promoter_expected", "nCpGs_promoter_expected",
#                                                     "proportion_methyl_promoter", "seqnames")))

names(methyl_rate_dt_ls) <- gsub(".methyl.*", "", basename(methyl_rate_files)) #add file names
methyl_rate_dt <- rbindlist(methyl_rate_dt_ls, idcol = "sample")
#head(methyl_rate_dt)
# methyl_rate_dt %>% group_by(sample, gene.id) %>% summarize(n())

#"key"
# geom_mean_promoter_methyl 
# methyl_promoter_entropy_Avg
# methyl_promoter_entropy_shann
# nCpGs_promoter_observed
# nGCs_promoter_expected
# nCpGs_promoter_expected
# proportion_methyl_promoter 


#make standard gene ids 
#methyl_rate_dt[expr_dt_melt, ]
# if(opt$orgn == "human"){
#     methyl_rate_dt[,`:=`(gene.id = gsub(".*ENSG","ENSG",key))][,`:=`(gene.id = gsub("\\+|\\-*$","",gene.id))]
# setkey(methyl_rate_dt, gene.id, sample)
# } else if (opt$orgn == "mouse"){
#     methyl_rate_dt[,`:=`(gene.id = gsub(".*ENSMUSG","ENSMUSG",key))][,`:=`(gene.id = gsub("\\+|\\-*$","",gene.id))] 
#     setkey(methyl_rate_dt, gene.id, sample)
# } else {
#     stop("organism not supported")
# }


##########################################################################################################################
######## read expresion data
#regardless there should be a column called ensgene or ensemble_gen
### To do: merge with expression data in different script


#deseq2 object
DE_dt <- fread(opt$Diff_Expr)
methyl_and_DE_dt <- methyl_rate_dt[DE_dt, on="gene.id", nomatch=NULL] #merge on methylation rate and gene expression
fwrite(methyl_and_DE_dt, file="methyl_and_DE_dt.tsv", sep="\t", row.names = FALSE, col.names = TRUE)
#methyl_and_DE_dt %>% dplyr::filter(gene.symbol %in% c("GAPDH", "ACTB", "B2M", "HPRT1", "RPL13A", "RPLP0", "TBP", "PPIA", "SDHA", "UBC"))

#mouse; all tissues
# https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-8-127/tables/1
methyl_and_DE_dt_houseKeeping <- methyl_and_DE_dt %>% dplyr::filter(gene.symbol %in% c("Eef2", "Rpl37", "Rpl38", "Leng8", "Eif3k","proteasome (prosome, macropain) 26S subunit", "Gapdh")) %>%
dplyr::select(c("sample","key","geom_mean_promoter_methyl",
            "methyl_promoter_entropy_shann", "nCpGs_promoter_observed_data",
            "nGCs_promoter_expected_ref",   "nCpGs_promoter_expected_ref",
                "prop_methyl_in_promoter_data","log2FoldChange", "pvalue","padj")) #%>% pivot_wider(names_from = "sample", values_from = "log2FoldChange")


write_tsv(methyl_and_DE_dt_houseKeeping, file="methyl_and_DE_dt_houseKeeping.tsv")
methyl_and_DE_dt_houseKeeping_ls <- split(methyl_and_DE_dt_houseKeeping, methyl_and_DE_dt_houseKeeping$sample)
#methyl_and_DE_dt_houseKeeping_ls[[1]][,-1]
#lapply(methyl_and_DE_dt_houseKeeping_ls, function(x){x <- x %>% select(c())})
saveRDS(methyl_and_DE_dt_houseKeeping_ls, file="methyl_and_DE_dt_houseKeeping_ls.rds")


methyl_and_DE_dt_houseKeeping_ls_split_genes <- split(methyl_and_DE_dt_houseKeeping, methyl_and_DE_dt_houseKeeping$key)
saveRDS(methyl_and_DE_dt_houseKeeping_ls_split_genes, file="methyl_and_DE_dt_houseKeeping_ls_split_genes.rds")

write_tsv(methyl_and_DE_dt_houseKeeping_ls_split_genes[["Gapdh_ENSMUSG00000057666.18-"]] %>% dplyr::select(!c(sample,key, pvalue,padj)), 
file="methyl_and_DE_dt_Gapdh_ENSMUSG00000057666_1.tsv")

# methyl_and_DE_dt %>% dplyr::filter(nCpGs_promoter_observed_data >= 5)
#remove genes with less than 10 CpGs
# methyl_and_DE_dt %>% group_by(sample, gene.id) %>% summarize(n())
#filter top genes
# methyl_and_DE_dt[log2FoldChange > 2 & padj > 0.05][order(log2FoldChange),by = sample][, head(.SD, 2), by = sample]
# methyl_and_DE_dt[log2FoldChange > 2 & padj > 0.05, ,by = "sample"]
plt_expr_methyl <- function(minObservedCpGs = 10, abs_l2fc=2, p_adj=0.01, nGenes=10){
methyl_and_DE_dt <- methyl_and_DE_dt[nCpGs_promoter_observed_data >= minObservedCpGs,]

topSig_DE_UP = methyl_and_DE_dt  %>% group_by(sample) %>% 
        dplyr::filter(log2FoldChange >= abs_l2fc & padj < p_adj) %>% 
        arrange(desc(log2FoldChange)) %>% dplyr::slice(1:nGenes) %>% 
        ungroup() %>%dplyr::mutate(group = "DE UP") %>% 
        as.data.table() #%>%pull(gene.symbol) ##what is methyllation activity of top expressed genes

topSig_DE_DwN =  methyl_and_DE_dt %>% group_by(sample) %>% 
dplyr::filter(log2FoldChange <= -(abs_l2fc) & padj < p_adj) %>%
 dplyr::slice_min(log2FoldChange, n=nGenes) %>% 
 ungroup() %>%dplyr::mutate(group = "DE DOWN") %>% 
 as.data.table() #what is methyllation activity of top expressed genes

top_Methyl = methyl_and_DE_dt  %>% group_by(sample) %>% 
arrange(desc(geom_mean_promoter_methyl)) %>% 
dplyr::slice(1:nGenes) %>% ungroup() %>% 
dplyr::mutate(group = "top methyl") %>% 
as.data.table() #%>%pull(gene.symbol) ##what is methyllation activity of top expressed genes

lowest_Methyl = methyl_and_DE_dt  %>% group_by(sample) %>% 
arrange(desc(geom_mean_promoter_methyl)) %>% 
dplyr::slice_min(geom_mean_promoter_methyl, n=nGenes) %>% 
ungroup() %>%dplyr::mutate(group = "lowest methyl") %>% 
as.data.table() #%>%pull(gene.symbol) ##what is methyllation activity of top expressed genes

de_me_dt <- rbindlist(list(topSig_DE_UP, topSig_DE_DwN, top_Methyl, lowest_Methyl), fill=TRUE)

library(ggrepel)
plt_de_me <- ggplot(data=de_me_dt, 
                aes(x=geom_mean_promoter_methyl, 
                    y= log2FoldChange, 
                    color=group, 
                    label=gene.symbol)) + 
            geom_point() + 
            facet_wrap(~sample) +   
            geom_text_repel() +
            theme_classic(base_size = 16) +
            labs(title = "Promoter Methyl Vrs DExpr ", 
                subtitle ="Each CpG supported by >= 5 reads in 75% of cohorts",
                x="Promoter hypergeometeric mean (>=5 CpGs in Promoter)", 
                y=paste0("log2FC;ACZT vs DMSO, (abs(log2FC) >= 2, padj<", p_adj))

ggsave(plt_de_me, file=paste0("figures_methyl/methylation_rate_promoters_DE_genes_N",nGenes,".pdf"), 
            width = 12, height=7)


plt_de_me_box <- ggplot(data=de_me_dt, 
                            aes(y=geom_mean_promoter_methyl, 
                                x=as.factor(group), fill=group,
                                group=as.factor(group))) + 
                geom_boxplot() + 
                geom_jitter(aes(group=as.factor(group))) + 
                facet_wrap(~sample) +
                labs(title = "Promoter Methyl ", 
                    subtitle ="Each CpG supported by >= 5 reads in 75% of cohorts",
                        y=paste0("Promoter hypergeometeric mean (>=", minObservedCpGs," CpGs in Promoter)")) +
geom_text_repel(aes(label=as.factor(gene.symbol), group=as.factor(group)), 
        segment.alpha = 0, position = position_dodge(width = 0.7)) 

ggsave(plt_de_me_box, file=paste0("figures_methyl/methylation_rate_promoters_DE_genes_box_N",nGenes,".pdf"), width = 12, height=7)  
}

nGenesSim <- seq(from=10, to=100, by=5)
lapply(nGenesSim, function(x) {plt_expr_methyl(minObservedCpGs = 10, abs_l2fc=2, p_adj=0.05, nGenes=x)})
# plt_expr_methyl(minObservedCpGs = 10, abs_l2fc=2, p_adj=0.01, nGenes=10)


##global methylation rate
plt_global_methyl <- ggplot(data=methyl_and_DE_dt, 
                            aes(y=geom_mean_promoter_methyl, 
                                x=as.factor(sample), 
                                # fill=sample,
                                group=as.factor(sample))) + 
                geom_boxplot() + 
                # geom_jitter(aes(group=as.factor(sample))) + 
                # facet_wrap(~sample) +
                labs(title = "Promoter Methyl ", 
                    subtitle ="Each CpG supported by >= 5 reads in 75% of cohorts",
                        y="Promoter hypergeometeric mean (>=10 CpGs in Promoter)") 
ggsave(plt_global_methyl, file=paste0("figures_methyl/global_methylation_rate_promoters",".pdf"), width = 12, height=7)  




#########################################################################################################
## below is code for transformed de
#########################################################################################################
methyl_rate_dt[,.N,by=sample]
# opt$input_expression
#expr_dt <- fread(opt$input_expression)
readExpData <- function(path=NULL, gene_col = "gene.id", methyl_rate_data=methyl_rate_dt,
                    rename_samples = TRUE, 
                    # old_names=c("P_1", "P_2","P_3","P_4","P_5","P_6"), 
                    # old_names=c("P_1", "P_2","P_3","P_4","P_5","P_6"), 
                    old_names=c("Parp_1", "Parp_2","Parp_3","Ctrl_4","Ctrl_5","Ctrl_6"),
                    new_names=c("BRCA_13135_Parp_1", "BRCA_13135_P_2","BRCA_13135_P_3","BRCA_13135_P_4","BRCA_13135_P_5","BRCA_13135_P_6")
                    ){

expr_dt <- fread(path)
#names(expr_dt)[6] <- "gene.id"
#print(gene_col)
# if(!any(names(expr_dt) %like% gene_col)){
#   print("renaming gene column as ensgene; thus esemble gene id.version number")
#   setnames(expr_dt, get(gene_col), "ensgene")
#   }
#optional rename sample  names
if(rename_samples){
    #expr_dt_melt[,`:=`(sample = paste0("BRCA_13135_", sample))]
    setnames(expr_dt, old_names, new_names, skip_absent=TRUE)
}

#print(head(expr_dt))
#convert expression data to long format
expr_dt_melt <- melt(expr_dt, id.vars = "gene.id", variable.name = "sample", value.name = "Gene_expr")
#unique(expr_dt_melt$sample)
setkey(expr_dt_melt, gene.id, sample)
#merge expression and methylation data


# print(head(expr_dt))
# print(head(methyl_rate_data))
methyl_rate_and_exprs_dt <- methyl_rate_data[expr_dt_melt, nomatch = NULL] #merge on methylation rate and gene expression
#methyl_rate_dt[198086,]
# key(methyl_rate_and_exprs_dt)
setkey(methyl_rate_and_exprs_dt, gene.id, sample)
return(methyl_rate_and_exprs_dt)
}

methyl_rate_and_exprs_dt <- readExpData(path=opt$input_expression, gene_col = "gene.id")


#find correlation between genes and methylation
# methyl_rate_and_exprs_dt[,.(exprXmedMethyl = cor(Gene_expr, median_promoter_methyl))]
# methyl_rate_and_exprs_dt[,.(exprXmedMethyl = cor.test(Gene_expr, median_promoter_methyl),
#                             exprXproportionMethyl = cor.test(Gene_expr, proportion_methyl_promoter))]
#cor.test(methyl_rate_and_exprs_dt$Gene_expr, methyl_rate_and_exprs_dt$median_promoter_methyl)

#get top expressed genes
methyl_rate_and_exprs_dt[,head(.SD, 2), by = list(sample)]
methyl_rate_and_exprs_dt %>% arrange(desc(Gene_expr)) %>% 
group_by(sample) %>% 
slice(1:2) %>% ungroup() %>% as.data.table() #what is methyllation activity of top expressed genes
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
                                                methyl_promoter_entropy_Avg,
                                                methyl_promoter_entropy_shann,
                                                geom_mean_promoter_methyl,
                                                # nCpGs_obs = nCpGs_promoter_observed, 
                                                # nGCs_exp=nGCs_promoter_expected, 
                                                # nCpGs_exp=nCpGs_promoter_expected, 
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
pdf(file=paste0("methylation_rate_promoters_ggpairs_BRCA_parp_inhibitor.pdf"),height = 12, width = 15)
print(plots_list_ggpairs)
dev.off()


library(ggpubr)
func_plt_expression_methyly <- function(DT, sample_name){
    plt_expression_methyly <- DT[,.(Gene_expr = log10(Gene_expr + 1), 
                                                geom_mean_promoter_methyl)] %>% 
                                ggscatter( x="Gene_expr", y="geom_mean_promoter_methyl",
   #add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) +     labs(title = paste0("Gene Promoter Methylation activity - ", sample_name), x="log10(Gene Expression)")

#    stat_cor(method = "pearson", label.x = 10000, label.y = 0.7) + 
       #coord_trans(x = "log10")+

    return(plt_expression_methyly)
}
# plots_list_expr_meth <- imap(methyl_rate_and_exprs_dt_ls, ~func_plt_expression_methyly(DT=.x , sample_name=.y))
# pdf(file=paste0("methylation_rate_promoters_and_Expression_BRCA_parp_inhibitor.pdf"),width = 12, height = 9)
# print(plots_list_expr_meth)
# dev.off()
#ggsave(plt_ggpairs, file="test_methylation_rate_promoters.pdf", width = 12, height=7)

