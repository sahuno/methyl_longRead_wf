## filter for sample
## adapted from siyu sun

#load packages
library(DESeq2) #main desdeq analysis
library(tidyverse)
library(EnhancedVolcano) # used for making volcano plots
library(viridis)
library(magrittr)
library("pheatmap")
library(data.table)
library("ggrepel")
library(ggfortify)
library(clusterProfiler)

#read in expresssion data
# work_dir = "/juno/work/greenbaum/projects/TRI_EPIGENETIC/results/RNA_seq"
#setwd(work_dir)

#load datasets from rna-seq results folder - sasha
annot <- fread('annot.tsv') 
counts_annot <- fread('CT/counts_annot.tsv')
cts <- fread('CT/counts.tsv')
#dim(annot)
# unique(annot)
sum(!is.na(annot$entrez.gene.id))
#get protein coding genes only
gs.coding_ids <- annot[gene.type=="protein_coding",.(gene.id)]
cts_coding <- cts[gs.coding_ids, on = "gene.id"] %>% as.data.frame()
cts_coding <- cts_coding %>% column_to_rownames(var = "gene.id") # #convert colnames to rowids

#### repeat landscape
gs.rep <- annot %>% dplyr::filter(chr == "REP") %>% 
    select(c(gene.id, gene.type)) %>% 
      mutate(tp = case_when(str_detect(gene.type,"^(Endogenous|ERV|LTR|Lomg_terminal_repeat|Endogenous_Retrovirus)") ~ "LTR", 
                            str_detect(gene.type,"^(Retroposon/SVA|MSAT)") ~ "SVA", 
                            str_detect(gene.type,"^(LINE|_LINE1|L1|L2|_LINE|RTE|CR1|R4|Penelope|Non-LTR_Retrotransposon)") ~ "LINE", 
                            str_detect(gene.type,"^tRNA") ~ "tRNA", 
                            str_detect(gene.type,"^rRNA") ~ "rRNA", 
                            str_detect(gene.type,"^SINE") ~ "SINE", 
                            str_detect(gene.type,"^(SAT|satellite|Satellite)") ~ "SAT", 
                            str_detect(gene.type,"^(DNA|hAT|_DNA|piggyBac|Mariner/Tc1|MuDR)") ~ "DNA", 
                            str_detect(gene.type,"^(X17_DNA)") ~ "DNA", 
                            str_detect(gene.type,"^(Simple_repeat|Low_complexity)") ~ "Simple_repeat", 
                            str_detect(gene.type,"^(scRNA|snRNA|srpRNA|RNA)") ~ "sRNA", 
                            str_detect(gene.type,"^Pseudogene") ~ "Pseudogene", 
                            TRUE ~ "Unclassified"))
########################################################################################





################################################################################################
### Step 1; read annotation data
##############################################################################################
#read coldata 
samples_df <- read.delim("/juno/work/greenbaum/projects/TRI_EPIGENETIC/metadata_triplicates.txt")
samples_df <- samples_df %>% mutate(sample = gsub("D","R",sample),
                                    sample = gsub("-",".",sample),
                                    condition = condition_short,
                                    condition = case_when(str_detect(condition, "DMSO") ~ "control", TRUE ~ "treatment"),
                                    type = factor(condition_short),
                                    condition = factor(condition)) %>% 
                                    dplyr::select(-c(condition_short,condition_long)) %>%
                                    column_to_rownames("sample")
#names(samples_df) <- c("type", "condition")
unique(samples_df$type)
# ncol(cts_coding)
# nrow(samples_df)
dds <- DESeqDataSetFromMatrix(countData = cts_coding,
                              colData = samples_df,
                              design = ~ type)
#reset levels
dds$condition <- relevel(dds$condition, ref = "control")
dds$type <- relevel(dds$type, ref = "DMSO")

##pre - filtering
keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
dds <- dds[keep,] #filter dds
head(counts(dds))


################################################################################################
############################### Step 2; Run deseq2 analysis
##############################################################################################
#run DESeq analysis - normalization and filtering
dds <- DESeq(dds)
res <- results(dds, alpha=0.05, contrast=c("type","AZCT","DMSO")) #use alpha of 0.05

summary(res)
resultsNames(dds) #double check what factors was used in the desing
#run for all contrasts
contrasts_ls <- lapply(resultsNames(dds)[-1], function(x){x <- gsub("_vs_", "_", x); str_vc <- unlist(str_split(x, pattern="_")); return(str_vc)})
paste(contrasts_ls[[1]], collapse = "_")

##shrink fold changes for MA plots
#resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2, type="apeglm")

pdf("MA_plot_lfcShrink_dmso_vrs_all_treated.pdf")
plotMA(resLFC, ylim=c(-2,2))
dev.off()
###########################################################################

################################################################################################
#### export transformed data for downstream analysis in R
################################################################################################
#get normalized for downstream analysis
dds <- estimateSizeFactors(dds); 
dds_normalized_cnts <- counts(dds, normalized=TRUE) 
dds_normalized_cnts <- dds_normalized_cnts %>% as.data.frame() %>% rownames_to_column("gene.id")
#head(dds_normalized_cnts)
setDT(dds_normalized_cnts)
fwrite(dds_normalized_cnts, file="DESeq2_data/dds_normalized_cnts.tsv", sep="\t")

deseq2rlog_df <- rlog(dds, blind=FALSE)
head(deseq2rlog_df)
head(assay(deseq2rlog_df))
deseq2rlog_assay_df <- assay(deseq2rlog_1) %>% as.data.frame() %>% rownames_to_column("gene.id") #%>% 
head(deseq2rlog_assay_df); dim(deseq2rlog_assay_df)
setDT(deseq2rlog_assay_df)
fwrite(deseq2rlog_assay_df, file="DESeq2_data/DESeqTransform_BlindFALSE.tsv", sep="\t")
################################################################################################




################################################################################################
############################### Step 3; ##annotation deseq results
##############################################################################################
DEResults_ls <- list()
run_multiple_contrasts <- function(contrast_input, dseqObject = dds, shrink = FALSE, padj = 0.05){
#run DESeq analysis - normalization and filtering
contrast_label <- paste(contrast_input, collapse = "_") #used for plot labels

results_per_contrast <- results(dseqObject, alpha=padj, contrast=contrast_input) #use alpha of 0.05
print(resultsNames(results_per_contrast))
# Shrink the log2 fold changes to be more accurate
if(shrink == TRUE){
results_per_contrast <- lfcShrink(dds, 
     contrast=contrast_input, 
     type = "apeglm")	 
     # The coef will be dependent on what your contras was. and should be identical to what
}

dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

saveRDS(results_per_contrast, file=paste0("data/Dseq2ResultsObject_",contrast_label,"_padjust.rds"))

DeSeq2Results_df <- as.data.frame(results_per_contrast) %>% rownames_to_column("gene.id")
# head(DeSeq2Results_df)
#DeSeq2Results_df <- DeSeq2Results_df %>% rownames_to_column(var = "ensembl.gene.id")
DeSeq2Results_df_annot <- left_join(DeSeq2Results_df, annot %>% dplyr::select(c(gene.id, gene.symbol, description,  entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results
# head(DeSeq2Results_df_annot)
# #sum(is.na(DeSeq2Results_df_annot$gene.symbol))
# head(DeSeq2Results_df_annot)

#### do volcano plots
plt_volc_pval <- EnhancedVolcano(DeSeq2Results_df_annot,
    lab = DeSeq2Results_df_annot$gene.symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = 0.05,
    title = paste0(contrast_label,' (p < 0.05)'),
    pointSize = 1, labSize = 5,
    colAlpha = 0.2) +
    theme(plot.subtitle = element_blank())
ggsave(plt_volc_pval, file = paste0("figures/volcano_",contrast_label,"_pvals.pdf"))

#with adjusted p - values
DeSeq2Results_df_annot_padjustOnly <- DeSeq2Results_df_annot %>% dplyr::select(!pvalue) %>% mutate(pvalue = padj) #recreate data frame for ploting
plt_volc_padjust <- EnhancedVolcano(DeSeq2Results_df_annot_padjustOnly,
    lab = DeSeq2Results_df_annot_padjustOnly$gene.symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = 0.05,
    title = paste0(contrast_label,' (BH < 0.05)'),
    widthConnectors = 0.75,
    labSize = 4.0,
    drawConnectors = TRUE,
    ylab = bquote(~Log[10]~ 'p adjusted'),
    pointSize = 1, 
    #labSize = 5,
    colAlpha = 0.2) + theme(plot.subtitle = element_blank())
    
ggsave(plt_volc_padjust, file = paste0("figures/volcano_",contrast_label,"_padjust.pdf"))



#### Identifying significant genes
# Subset to return genes with padj < 0.05
sigLRT_genes <- DeSeq2Results_df_annot %>% 
  dplyr::filter(padj < padj)

# Get number of significant genes
nrow(sigLRT_genes)

# Compare to numbers we had from Wald test
# nrow(sigOE)
# nrow(sigKD)
return(results_per_contrast) #return results for specific contrasts
}


#run for all contrasts
DEResults_ls <- lapply(contrasts_ls, function(x) run_multiple_contrasts(contrast_input=x, dseqObject = dds))
names_contrasts_ls <- unlist(lapply(contrasts_ls, function(x) paste0(x, collapse = "_")))
names(DEResults_ls) <- names_contrasts_ls

#lapply(names_contrasts_ls, function(x) print(paste0(x, " ", dim(DEResults_ls[[x]]))))
################################################################################################
############################### Step #; GSEA
##############################################################################################
gsea_analysis <- function(name_de_res, specie_type = "Mus musculus", category_tag = "H", nCategory_2show = 10, nShowBar=100){
  require(DOSE)
  library(msigdbr)

de_res_df <- DEResults_ls[[name_de_res]] %>% as.data.frame() %>% rownames_to_column(var = "gene.id")

  DEResults_type_trt_ctrl_annt  <- left_join(de_res_df, annot %>% dplyr::select(c(ensembl.gene.id, gene.id, gene.symbol, description, entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results
#gene set enrichment analysis
## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(DEResults_type_trt_ctrl_annt, entrez.gene.id != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez.gene.id) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez.gene.id

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)


#get data base to use
m_t2g <- msigdbr(species = specie_type, category = substitute(category_tag)) %>% 
  dplyr::select(gs_name, entrez_gene) #hall mark genes

print(m_t2g)
# Run GSEA analysis
msig_GSEA_obj <- GSEA(foldchanges, TERM2GENE = m_t2g, verbose = FALSE)

#, category = "C2"
#make dotplot
dotplt_gsea_out <- dotplot(msig_GSEA_obj, showCategory=nCategory_2show, split=".sign") + facet_grid(.~.sign)
ggsave(dotplt_gsea_out, file = paste0("figures/dotPlot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_",nCategory_2show,".pdf"))


# gsea_plot2 <-  enrichplot::gseaplot2(msig_GSEA_H, geneSetID="HALLMARK_MYC_TARGETS_V1")
# ggsave(gsea_plot2, file = paste0("figures/gsea_msigdbr_",name_de_res,"_type_CKI_DMSO_top_",nCategory_2show,".pdf"))

#ridge plot; enrichment distribution
enrich_dist_plot_out <- ridgeplot(msig_GSEA_obj) + labs(x = "enrichment distribution")
ggsave(enrich_dist_plot_out, file = paste0("figures/ridgePlot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_",nCategory_2show,".pdf"), width=18, height=27)


# msig_GSEA_results <- msig_GSEA_obj@result
# #rank gene sets by normalized  enrichment scores and q-values
# msig_GSEA_results_sig_qval <- msig_GSEA_results %>% dplyr::filter(qvalues < 0.05) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA

# lollipop_gsea <- msig_GSEA_results_sig_qval %>% top_n(10) %>% ggplot(aes(x=ID, y = NES)) +   geom_segment(
#     aes(x=ID, xend=ID, y=0, yend=NES)) +
#   geom_point()  +  coord_flip() +
#   labs(title = "normalized enrichment scores (NES) for significant gene sets (q-value < 0.05)")
# ggsave(lollipop_gsea, file = paste0("figures/lollipop_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_",nCategory_2show,".pdf"), width=18, height=27)

#return(msig_GSEA_obj)
}

#paste0("figures/gsea_msigdbr_H_type_CKI_DMSO_EnrichDistrib.pdf")
#run gsea for all contrasts
lapply(names_contrasts_ls, function(x) gsea_analysis(x))










DEResults_type_trt_ctrl <- as.data.frame(DEResults_ls[["type_CKI_DMSO"]]) %>% rownames_to_column(var = "gene.id")
DEResults_type_trt_ctrl_annt  <- left_join(DEResults_type_trt_ctrl, annot %>% dplyr::select(c(ensembl.gene.id, gene.id, gene.symbol, description, entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results
head(DEResults_type_trt_ctrl_annt)
#sum(!is.na(DEResults_type_trt_ctrl_annt$entrez.gene.id))

#gene set enrichment analysis

## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(DEResults_type_trt_ctrl_annt, entrez.gene.id != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez.gene.id) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez.gene.id


## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)
head(foldchanges)



# DO NOT RUN
#BiocManager::install("msigdbr")
library(msigdbr)
msigdbr_species()

msigdbr(species = "Mus musculus") %>% group_by(gs_cat) %>% summarise(n_gs = n()) %>% arrange(desc(n_gs))

# Use a specific collection; example C6 oncogenic signatures, "Mus musculus" for mouse
#for cancer_hallmark gene sets use "H1"
#adapt: landau methods; "C2" CGP-expression signatures of genetic and chemical pertubaions and CP- canonical pathways derived from kegg, reactome, and biocarta
m_t2g_C2 <- msigdbr(species = "Mus musculus", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g_H1 <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene) #hall mark genes
# Run GSEA
msig_GSEA_C2 <- GSEA(foldchanges, TERM2GENE = m_t2g_C2, verbose = FALSE)
msig_GSEA_H <- GSEA(foldchanges, TERM2GENE = m_t2g_H1, verbose = FALSE)

#, category = "C2"
#make dotplot
require(DOSE)
dp_out_C2 <- dotplot(msig_GSEA_C2, showCategory=10, split=".sign") + facet_grid(.~.sign)
dp_out_H <- dotplot(msig_GSEA_H, showCategory=10, split=".sign") + facet_grid(.~.sign)

ggsave(dp_out_C2, file = paste0("figures/gsea_msigdbr_C2_type_CKI_DMSO_top10.pdf"))
ggsave(dp_out_H, file = paste0("figures/gsea_msigdbr_H_type_CKI_DMSO_top10.pdf"))

#gsea plot
gsea_plot <- gseaplot(msig_GSEA_H, geneSetID="HALLMARK_MYC_TARGETS_V1")
gsea_plot2 <-  enrichplot::gseaplot2(msig_GSEA_H, geneSetID="HALLMARK_MYC_TARGETS_V1")
# gsea_rank_plot <- enrichplot::gsearank(msig_GSEA_H1, geneSetID="KEGG_RIBOSOME") #not necessary

ggsave(gsea_plot, file = paste0("figures/gsea_msigdbr_H_type_CKI_DMSO_ES.pdf"))
ggsave(gsea_plot2, file = paste0("figures/gsea_msigdbr_H_type_CKI_DMSO_ES_2.pdf"))
# ggsave(gsea_rank_plot, file = paste0("figures/gseaRankPlot_msigdbr_H1_type_CKI_DMSO.pdf")) #not necessary

#ridge plot; enrichment distribution
enrich_dist_plot <- ridgeplot(msig_GSEA_H) + labs(x = "enrichment distribution")
ggsave(enrich_dist_plot, file = paste0("figures/gsea_msigdbr_H_type_CKI_DMSO_EnrichDistrib.pdf"), width=18, height=27)



#save some results
msig_GSEA_H_results <- msig_GSEA_H@result
head(msig_GSEA_H_results)
#rank gene sets by normalized  enrichment scores and q-values
msig_GSEA_H_results %>% arrange(qvalues, desc(abs(NES))) %>% head(10)
msig_GSEA_H1_results_sig_qval <- msig_GSEA_H1_results %>% dplyr::filter(qvalues < 0.05) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA

lollipop_gsea <- msig_GSEA_H1_results_sig_qval %>% top_n(100) %>% ggplot(aes(x=ID, y = NES)) +   geom_segment(
    aes(x=ID, xend=ID, y=0, yend=NES), 
    # color=ifelse(data$x %in% c("A","D"), "orange", "grey"), 
    # size=ifelse(data$x %in% c("A","D"), 1.3, 0.7)
  ) +
  geom_point(
    # color=ifelse(data$x %in% c("A","D"), "orange", "grey"), 
    # size=ifelse(data$x %in% c("A","D"), 5, 2)
  )  +  coord_flip() +
  labs(title = "normalized enrichment scores (NES) for significant gene sets (q-value < 0.05)")
ggsave(lollipop_gsea, file = paste0("figures/gsea_msigdbr_H1_type_CKI_DMSO_lollipop.pdf"), width=18, height=27)

#msig_GSEA_H1_results %>% arrange(desc(abs(NES)), qvalues) %>% head(10)


# str(msig_GSEA_H1)
slotNames(msig_GSEA_H1)
head(msig_GSEA_H1@result)
names(head(msig_GSEA_H1@result))

# categorySize can be either 'pvalue' or 'geneNum'
#cnetplot(dp_out_H1, categorySize="pvalue", foldChange=foldchanges)


#do enrichment maps
gsea_plot <- function(x){
emp <- emapplot(x)
return(emp)
}









## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
              organism = "mmu", # supported organisms listed below
              nPerm = 1000, # default number permutations
              minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
              pvalueCutoff = 0.05, # padj cutoff value
              verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
# Write results to file
write.csv(gseaKEGG_results, "data/gseaOE_kegg.csv", quote=F)


## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa03008')

#visualize pathways
detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts

## Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
              pathway.id = "hsa03008",
              species = "hsa",
              limit = list(gene = 2, # value gives the max/min limit for foldchanges
              cpd = 1))
## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
   pathview(gene.data = foldchanges, 
            pathway.id = gseaKEGG_results$ID[x], 
            species = "hsa",
            limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), 
           get_kegg_plots)



























################################################################################################
############################### Step #; ##GSEA of deseq2 results








#######################################
################ Visualizations
#########################################
#make MA plots ######
pdf("MA_plot.pdf")
plotMA(res)
dev.off()

#plot genes
# pdf("ENSG00000081277_genePlot_QC.pdf")
# plotCounts(dds, gene="ENSG00000081277.13", intgroup = "condition")
# dev.off()

#volcano plots

########################################################################
###########generating normlaized counts
# Transform count data using the variance stablilizing transform
deseq2VST_1 <- vst(dds, blind=TRUE)


#head(assay(deseq2VST_1))

# Convert the DESeq transformed object to a data frame
deseq2VST_assay <- assay(deseq2VST_1)
deseq2VST <- as.data.frame(deseq2VST_assay)
deseq2VST$`gene.id` <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(DeSeq2Results_df[DeSeq2Results_df$padj <= .05 & abs(DeSeq2Results_df$log2FoldChange) > 1, ])
DeSeq2Results_df %>% dplyr::filter(padj <= .05 & abs(log2FoldChange) > 1)


DeSeq2Results_significantOnlyDf <- DeSeq2Results_df_annot %>% dplyr::filter(!gene.symbol == "" & padj <= .05 & abs(log2FoldChange) > 1) 

#join annotation and normalized counts of signifcant genes with lfc >= 1
deseq2VST_sig <- left_join(DeSeq2Results_significantOnlyDf %>% 
                              dplyr::select(c(gene.id,gene.symbol)), 
                          deseq2VST, by="gene.id") #merge gene annotations and deseq results
head(deseq2VST_sig)
deseq2VST_sig_longer <- deseq2VST_sig %>% pivot_longer(!c("gene.id","gene.symbol"), names_to = "samples", values_to = "normalized_counts")
head(deseq2VST_sig_longer)


 
#deseq2VST_sig_longer %>% summarise(tnorm = log2(normalized_counts)) %>% summarize(range(tnorm))
# Make a heatmap
heatmap_plot <- ggplot(deseq2VST_sig_longer, 
                          aes(x=samples, y=gene.symbol, fill=normalized_counts)) + 
                            geom_raster() + scale_fill_viridis() + 
                            labs(title = "Heatmap of Differential Gene Expression (All treatment vrs Ctrl)") + 
                              theme(axis.text.x=element_text(angle=65, hjust=1), 
                                    legend.position = "top")
#heatmap_plot
ggsave(heatmap_plot, filename = "heatmap_rna_seq.pdf")


########################################################################
###############PCA plots of transformed counts
#pdf("/juno/work/greenbaum/users/ahunos/rotation/figures/rna_seq/pca_samples.pdf")
ppca <- plotPCA(deseq2VST_1, intgroup=c("condition"))
pp2 <- ppca + 
labs(title = "PCA of samples (treatment. vrs Ctrl)") 
ggsave(pp2, filename= "pca_samples.pdf", width=7, height=5)
#dev.off()

pca <- prcomp(t(assay(deseq2VST_1)), scale. = TRUE)
pdf("pc_ggfortify.pdf")
autoplot(pca, data= samples_df, shape=FALSE, colour = 'condition')
dev.off()





deseq2VST_sig4Clustering <- deseq2VST_sig %>% dplyr::select(!c(gene.id))
rownames(deseq2VST_sig4Clustering) <- make.names(deseq2VST_sig4Clustering[,1], unique = TRUE)
deseq2VST_sig4Clustering <- deseq2VST_sig4Clustering %>% dplyr::select(!"gene.symbol")
deseq2VST_sig4Clustering_mat <- data.matrix(deseq2VST_sig4Clustering)
class(deseq2VST_sig4Clustering_mat)
annot_coldf <- as.data.frame(colData(dds)[,"type", drop=FALSE])


pdf("clusteringAnalysis.pdf")
pheatmap(deseq2VST_sig4Clustering_mat, main = "Differentially Expressed trt Vrs Ctrl)", annotation_col=annot_coldf)
dev.off()

#save.image(file="/juno/work/greenbaum/users/ahunos/rotation/data/RNA_seq_removing_outlierLowQualitySasmple.RData") 