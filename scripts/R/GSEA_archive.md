#samuel ahuno


#load libraries
library(clusterProfiler)

#gene set enrichment analysis

## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrezid


## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)
head(foldchanges)


## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
              organism = "hsa", # supported organisms listed below
              nPerm = 1000, # default number permutations
              minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
              pvalueCutoff = 0.05, # padj cutoff value
              verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
# Write results to file
write.csv(gseaKEGG_results, "results/gseaOE_kegg.csv", quote=F)


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


# DO NOT RUN
BiocManager::install("msigdbr")
library(msigdbr)
library(msigdbr)
msigdbr_show_species()


# Use a specific collection; example C6 oncogenic signatures, "Mus musculus" for mouse
#for cancer_hallmark gene sets use "H1"
#adapt: landau methods; "C2" CGP-expression signatures of genetic and chemical pertubaions and CP- canonical pathways derived from kegg, reactome, and biocarta
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

# Run GSEA
msig_GSEA <- GSEA(foldchanges, TERM2GENE = m_t2g, verbose = FALSE)

#make dotplot
require(DOSE)
dotplot(msig_GSEA, showCategory=10, split=".sign") + facet_grid(.~.sign)