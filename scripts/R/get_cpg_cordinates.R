#given a grangeslist, get the cpgs in the granges

getMotif_sites <- function(ref_genome=NULL, grList=NULL, context="CG"){
        library(Biostrings)
  # Load the required genome
  genome_pkg <- switch(ref_genome,
                   mm10 = "BSgenome.Mmusculus.UCSC.mm10",
                   hg38 = "BSgenome.Hsapiens.UCSC.hg38",
                   stop("Invalid ref_genome specified")
  )

#load ref genome
print(genome)

    # Load the actual genome data
    library(genome_pkg, character.only = TRUE)
    genome_data <- get(genome_pkg)


#get seq of the granges in question x=grList[[1]]
    seqGrList <- lapply(grList, function(x) {
        # print(seqnames(x))
        seqGr <- getSeq(genome_data, x) #get seqs from ref genome
        motif_locations <- Biostrings::vmatchPattern(context, seqGr) #get cpg sites


        # Calculate the starts based on the ends and the width
ends <- motif_locations@ends[[1]]
width <- nchar(context)
starts <- ends - width + 1

# Convert to GRanges
        gr <- GRanges(seqnames = seqnames(x),
                      ranges = IRanges(start = starts, end = ends))

        })
        
        return(GRangesList(seqGrList))
}


#test run
library(plyranges)
test_grList <- lapply(dt_ov_ls[1:2], function(x) {
Null_ov <- x[, .(prom_chrom, prom_start,  prom_end, mod_base = "*", score = 0, prom_strand, key, Nvalid_cov=0,fraction_modified=0.0)] %>% head(1) %>% as.data.frame()
setnames(Null_ov, c("prom_chrom", "prom_start", "prom_end", "prom_strand"), c("seqnames", "start", "end", "strand"))
rng <- Null_ov %>% as_granges()
return(rng)
})


out_cps <- getMotif_sites(ref_genome="mm10", grList=test_grList, context="CG")
out_cps[[1]]

# dput(out_cps[[1]])
# structure(out_cps[[1]])
# as.granges(out_cps[[1]])
# GRangesList(out_cps)

