#given a grangeslist, get the cpgs in the granges

getMotif_sites <- function(ref_genome=NULL, grList=NULL, context="CG"){
  # Load the required genome
  genome <- switch(ref_genome,
                   mm10 = BSgenome.Mmusculus.UCSC.mm10,
                   hg38 = BSgenome.Hsapiens.UCSC.hg38,
                   stop("Invalid ref_genome specified")
  )

#get seq of the granges in question 
    seqGrList <- lapply(grList, function(x) {
        # print(seqnames(x))
        seqGr <- getSeq(genome, x) #get seqs from ref genome
        motif_locations <- Biostrings::vmatchPattern(context, seqGr) #get cpg sites
# print(str(seqGr))
        # Calculate the starts based on the ends and the width
ends <- motif_locations@ends[[1]]
width <- motif_locations@width0
starts <- ends - width + 1

# Convert to GRanges
gr <- GRanges(seqnames = seqnames(x), # Assuming a dummy sequence name "seq1"
              ranges = IRanges(start = starts, end = ends))

        return(gr)
        })


        return(GRangesList(seqGrList))
}

# dput(out_cps[[1]])
# structure(out_cps[[1]])
# as.granges(out_cps[[1]])
# GRangesList(out_cps)
#test run
# test_grList <- lapply(dt_ov_ls[1:2], function(x) {
# Null_ov <- x[, .(prom_chrom, prom_start,  prom_end, mod_base = "*", score = 0, prom_strand, key, Nvalid_cov=0,fraction_modified=0.0)] %>% head(1) %>% as.data.frame()
# setnames(Null_ov, c("prom_chrom", "prom_start", "prom_end", "prom_strand"), c("seqnames", "start", "end", "strand"))
# rng <- Null_ov %>% as_granges()
# return(rng)
# })


# out_cps <- getMotif_sites(ref_genome="mm10", grList=test_grList, context="CG")
# out_cps[[1]]

