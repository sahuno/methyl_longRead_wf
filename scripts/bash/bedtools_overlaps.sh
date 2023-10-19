#!/bin/bash
#find overlaps between methyl and gene promoters
#input should sorted bed files
#input $1 = methylation file
#input $2 = gene promoters file
#input $3 = output file


module load bedtools
bedtools intersect -a $1 -b $2 -wb >> $3 

# #Test case: bedtools genes and files
# sort -k1,1 -k2,2n /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/sandbox/D-A-3.chr11.bed > D-A-3.chr11.sorted.bed
# sort -k1,1 -k2,2n /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/sandbox/gene_promoters_encode1kb_mm10_head.tsv > gene_promoters_encode1kb_mm10_head_sorted.tsv
# bedtools intersect -a D-A-3.chr11.sorted.bed -b gene_promoters_encode1kb_mm10_head_sorted.tsv > D-A-3.chr11.sorted.intersect.bed 
# bedtools intersect -a D-A-3.chr11.sorted.bed -b gene_promoters_encode1kb_mm10_head_sorted.tsv -wb > D-A-3.chr11.sorted.intersect_wb.bed 

# sort -k1,1 -k2,2n /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/gene_promoters_encode1kb_mm10.tsv > /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/gene_promoters_encode1kb_mm10.sorted.bed


### bedtools operations

#how many genes have valid cpgs
#bedtools groupby /juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/sandbox/D-A-3.chr11.sorted.intersect_wb.bed -g 19-21 -c 10 -o count
