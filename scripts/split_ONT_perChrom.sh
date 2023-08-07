#!/bin/sh
baseName=$1
BigFile=$2
dest=$3
# tmp_noheader=$4
# tmp_standardChrom=$5

#input=/juno/work/greenbaum/users/lihmj/nanopore_methyl/around_gene_body/PAM57786.per_read_modified_base_calls.gene5kb.txt
#input=/work/greenbaum/projects/ont_pipeline/projects/DIVYA_BRCA/results/methylation/megalodon/13135_P_2/per_read_modified_base_calls.txt
#input=/juno/work/greenbaum/users/ahunos/rotation/data/brca/per_read_modified_base_calls/BRCA_13135_P_1/BRCA_13135_P_1_per_read_modified_base_calls.txt

#set dir
FILE=${dest}/${baseName}

# tail -n +2 ${BigFile} > $tmp_noheader #remove header file
# awk '$2 !~ /_/' $tmp_noheader > ${tmp_standardChrom} #keep standard chrom only
# sed -ie 's/^chr//' ${tmp_standardChrom} #remove `chr` prefix in naming 


#use directories instead
tail -n +2 ${BigFile} > ${dest}/${baseName}.tmp_noheader.txt #remove header file
awk '$2 !~ /_|E|M/' ${dest}/${baseName}.tmp_noheader.txt > ${dest}/${baseName}.tmp_standardChrom.txt #keep standard chrom only
#sed -ie 's/^chr//' ${dest}/${baseName}.tmp_standardChrom.txt #remove `chr` prefix in naming 


echo "spliting ${BigFile} to ${FILE}"
awk -v fileName=$FILE 'BEGIN { FS = "\t" }; {print > fileName"."$2".per_read_modified_base_calls.txt"}' ${dest}/${baseName}.tmp_standardChrom.txt
# (printf "read_id\tchrm\tstrand\tpos\tmod_log_prob\tcan_log_prob\tmod_base\n" && cat $i) > $i.tmp && mv $i.tmp $i

rm ${dest}/${baseName}.tmp_noheader.txt
rm ${dest}/${baseName}.tmp_standardChrom.txt

#(printf "read_id\tchrm\tstrand\tpos\tmod_log_prob\tcan_log_prob\tmod_base\n" && cat $i) > $i.tmp && mv $i.tmp $i