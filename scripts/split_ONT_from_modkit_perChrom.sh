#!/bin/sh

## # BSUB -J splitChrom
## BSUB -q general
# #BSUB -gpu num=2
## BSUB -n 4
## BSUB -R rusage[mem=8]
## BSUB -W 0:30
## BSUB -e splitChrom%J_%I.err
## # BSUB -o splitChrom%J_%I.out

#samuel ahuno
#tested this works
baseName=$1 #sampel name
Big_bed_file=$2 #input file
dest=$3 #output directory

#test
# baseName=D-0-1
# # BigFile=/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-A-3/D-A-3_modpileup_5mC.bed
# Big_bed_file=/juno/work/greenbaum/users/ahunos/sandbox/results_modkit/results/modkit/D-0-1/D-0-1_modpileup_5mC.bed
# dest=/juno/work/greenbaum/users/ahunos/apps/methyl_longRead_wf/sandbox

#set dir
DIR_SAVE=${dest}/${baseName}
echo $DIR_SAVE
echo $BigFile
echo $dest


#remove all non-standard chromosomes
# awk '$1 !~ /_|E|M/' ${BigFile} > ${dest}/${baseName}.tmp_standardChrom.txt #keep standard chrom only
# awk -v fileName=$FILE 'BEGIN { FS = "\t" }; {print > fileName"."$1".bed"}' ${dest}/${baseName}.tmp_standardChrom.txt

echo "spliting ${Big_bed_file} to chromosomal files"
#keep standard chrom only
awk '$1 !~ /_|E|M/' ${Big_bed_file} | awk -v fileName=$DIR_SAVE 'BEGIN { FS = "\t" }; {print > fileName"."$1".bed"}'

# (printf "read_id\tchrm\tstrand\tpos\tmod_log_prob\tcan_log_prob\tmod_base\n" && cat $i) > $i.tmp && mv $i.tmp $i

# echo "removing tmp files"
# rm ${dest}/${baseName}.tmp_noheader.txt
#rm ${dest}/${baseName}.tmp_standardChrom.txt
# echo "done removing tmp files"

echo "done splitting ${Big_bed_file} to chromosomal files"