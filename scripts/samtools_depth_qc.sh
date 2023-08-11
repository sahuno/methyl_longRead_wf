#!bin/bash

#to run
#sh samtools_depth_qc.sh '/juno/work/greenbaum/projects/TRI_EPIGENETIC/megalodon/' /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/CpG_sites_0based_mm10.bed
#cd /juno/work/greenbaum/projects/TRI_EPIGENETIC/megalodon/

find $1 -name "mappings.sorted.bam" > megalodon_mapping_bam_files.txt
samtools depth -aa \
-f megalodon_mapping_bam_files.txt \
-b $2 \
-o out_megalodon_mapping_depth_triplicate_mm10.txt 


#plot depth
Rscript scripts/R/plot_samtools_depth_coverage.r --depth_file out_megalodon_mapping_depth_triplicate_mm10.txt  --paths_bams_txt megalodon_mapping_bam_files.txt


#-g SECONDARY,QCFAIL,DUP \

#add header to depth file
# header=`sed 's|.*/megalodon/\(.*\)/mappings.sorted.bam|\1|' megalodon_mapping_bam_files.txt  | tr '\n' '\t' | sed 's/\t$/\n/' | sed 's/^/chr\tpos\t/'`
# echo "$header"
# (echo $header && cat out_megalodon_mapping_depth_triplicate_mm10.txt ) > out_megalodon_mapping_depth_triplicate_mm10.txt1 && mv out_megalodon_mapping_depth_triplicate_mm10.txt1 out_megalodon_mapping_depth_triplicate_mm10.txt 
#sed '1d' file.txt > tmpfile; mv tmpfile file.txt


