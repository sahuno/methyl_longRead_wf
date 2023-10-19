(base) bash-4.2$ cat results/build_consensus/consensus_HighQual_motif_sites.bed | cut -f 1| sort | uniq -c | sort -nr
2996084 chr2
2899691 chr1
2646097 chr5
2623161 chr4
2420289 chr7
2327420 chr6
2323884 chr3
2309669 chr11
2214641 chr8
2184775 chr10
2116031 chr9
1925251 chr13
1884879 chr12
1777787 chr14
1728712 chr15
1686730 chr17
1480446 chr16
1419712 chr18
1091594 chr19
1002338 chrX
   2701 chrY
    572 chrM
(base) bash-4.2$ 

cat results/build_consensus/consensus_HighQual_motif_sites.bed | cut -f 4| sort | uniq -c | sort -nr


bedtools intersect -a results/splitChrom/D-A-3_modpileup_5mC/data/D-A-3_modpileup_5mC.chr13.bed -b results/build_consensus/consensus_HighQual_motif_sites.bed > results/filter_consensus/D-A-3_modpileup_5mC/data/D-A-3_modpileup_5mC.chr13.filtered_Consensus.bed
bedtools intersect -a results/sortBed/D-A-3_modpileup_5mC/data/D-A-3_modpileup_5mC.chr13.sorted.bed -b results/build_consensus/consensus_HighQual_motif_sites.bed > results/filter_consensus/D-A-3_modpileup_5mC/data/D-A-3_modpileup_5mC.chr13.filtered_Consensus.bed



sort -k1,1 -k2,2n results/splitChrom/D-A-3_modpileup_5mC/data/D-A-3_modpileup_5mC.chr13.bed > sandbox/test_sort_D-A-3_modpileup_5mC.chr13.bed


bedtools intersect -a sandbox/test_sort_D-A-3_modpileup_5mC.chr13.bed -b results/build_consensus/consensus_HighQual_motif_sites.bed


grep "" results/sortBed/D-0-1_modpileup_5mC/data/D-0-1_modpileup_5mC.chr18.sorted.bed