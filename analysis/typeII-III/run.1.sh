#!/usr/bin/sh

# reorganize

# ln
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt typeIIAndIII.txt
# PAS read count
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/h/PAusage.bed6+ PAusage.bed6+

# type_II_III list
tail -n +2 output/final_list/typeIIAndIII.txt | cut -f 2 -d ',' > gene.lst

grep -f gene.lst -w PAusage.bed6+ > PAusage.fil.bed6+

# proximal end
cat <(less typeIIAndIII.txt | awk -v OFS='\t' -F',' '$5=="'+'"{print $4,$7-1,$7,".",".",$5}'|\
bedtools intersect -a PAusage.fil.bed6+ -b - -wa -s) \
<(less typeIIAndIII.txt | awk -v OFS='\t' -F',' '$5=="'-'"{print $4,$6-1,$6,".",".",$5}'|\
bedtools intersect -a PAusage.fil.bed6+ -b - -wa -s) > p_e_r_c.txt

# intersect with public annotation (39 genes are validated)
bedtools intersect -a p_e_r_c.txt -b ../supple/benchmark/atlas.upDown30.pos.comChr.bed -s -wa -wb | awk '{print $0,$1":"$10-30":"$6}' | \
awk '{gsub(/chr/,"",$14);print $0}' > 3end_validated_typeII_III.txt
join -t $'\t' -1 14 -2 4 <(sort -k 14,14 3end_validated_typeII_III.txt) \
<(sort -k 4,4 ../supple/benchmark/atlas.clusters.2.0.GRCh38.96.bed)|cut -f 1,2-8,20-24 > 3end_validated_typeII_III.merge.txt
echo -e "Chr\tStart\tEnd\tGene_name\tRead_counts\tStrand\tPA_usage\tPercentage_of_samples\tNumber_of_protocols\tAverage_expression\tCluster_annotation\tPolyA_signals"

# pas_supported_sample
# pas_supported_sample

# merge in R with a high-confidence type II or type III PAS