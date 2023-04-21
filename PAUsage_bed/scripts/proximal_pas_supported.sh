#!/usr/bin/sh

################################################
#File Name: scripts/proximal_pas_supported.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 18 Apr 2023 09:27:07 PM CST
################################################

set -eou pipefail

disrupt_cds=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/disrupt_cds_pro_pa.txt
PAS_3_atlas=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/benchmark/atlas.pos.bed
PAS_3_atlas_all=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/benchmark/atlas.clusters.2.0.GRCh38.96.bed

[ -d output/proximal_pas_supported ] || mkdir output/proximal_pas_supported && pushd output/proximal_pas_supported

# function
select_com_chr(){
    grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]'
}
updown30(){
awk -v OFS='\t' '{if($2-30>=0)print $1,$2-30,$3+30,$4,$5,$6;else print $1,1,$3,$4,$5,$6}'|sort -k1,1 -k2,2n
}

less $PAS_3_atlas|awk -v OFS='\t' '{print "chr"$1,$2,$3,$4,$5,$6}'|select_com_chr|updown30 > PAS_3_atlas_updown30.bed
less $disrupt_cds | tail -n +2 |awk -v FS='\t' -v OFS='\t' '{print $6,$7,$8,$1,$9,$10,$11}' > p_e_r_c.txt
# intersect with public annotation
bedtools intersect -a p_e_r_c.txt -b PAS_3_atlas_updown30.bed -s -wa -wb | awk -v OFS='\t' '{print $0,$1":"$10-30":"$6}' | \
awk -v OFS='\t' '{gsub(/chr/,"",$14);print $0}' > 3end_validated.txt
join -t $'\t' -1 14 -2 4 <(sort -k 14,14 3end_validated.txt) \
<(sort -k 4,4 $PAS_3_atlas_all)|cut -f 1,2-8,20-24 > 3end_validated.merge.txt
echo -e "ID\tChr\tStart\tEnd\tGene_name\tRead_counts\tStrand\tPA_usage\tPercentage_of_samples\tNumber_of_protocols\tAverage_expression\tCluster_annotation\tPolyA_signals" > header.txt

cat header.txt 3end_validated.merge.txt > 3end_validated.merge.h.txt