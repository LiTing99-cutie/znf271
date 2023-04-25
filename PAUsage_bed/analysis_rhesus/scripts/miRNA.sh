#!/usr/bin/sh

################################################
#File Name: scripts/miRNA.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 20 Apr 2023 08:37:32 PM CST
################################################

set -eou pipefail

[ -d output/miRNA ] || mkdir -p output/miRNA 
pushd output/miRNA 
bed_path=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output
miRNA_site=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/miRNA/Predicted_Target_Locations.default_predictions.hg38.bed

# APA and miRNA-binding sites correlation
for type in te pro dis;do
bedtools intersect -a $bed_path/$type.bed -b $miRNA_site -wo | awk -v OFS='\t' '$20=$19/($9-$8)' | awk '$20>0.5'> hg38.$type.miRNA.txt
done

