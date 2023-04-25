#!/usr/bin/sh

################################################
#File Name: scripts/Step_1_te_backbone.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 09:41:17 PM CST
################################################

set -eo pipefail

gpe=$1
gtf=$2
output_path=output
do_te_per=output/do_te_per.bed
do_te=output/do_te.bed
exon=output/exon.bed6
echo "extract dominant te"
less $gpe | \
awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12,$1}
else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12,$1}}' | \
sort -k4,4 > $output_path/terminal_exon.txt
Rscript ../scripts/dominant_te.R "$output_path/terminal_exon.txt" "$output_path"
echo "extract all exons"
bash scripts/all_exon.sh $gtf $output_path
echo "exclude te overlapped with other exons"
Rscript ../scripts/dominant_te_inter_other_trans_other_gene_fil.R "$do_te_per" "$do_te" "$exon" "$output_path"