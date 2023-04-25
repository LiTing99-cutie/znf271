#!/usr/bin/sh

################################################
#File Name: scripts/Step_1_te_backbone.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 09:41:17 PM CST
################################################

set -eo pipefail

project_path=/home/user/data2/lit/project/ZNF271/02-APA-1/

pushd $project_path/
gtfToGenePred -genePredExt -geneNameAsName2 annotation/gencode.v41.basic.annotation.gtf.gz annotation/gencode.v41.basic.annotation.gpe
less annotation/gencode.v41.basic.annotation.gpe | \
awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12,$1}
else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12,$1}}' | \
sort -k4,4 > annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt
popd

Rscript scripts/dominant_te.R
bash scripts/all_exon.sh
Rscript scripts/dominant_te_inter_other_trans_other_gene_fil.R