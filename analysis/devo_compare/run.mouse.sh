#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 14 Feb 2023 09:59:25 PM CST
################################################

set -eou pipefail

# mkdir 
mkdir -p mouse/annotation/terminal_exon/
mkdir -p mouse/iso_seq
mkdir script
mkdir log
# variable
gtf=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.gtf.gz
universal_script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/ref_based_develop_diff_pa_usage.universal.sh
terminal_exon_script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/terminal_exon_annotation.loose.py
# script
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/PAusage.m.c_anno.bed6+ mouse/iso_seq/PAusage.bed6+
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/bin/proximal_distal_gtf.sh script/
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/bin/stringtie.universal.forOtherSpecies.sh script/stringtie.universal.sh
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/bin/wilcox_test.forOtherSpecies.R script/wilcox_test.R
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/bin/ref_based_develop_diff_pa_usage.universal.forOtherSpecies.sh script/universal.sh
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/bin/develop_case.forOtherSpecies.R script/develop_case.R

# 1. get terminal exons
    # extract chr,strand,terminal exon start,terminal exon end,ensembl gene id
    gtfToGenePred -genePredExt -geneNameAsName2 $gtf mouse/annotation/annotation.gpe
    less mouse/annotation/annotation.gpe | \
    awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12}
    else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12}}' | \
    sort -k4,4 > mouse/annotation/terminal_exon.genePredExt

echo "Generate annotation and call rpkm of common and extended region" 
nohup bash -c "time bash script/universal.sh \
ref_based_all_1_loose \
mouse/iso_seq/PAusage.bed6+ \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/devo_compare \
$terminal_exon_script \
mouse/annotation/terminal_exon.genePredExt \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/devo_compare/script \
mouse" &> log/mouse.log &