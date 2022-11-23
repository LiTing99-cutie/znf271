#!/usr/bin/sh

################################################
#File Name: coverage.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 28 Sep 2022 02:51:47 PM CST
################################################

set -eou pipefail

gene_name=$1
export gene_name=$1
bam_lst=$2

# zless /home/user/data2/lit/project/ZNF271/02-APA/annotation/gencode.v41.annotation.gtf.gz | grep -i $gene_name | awk '$3=="gene"' | cut -f 1-8  \
# > /home/user/data2/lit/project/ZNF271/02-APA/annotation/target.$gene_name.bed

cat $bam_lst | parallel \
-j 40 --plus --joblog log/target_bam.$gene_name.log \
'sample={/.};bedtools intersect -a {} -b /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/target.$gene_name.bed > \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/$sample.target.$gene_name.bam' 
