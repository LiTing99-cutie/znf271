#!/usr/bin/sh

################################################
#File Name: scripts/Step_3_RNA_seq.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 11:10:22 PM CST
################################################

set -eo pipefail
script_path=/home/user/data2/lit/project/ZNF271/02-APA-1/
## GTF
tail -n +2 output/pro_dis_bin.txt >  output/pro_dis_bin.nh.txt
bash $script_path/bin/proximal_distal_gtf.sh output/pro_dis_bin.nh.txt
## stringtie
bash -c "time bash $script_path/bin/stringtie.universal.sh \
output/pro_dis_bin.nh.gtf \
output/stringtie \
frag_filter"
echo "stringtie done"
## Wilcoxon test 
echo start: $(date)
Rscript scripts/wilcox_test.R
echo end: $(date)
echo "wilcox done"