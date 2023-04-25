#!/usr/bin/sh

################################################
#File Name: scripts/Step_3_RNA_seq.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 11:10:22 PM CST
################################################

set -eo pipefail

meta_c=$1
frag_score=$2
uniq_bam_path=$3
stringtie_script=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/analysis_mouse/scripts/stringtie.universal.sh
rpkm_file=output/stringtie/stringtie.rpkm.txt
output_path=output
## GTF
tail -n +2 output/pro_dis_bin.txt >  output/pro_dis_bin.nh.txt
bash /home/user/data2/lit/project/ZNF271/02-APA-1/bin/proximal_distal_gtf.sh output/pro_dis_bin.nh.txt
## stringtie
bash -c "time bash $stringtie_script \
output/pro_dis_bin.nh.gtf \
output/stringtie \
no \
0.885 \
$uniq_bam_path \
$frag_score"
echo "stringtie done"
## Wilcoxon test 
echo start: $(date)
Rscript ../scripts/wilcox_test.R "$rpkm_file" "$meta_c" "$frag_score" "$output_path"
echo end: $(date)
echo "wilcox done"