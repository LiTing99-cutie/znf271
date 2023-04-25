#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 25 Apr 2023 08:10:51 PM CST
################################################

set -eo pipefail

source ./config.txt

[ -d log ] || mkdir log
[ -d output ] || mdkir output

echo "Step_1"
bash scripts/Step_1_te_backbone.sh $gpe $gtf

echo "Step_1_1"
bash scripts/Step_1_1_cutoff_filter_compare.sh PAusage_bed/PAusage.bed6+ output PAusage_bed no

echo "Step_2"
bash scripts/Step_2_pas_map_to_backbone.sh PAusage_bed/PAusage.fil.bed6+ output/do_te_fil.bed output

echo "Step_3"
bash scripts/Step_3_RNA_seq.sh $md $frag $uniq_bam_path &>log/Step_3_RNA_seq.log

echo "Step_3_1"
time bash scripts/Step_3_1_disrupt_cds.sh $gtf $gpe $archive $species no &> log/Step_3_1_disrupt_cds.log

echo "Step_4"
bash scripts/Step_4_case.sh $md $frag $case