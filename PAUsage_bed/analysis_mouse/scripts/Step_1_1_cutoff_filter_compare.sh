#!/usr/bin/sh

################################################
#File Name: scripts/Step_1_1_cutoff_filter_compare.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 11:32:33 PM CST
################################################

set -eo pipefail

PAUsage_bed=$1
output_path=$2
PAUsage_fil_bed_path=$3
compare=$4
## select cutoff 
Rscript ../scripts/cutoff_pas.R "$PAUsage_bed" "$output_path"

## filter
Rscript ../scripts/fil.PAS.R "$PAUsage_bed" 0 0 2 "$PAUsage_fil_bed_path"

## compare with 3'-seq and dapars2
if [ $compare == "yes" ];then
	bash scripts/compare_3_seq_dapars2.sh
	Rscript scripts/compare_3_seq_dapars2.R
fi