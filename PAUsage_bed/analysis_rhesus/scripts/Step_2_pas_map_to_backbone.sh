#!/usr/bin/sh

################################################
#File Name: scripts/Step_2_pas_map_to_backbone.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 09:50:40 PM CST
################################################

set -eo pipefail

PAUsage_fil_bed=$1
do_te=$2
output_path=$3
Rscript ../scripts/pas_map.R "$PAUsage_fil_bed" "$do_te" "$output_path"