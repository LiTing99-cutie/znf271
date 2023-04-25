#!/usr/bin/sh

################################################
#File Name: scripts/Step_1_1_cutoff_filter_compare.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 11:32:33 PM CST
################################################

set -eo pipefail


## select cutoff 
Rscript scripts/cutoff_pas.R

## filter
Rscript scripts/fil.PAS.R

## compare with 3'-seq and dapars2
bash scripts/compare_3_seq_dapars2.sh
Rscript scripts/compare_3_seq_dapars2.R