#!/usr/bin/sh

################################################
#File Name: scripts/Step_4_description_function.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 11:29:32 PM CST
################################################

set -eo pipefail

md=$1
frag_score=$2
gene_name=$3

Rscript ../scripts/develop_case.R "$PWD/output/stringtie/stringtie.rpkm.txt" \
"$md" \
"$frag_score" \
"$gene_name" \
"$PWD/output/"