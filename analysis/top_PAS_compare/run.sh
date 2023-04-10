#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 07 Feb 2023 07:39:45 PM CST
################################################

set -eou pipefail

mkdir data/
for spe in h r m;do
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/$spe/PAusage.bed6+ data/PAusage.$spe.bed6+
done

# cp from galaxy
# /rd1/user/lit/project/271/PA/analysis/crossSpecies/top_PAS_compare/analysis/compare
# res_hm.rds
# res_hr.rds