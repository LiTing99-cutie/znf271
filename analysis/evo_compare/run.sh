#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 06 Mar 2023 02:54:16 PM CST
################################################

set -eou pipefail

# generate gene_name/per/cnt
for spe in h m r;do
bash test.sh $spe
Rscript sknewness.R /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/$spe/$spe.cnt.bed /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/$spe/sknewness.txt
done

[ -d log ] || mkdir log
Rscript pairwise.uni.R /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h/h.cnt.bed \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m/m.cnt.bed \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hm.txt \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hm_diff.h_g.txt >log/hm.log

Rscript pairwise.uni.R /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h/h.cnt.bed \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/r/r.cnt.bed \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hr.txt \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hr_diff.h_g.txt >log/hr.log

Rscript pairwise.uni.R /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/r/r.cnt.bed \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m/m.cnt.bed \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/rm.txt \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/rm_diff.h_g.txt >log/rm.log

# 23/3/27 summarize all 