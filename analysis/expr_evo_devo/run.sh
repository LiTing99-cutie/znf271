#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 02 Mar 2023 11:29:08 AM CST
################################################

set -eou pipefail

# rpkm for each species

## gtf
mkdir -p anno/gtf
ln -s /home/user/data/lit/database/in_house/rheMac10Plus/rheMac10Plus.addgeneName.gtf anno/gtf/rhesus.gtf
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.gtf anno/gtf/mouse.gtf
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gtf anno/gtf/human.gtf

## stringtie loop
mkdir log
nohup bash scripts/stringtie.sh &> log/stringtie.log &

# rhesus gtf shouled be RheMac10Plus not rheMac8
# mouse gtf should be unzipped

nohup bash scripts/stringtie.rerun.sh &> log/stringtie.rerun.log &

nohup bash -c "time Rscript one_way_anova_test.R" &> anova.log &

# only test on evo_list 14942 genes
# input "output/anova.rds"
nohup bash -c "time Rscript devo.R" &> devo.log &
# 100 genes 
# real    0m29.984s
# user    1m37.309s
# sys     0m13.817s

nohup bash -c "time Rscript devo_parallel.R" &> devo_parallel.log &
# 100 genes
# real    3m57.025s
# user    2m57.788s
# sys     0m35.879s