#!/usr/bin/sh

################################################
#File Name: scripts/te_pro_dis_bed.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 20 Apr 2023 08:24:18 PM CST
################################################

set -eou pipefail

te=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/pro_dis_bin.nh.txt

less $te | awk -v OFS='\t' '{print $2,$4,$5,$1,".",$3}' > output/pro.bed
less $te | awk -v OFS='\t' '{print $2,$6,$7,$1,".",$3}' > output/dis.bed
awk -v OFS='\t' '{if($3=="+"){print $2,$4,$7,$1,".",$3}else{print $2,$6,$5,$1,".",$3}}' $te> output/te.bed