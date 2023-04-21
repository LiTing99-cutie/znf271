#!/usr/bin/sh

################################################
#File Name: scripts/fragmention.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 20 Apr 2023 09:29:13 PM CST
################################################
set -eo pipefail

output_path=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output_m
for frag in $(ls $output_path/tmp_frag/*frag.txt);do
perl /home/user/data2/lit/project/ZNF271/src/calRnaSeqRbScore.modi.pl --fragmentation $frag
done | cat > $output_path/fragmentation.score.txt

echo "all done"