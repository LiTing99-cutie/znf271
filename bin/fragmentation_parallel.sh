#!/usr/bin/sh

################################################
#File Name: fragmentation_parallel.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 03 Oct 2022 09:52:35 PM CST
################################################

set -eou pipefail

# absolute path of bams (uniq bams)
bam_lst=$1
# parallel job number
job_n=$2
# parallel log
log=$3
# output dir 
output_dir=$4

export output_dir=$output_dir

cat $bam_lst | parallel \
-j $job_n --plus --joblog $log \
'sample={/.};
python2 /home/user/data2/lit/project/ZNF271/src/fragmentation.py \
--bam {} \
--gpe /home/user/data2/lit/project/ZNF271/02-APA/annotation/gencode.v41.annotation.gpe > $output_dir/$sample.frag.txt' 



