#!/usr/bin/sh

################################################
#File Name: uniq.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 03 Oct 2022 06:41:02 PM CST
################################################

set -eou pipefail

# for sorted bams

# absolute path of bams (mapped bams)
bam_lst=$1
# thread number
thread_n=$2
# output dir of uniq bams (absolute)
output_dir=$3

for bam in $(cat $bam_lst);do
sample=$(basename $bam|sed 's/.bam//g')
echo "extract uniq reads of ${sample}..."
# samtools view -@ $thread_n -q 60 -bhu $bam | samtools sort -@ $thread_n -o $output_dir/$sample.uniq.sorted.bam -
# for m and r uniq,  no need to sort
samtools view -@ $thread_n -q 60 -bhu -o $output_dir/$sample.uniq.bam $bam
done