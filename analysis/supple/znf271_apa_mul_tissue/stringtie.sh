#!/usr/bin/sh

################################################
#File Name: stringtie.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 25 Mar 2023 02:28:48 PM CST
################################################

set -eou pipefail

bam_lst=$1
annotation=$2
output_dir=$3

time cat bam_lst/human.uniq_bam.lst|while read bam;do
	sample=$(basename $bam|sed 's/.sorted.uniq.bam//g')
	echo $sample
	stringtie -p 60 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
done
