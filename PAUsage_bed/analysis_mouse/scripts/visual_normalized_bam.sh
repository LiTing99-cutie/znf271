#!/usr/bin/sh

################################################
#File Name: scripts/visual_normalized_bam.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 24 Apr 2023 04:11:17 PM CST
################################################

set -eo pipefail

# variable
bamCoverage=/home/user/BGM/lit/anaconda3/envs/devtools/bin/bamCoverage
bam_path=$1
output_path=$2
output_bam_path=$3
metadata_path=$4
[ -d $output_bam_path ] || mkdir -p $output_bam_path
[ -d $output_path ] || mkdir -p $output_path
# for two stage
for stage in embryo postnatal;do
SN=$metadata_path/$stage.sampleName.txt
ls $bam_path/*bam | egrep -f $SN >$output_path/$stage.lst
echo "samtools merge"
time samtools merge -@ 80 -f -o $output_bam_path/$stage.bam -b $output_path/$stage.lst 
echo "samtools index"
time samtools index -@ 80 $output_bam_path/$stage.bam
echo "bamCoverage"
time $bamCoverage --normalizeUsing RPKM -p 80 -b $output_bam_path/$stage.bam -o $output_path/$stage.bw -bs 100
done