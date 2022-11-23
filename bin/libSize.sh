#!/usr/bin/sh

################################################
#File Name: bin/lib.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 05 Oct 2022 02:15:46 PM CST
################################################

set -eou pipefail

# absolute path of bams (uniq bams) *uniq.sorted.bam
bam_lst=$1
# parallel job number
job_n=$2
# parallel log
log=$3
# output dir 
output_dir=$4
# pair-end or single-end
end=$5

export output_dir=$output_dir
[ -d $output_dir ] || mkdir -p $output_dir
[ -f $output_dir/library_size.txt ] && rm -rf $output_dir/library_size.txt

# bam_stat all the bams
cat $bam_lst | parallel \
-j $job_n --plus --joblog $log \
'sample={/...} ; echo $sample;bam_stat.py -q 60 -i {} > $output_dir/$sample.txt'

if [ $end = "se" ];then
    for file in $(ls $output_dir/*.txt);
    do 
    sample=$(basename $file | sed 's/.txt//g' )
    size=$(grep "Total" $file | awk '{print $NF}')
    echo -e "$sample\t$size" >> $output_dir/library_size.txt
    done
else
    for file in $(ls $output_dir/*.txt);
    do 
    sample=$(basename $file | sed 's/.txt//g' )
    size=$(grep "proper pairs" $file | awk '{print $NF/2}')
    echo -e "$sample\t$size" >> $output_dir/library_size.txt
    done
fi