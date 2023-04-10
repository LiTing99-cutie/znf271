#!/usr/bin/sh

################################################
#File Name: stringtie.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 10 Oct 2022 04:51:05 PM CST
################################################

set -eou pipefail
annotation=$1
output_dir=$2
frag_filter_or_not=$3
spe=$4
echo $annotation
echo $output_dir
[ -d $output_dir ] || mkdir -p $output_dir


for bam in $(ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe/*bam);do
    sample=$(basename $bam|sed 's/.sorted.*bam//g')
    echo $sample
    stringtie -p 60 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
done

[ -f $output_dir/stringtie.rpkm.txt ] && rm -rf $output_dir/stringtie.rpkm.txt

for tab in $(ls $output_dir/stringtie.E-MTAB*tab);do
sample=$(basename $tab|sed -r -e 's/stringtie.E-MTAB-[[:digit:]]{4}.//g' -e 's/.tab//g')
cat $tab | \
egrep 'distal|proximal' | awk -v OFS='\t' '{print $1,$2,$8,"'$sample'"}' >> $output_dir/stringtie.rpkm.txt
done