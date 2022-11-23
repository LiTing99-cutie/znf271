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
echo $annotation
echo $output_dir
[ -d $output_dir ] || mkdir -p $output_dir

if [ $frag_filter_or_not ];then
    for bam in $(cat /home/user/data2/lit/project/ZNF271/02-APA/data/lst/frag_fil_develop_bam.lst);do
        sample=$(basename $bam|sed 's/.sorted.uniq.sorted.bam//g')
        echo $sample
        stringtie -p 60 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
    done
else
    for bam in $(ls /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/*bam);do
        sample=$(basename $bam|sed 's/.sorted.uniq.sorted.bam//g')
        echo $sample
        stringtie -p 60 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
    done
fi 


[ -f $output_dir/stringtie.rpkm.txt ] && rm -rf $output_dir/stringtie.rpkm.txt

for tab in $(ls $output_dir/stringtie.E-MTAB-6814.*.tab);do
sample=$(basename $tab|sed -e 's/stringtie.E-MTAB-6814.//g' -e 's/.tab//g')
cat $tab | \
egrep 'distal|proximal' | awk -v OFS='\t' '{print $1,$2,$8,"'$sample'"}' >> $output_dir/stringtie.rpkm.txt
done