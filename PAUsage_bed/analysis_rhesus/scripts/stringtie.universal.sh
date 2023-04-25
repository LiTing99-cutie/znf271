#!/usr/bin/sh

################################################
#File Name: stringtie.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 10 Oct 2022 04:51:05 PM CST
################################################

set -eo pipefail
annotation=$1
output_dir=$2
frag_filter=$3
cutoff=$4
uniq_bam_path=$5
frag_score=$6

echo $annotation
echo $output_dir
echo "whether calculate only on filtered bams: $frag_filter"
echo "fragmentation score cutoff: $cutoff"
echo $uniq_bam_path
echo $frag_score

[ -d $output_dir ] || mkdir -p $output_dir

# stringite call region rpkm
source activate base
if [ $frag_filter == "yes" ];then
    for bam in $(ls $uniq_bam_path/*bam|grep -f $output_dir/fil.sample.lst);do
        # generate filtered sample list
        less $frag_score |awk '$2>"'$cutoff'"'| cut -f 1 > $output_dir/fil.sample.lst
        sample=$(basename $bam|sed 's/.sorted.uniq.sorted.bam//g')
        echo $sample
        stringtie -p 60 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
    done
else
    for bam in $(ls $uniq_bam_path/*bam);do
        sample=$(basename $bam|sed -r 's/.sorted.uniq.sorted.bam//g;s/.uniq.sorted.bam//g;s/.sorted.uniq.bam//g;s/E-MTAB-[0-9]+.//g')
        echo $sample
        stringtie -p 60 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
    done
fi 


[ -f $output_dir/stringtie.rpkm.txt ] && rm -rf $output_dir/stringtie.rpkm.txt

for tab in $(ls $output_dir/*tab);do
sample=$(basename $tab|sed -e 's/stringtie.//g' -e 's/.tab//g')
cat $tab | \
egrep 'distal|proximal' | awk -v OFS='\t' '{print $1,$2,$8,"'$sample'"}' >> $output_dir/stringtie.rpkm.txt
done
