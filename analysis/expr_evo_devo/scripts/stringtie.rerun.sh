#!/usr/bin/sh

################################################
#File Name: test.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 02 Mar 2023 02:25:40 PM CST
################################################

set -eou pipefail

for spe in rhesus mouse;do
[ -d rpkm/$spe ] || mkdir -p rpkm/$spe
annotation=anno/gtf/$spe.gtf
output_dir=rpkm/$spe
for bam in $(ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe/*bam);do
sample=$(basename $bam|sed 's/.sorted.*bam//g')
echo "$spe:$sample"
stringtie -e -p 80 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
done
done




for spe in human rhesus mouse;do
output_dir=rpkm/$spe
[ -f $output_dir/stringtie.rpkm.txt ] && rm -rf $output_dir/stringtie.rpkm.txt
for tab in $(ls $output_dir/stringtie.*tab);do
sample=$(basename $tab|sed -r -e 's/stringtie.E-MTAB-[[:digit:]]{4}.//g' -e 's/.tab//g')
cat $tab | tail -n +2|awk -v OFS='\t' '{print $2,$8,"'$sample'"}' >> $output_dir/stringtie.rpkm.txt
done
done