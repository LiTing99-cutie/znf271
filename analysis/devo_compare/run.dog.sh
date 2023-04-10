#!/usr/bin/sh

################################################
#File Name: run.dog.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 13 Mar 2023 08:09:14 PM CST
################################################

set -eou pipefail



annotation=./dog_stringtie_anno/anno.gtf
spe=dog
output_dir=./$spe/analysis/stringtie_whole_genome_ref_based_all_1_loose
echo $annotation
echo $output_dir
[ -d $output_dir ] || mkdir -p $output_dir

for bam in $(ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe/*bam);do
    sample=$(basename -s ".uniq.sorted.bam" $bam)
    echo $sample
    stringtie -p 60 -G $annotation -o $output_dir/stringtie.$sample.gtf -A $output_dir/stringtie.$sample.tab --rf $bam
done

[ -f $output_dir/stringtie.rpkm.txt ] && rm -rf $output_dir/stringtie.rpkm.txt

for tab in $(ls $output_dir/stringtie.SRR*tab);do
sample=$(basename $tab|cut -f 2 -d '.')
cat $tab | \
egrep 'distal|proximal' | awk -v OFS='\t' '{print $1,$2,$8,"'$sample'"}' >> $output_dir/stringtie.rpkm.txt
done

Rscript script/develop_case.R dog ENSCAFG00845006933

cp script/wilcox_test.R script/wilcox_test_dog.R

Rscript script/wilcox_test_dog.R dog ENSCAFG00845006933