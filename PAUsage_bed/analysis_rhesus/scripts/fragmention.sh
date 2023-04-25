#!/usr/bin/sh

################################################
#File Name: scripts/fragmention.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 20 Apr 2023 09:29:13 PM CST
################################################

set -eo pipefail

# uniq bam file path
bam_path=$1
# threads
thread=$2
# prefix
prefix=$3
# output_path
output_path=$4
# gpe
gpe=$5

[ -d $output_path/tmp_frag/ ] || mkdir -p $output_path/tmp_frag/
script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/fragmentation_parallel.sh

# generate bam_list and index if necessary
ls $bam_path/*bam | egrep -i '[^/]brain[^/]|[^/]cerebellum[^/]'> $output_path/tmp_frag/bam.lst
bai_number=$(find $bam_path -name *bai | wc -l)
bam_number=$(less $output_path/tmp_frag/bam.lst|wc -l)
if [ $bai_number -lt $bam_number ];then
cat $output_path/tmp_frag/bam.lst | parallel --joblog log/index.$prefix.log -j 20 samtools index {}
fi
echo "index done"

bash -c "time bash $script $output_path/tmp_frag/bam.lst $thread log/frag.$prefix.log $output_path/tmp_frag $gpe"

for frag in $(ls $output_path/tmp_frag/*frag.txt);do
perl /home/user/data2/lit/project/ZNF271/src/calRnaSeqRbScore.modi.pl --fragmentation $frag
done | cat > $output_path/fragmentation.score.txt

less $output_path/fragmentation.score.txt | sed -r 's/.*E-MTAB-[0-9]+.//g;s/.sorted.uniq.frag.txt//g' > $output_path/fragmentation.score.clean.txt
