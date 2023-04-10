#!/usr/bin/sh

################################################
#File Name: bin/merge_stringtie_tab.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 25 Mar 2023 01:16:01 PM CST
################################################

set -eou pipefail

output_dir=$1
[ -f $output_dir/stringtie.rpkm.txt ] && rm -rf $output_dir/stringtie.rpkm.txt
for tab in $(ls $output_dir/stringtie.E-MTAB-6814.*.tab);do
sample=$(basename $tab|sed -e 's/stringtie.E-MTAB-6814.//g' -e 's/.tab//g')
cat $tab | \
egrep 'distal|proximal' | awk -v OFS='\t' '{print $1,$2,$8,"'$sample'"}' >> $output_dir/stringtie.rpkm.txt
done