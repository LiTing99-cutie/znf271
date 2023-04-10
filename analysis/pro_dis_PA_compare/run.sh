#!/usr/bin/sh

################################################
#File Name: pro_dis_PA_compare/run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 18 Feb 2023 03:26:07 PM CST
################################################

set -eou pipefail

# 1.ortholog
mkdir ortholog
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/human_macaque_ortholog.txt ortholog/human_macaque_ortholog.txt
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/macaque_mouse_ortholog.txt ortholog/macaque_mouse_ortholog.txt
# ensembl 107
# human_mouse_ortholog.txt

cd ortholog
join <(awk -F' ' '$3~/ortholog_one2one/' human_macaque_ortholog.txt|sort -k1,1) \
<(awk -F' ' '$3~/ortholog_one2one/' human_mouse_ortholog.txt|sort -k1,1) | cut -f 1,2,4 -d ' ' > hrm.ortholog.txt

# 2023/3/6
awk -F' ' '$3~/ortholog_one2one/' human_macaque_ortholog.txt|sort -k1,1|cut -f 1,2 > hr.txt
awk -F' ' '$3~/ortholog_one2one/' human_mouse_ortholog.txt|sort -k1,1|cut -f 1,2 > hm.txt
awk -F' ' '$3~/ortholog_one2one/' macaque_mouse_ortholog.txt|sort -k1,1|cut -f 1,2 > rm.txt
# 2.terminal exon and PA file
cd ..
for spe in h m r;do
mkdir $spe
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/$spe/PAusage.bed6+ $spe
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/$spe/terminal_exon.bed $spe
done

for spe in h m r;do
pushd $spe
bedtools intersect -a terminal_exon.bed -b <(grep -v ERCC PAusage.bed6+) -wb -s| awk '$4==$10 && $11>1' | cut -f 7-12 > PAusage.te.bed6+
popd
done

# 3.output
mkdir output