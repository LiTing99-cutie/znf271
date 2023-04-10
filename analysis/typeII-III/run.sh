#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 01 Feb 2023 07:54:09 PM CST
################################################

set -eou pipefail

# type_II_III list
tail -n +2 output/final_list/typeIIAndIII.txt | cut -f 2 -d ',' > gene.lst

# PAS read count
/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/h/PAusage.bed6+


grep -f gene.lst -w /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/h/PAusage.bed6+ > PAusage.fil.bed6+


# proximal end
cat <(less output/final_list/typeIIAndIII.txt | awk -v OFS='\t' -F',' '$5=="'+'"{print $4,$7-1,$7,".",".",$5}'|\
bedtools intersect -a PAusage.fil.bed6+ -b - -wa -s) \
<(less output/final_list/typeIIAndIII.txt | awk -v OFS='\t' -F',' '$5=="'-'"{print $4,$6-1,$6,".",".",$5}'|\
bedtools intersect -a PAusage.fil.bed6+ -b - -wa -s) > p_e_r_c.txt

# distal end
cat <(less output/final_list/typeIIAndIII.txt | awk -v OFS='\t' -F',' '$5=="'+'"{print $4,$9-1,$9,".",".",$5}'|\
bedtools intersect -a PAusage.fil.bed6+ -b - -wa -s) \
<(less output/final_list/typeIIAndIII.txt | awk -v OFS='\t' -F',' '$5=="'-'"{print $4,$8-1,$8,".",".",$5}'|\
bedtools intersect -a PAusage.fil.bed6+ -b - -wa -s) > d_e_r_c.txt

# paste
paste -d '\t' <(cut -f 4-5 p_e_r_c.txt) <(cut -f 5 d_e_r_c.txt) > g_p_d_r_c.txt

res=/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.with_geneName.txt

join -1 1 -2 3 <(sort -k1,1 g_p_d_r_c.txt) <(sort -k3,3 $res) | tr ' ' '\t'> g_p_d_r_c_pau_expr.txt

join -t $'\t' -1 1 -2 2 <(sort -k1,1 g_p_d_r_c_pau_expr.txt) \
<(less /home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt|tr , '\t'|sort -k2,2)> g_p_d_r_c_pau_expr_type.tx

awk -F='\t' -v OFS='\t' '{print $1,$5,$2,$3,$9,$23}' g_p_d_r_c_pau_expr_type.txt |sed 's/non_coding/no_protein/g;s/another protein product/alternative_protein/g' \
> supple.nh.txt

cat <(echo -e "Gene Name\tGene Type\tReads Count of Proximal PAS\tReads Count of Distal PAS\tDelta PAU\tProtein Corresponding To Proximal PAS") supple.nh.txt \
> supple.h.txt