#!/usr/bin/sh

################################################
#File Name: test.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 06 Mar 2023 02:40:57 PM CST
################################################

set -eou pipefail

# variables
spe=$1
alias g_c="cut -f4|sort|uniq -c|wc -l"

# softlinks
[ -d $spe ] || mkdir $spe 
[ -f $spe/terminal_exon.bed ] || ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/$spe/terminal_exon.bed $spe
[ -f $spe/PAusage.bed6+ ] || ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/$spe/PAusage.bed6+ $spe

# extend end for 30bp
# 23-3-20 add |sort -k1,1 -k2,2n
cd $spe
awk -v OFS='\t' '{if($6=="+")print $1,$2,$3+30,$4,$5,$6;else print $1,$2-30,$3,$4,$5,$6}' terminal_exon.bed |
awk -v OFS='\t' '{if($2<0)print $1,1,$3,$4,$5,$6;else print $1,$2,$3,$4,$5,$6}'|sort -k1,1 -k2,2n> terminal_exon.ext.bed

# intersect with PAusage.bed6+
# gene_name; PAS localization; cnt
# 23-3-20 add |sort -k1,1 -k2,2n and -nonamecheck
# bedtools intersect -a terminal_exon.ext.bed -b <(grep -v ERCC PAusage.bed6+) -wa -wb -s |awk '$4==$10 && $11>1'|\
bedtools intersect -nonamecheck -a terminal_exon.ext.bed -b <(sort -k1,1 -k2,2n PAusage.bed6+) -wa -wb -s |awk '$4==$10 && $11>1'|\
awk -v OFS='\t' '{if($6=="+")print $4,$9-$2,$11;else print $4,$3-$9,$11}' > $spe.cnt.bed
# 2023-4-7 add
bedtools intersect -nonamecheck -a terminal_exon.ext.bed -b <(sort -k1,1 -k2,2n PAusage.bed6+) -wa -wb -s |awk '$4==$10 && $11>1' > PAusage.te.bed6+