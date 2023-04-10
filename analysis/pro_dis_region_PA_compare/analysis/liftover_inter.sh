#!/usr/bin/sh

################################################
#File Name: liftover_inter.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 03 Mar 2023 02:06:51 PM CST
################################################

set -eou pipefail

# PA usage bed
mkdir iso_seq
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/m/PAusage.bed6+ iso_seq/Mm10.bed6+
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/r/PAusage.bed6+ iso_seq/RheMac8.bed6+
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/h/PAusage.bed6+ iso_seq/hg38.bed6+

# human terminal exons 7468
te=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.txt
# here -1 to make sure everything is right
# but before we do not do this? because the terminal exon annotation is based on loose_for_m_r_h
# manually checked
less $te | \
awk -v OFS='\t' '{if($3=="+"){print $2,$4,$5,$1,".",$3}else{print $2,$4-1,$5,$1,".",$3}}' > hg38.pro.bed
less $te | \
awk -v OFS='\t' '{if($3=="+"){print $2,$6,$7,$1,".",$3}else{print $2,$6-1,$7-1,$1,".",$3}}' > hg38.dis.bed

# liftover and intersect 
liftOver_path=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/liftover
for type in pro dis;do
for spe in Mm10 RheMac8;do
liftOver -minMatch=0.5 hg38.$type.bed \
$liftOver_path/hg38/hg38To$spe.over.chain.gz \
hg38.$type.To$spe.bed hg38.$type.To$spe.unmapped
bedtools intersect -a hg38.$type.To$spe.bed -b <(grep -v ERCC iso_seq/$spe.bed6+|awk '$5>1') -wa -wb> hg38.$type.$spe.bed6+
done
done

# intersect for human only 
for type in pro dis;do
bedtools intersect -a hg38.$type.bed -b <(grep -v ERCC iso_seq/hg38.bed6+|awk '$5>1') -wa -wb | awk '$4==$10'> hg38.$type.bed6+
done

# terminal exons in human cannot liftover
for spe in Mm10 RheMac8;do
cat <(less hg38.dis.To$spe.unmapped | grep -v "^#" | cut -f 4) \
<(less hg38.pro.To$spe.unmapped | grep -v "^#" | cut -f 4) |sort|uniq > hg38.To$spe.unmapped.gene.lst
done

cat hg38.ToRheMac8.unmapped.gene.lst hg38.ToMm10.unmapped.gene.lst|sort|uniq > hg38.unmapped.gene.lst
