#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 31 Mar 2023 02:45:22 PM CST
################################################

set -eou pipefail


# PASs in terminal exon in each species
# Normalization on coverage to make cutoff on PAS read coverage comparable
mkdir h m r
for spe in h m r ;do
if [ $spe == h ];then
sf=22.2904
elif [ $spe == m ];then
sf=12.7089
else [ $spe == r ]
sf=1
fi
bedtools intersect -nonamecheck -a ../evo_compare/$spe/terminal_exon.ext.bed -b <(sort -k1,1 -k2,2n ../evo_compare/$spe/PAusage.bed6+) -wa -wb -s | tee $spe/tmp.bed|\
awk -v sf=$sf '$13=$11/sf'|awk '$4==$10&&$13>=2'|awk '{print $7,$8,$9,$10":"$9,$13,$12}' > $spe/PAS.te.cnt.scale.bed
Rscript bin/PAusage.R $spe/PAS.te.cnt.scale.bed $spe/PAS.te.PAusage.scale.bed
less $spe/tmp.bed | awk '$4==$10 && $11>1' | awk '{print $7,$8,$9,$10":"$9,$11,$12}' > $spe/PAS.te.cnt.bed
Rscript bin/PAusage.R $spe/PAS.te.cnt.bed $spe/PAS.te.PAusage.bed
done

# Only include those genes with 1:1 h-r orthologs and replace r gene name with h gene name
Rscript bin/ortholog_name_replace.R r/PAS.te.PAusage.scale.bed r/PAS.te.PAusage.scale.rn.bed

# Liftover to other species
liftOver_path=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/liftover/
liftOver -minMatch=0.75 h/PAS.te.cnt.scale.bed \
$liftOver_path/hg38/hg38ToRheMac8.over.chain.gz \
h/PAS.te.cnt.scale.ToRheMac8.bed \
h/PAS.te.cnt.scale.ToRheMac8.unmapped
liftOver -minMatch=0.75 h/PAS.te.cnt.scale.ToRheMac8.bed \
$liftOver_path/rheMac8/rheMac8ToHg38.over.chain.gz \
h/PAS.te.cnt.scale.ToRheMac8.Tohg38.bed \
h/PAS.te.cnt.scale.ToRheMac8.Tohg38.unmapped

# Reciprocol mapped PAS (exclude genes with PAS unmapped or not reciprocol)
join.pl -i1 h/PAS.te.cnt.scale.bed -f1 4 -i2 h/PAS.te.cnt.scale.ToRheMac8.Tohg38.bed -f2 4 |awk '$1==$7 && $2==$8 && $6==$12'|cut -f4|\
join.pl -i1 h/PAS.te.cnt.scale.ToRheMac8.bed -f1 4 -o1 |\
awk '{split($4,a,":");print a[1],$0}' |join.pl -i2 <(cat <(grep -v "#" h/PAS.te.cnt.scale.ToRheMac8.unmapped|cut -f4|cut -f1 -d ":") <(join.pl -i1 h/PAS.te.cnt.scale.bed -f1 4 -i2 h/PAS.te.cnt.scale.ToRheMac8.Tohg38.bed -f2 4 | \
awk '$1!=$7 || $2!=$8 || $6!=$12'|cut -f4 |sort|uniq)) -v -o1|cut -f2- >h/PAS.te.cnt.scale.ToRheMac8.forCmp.bed6

# common PAS
ls h/PAS.te.cnt.scale.ToRheMac8.forCmp.bed6|while read file;do
	awk -v OFS="\t" '{if($2>=30){print $1,$2-30,$3+30,$4,$5,$6}else{print $1,"0",$3+30,$4,$5,$6}}' $file|bedtools intersect -a stdin -b r/PAS.te.PAusage.scale.rn.bed -wo|cut -f4,7-12|join.pl -i1 h/PAS.te.PAusage.scale.bed -f1 4 -f2 1|cut -f1-6,8-|awk '{split($4,a,":");split($10,b,":");if(a[1]==b[1]){print $0}}' >h/PAS.te.PAusage.scale.CommonPA.H.Rpos.tsv;
done;

