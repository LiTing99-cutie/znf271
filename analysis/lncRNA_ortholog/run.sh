#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 09 Mar 2023 09:29:40 PM CST
################################################

set -eou pipefail

## 0. get orthologous region based on gentree ##
	mkdir synteny
	cp /tmp/axt_synteny.txt ./synteny

	mkdir gtf
	pushd gtf
	wget https://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz &
	wget https://ftp.ensembl.org/pub/release-73/gtf/macaca_mulatta/Macaca_mulatta.MMUL_1.73.gtf.gz &
	wget https://ftp.ensembl.org/pub/release-73/gtf/mus_musculus/Mus_musculus.GRCm38.73.gtf.gz &
	popd
	mv gtf gtf_v73
#

## Download rbest chain and test ##
	ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/liftover/ ./liftover
	cd liftover
	mkdir rbest
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsRheMac8/reciprocalBest/hg38.rheMac8.rbest.chain.gz -P rbest &
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm10/reciprocalBest/hg38.mm10.rbest.chain.gz -P rbest &

	# test
	awk -v OFS='\t' '{print $1,$4,$5}' /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/target.znf271.bed > target.bed

	# chr18   23992530        24039313
	liftOver -minMatch=0 target.bed liftover/rbest/hg38.rheMac8.rbest.chain.gz hg38.rheMac8.Mapped hg38.rheMac8.unMapped
	liftOver -minMatch=0 target.bed liftover/rbest/hg38.mm10.rbest.chain.gz hg38.mm10.Mapped hg38.mm10.unMapped
##
## 1. liftover based on gene ##

	# gtf
	ln -s ../expr_evo_devo/anno/gtf gtf

	for spe in human mouse;do
	zless gtf/$spe.gtf |awk -v OFS='\t' '$3=="gene" {split($14,z,";");gene_symbol=z[1]; \
	gsub("\"","",gene_symbol);print $1,$4-1,$5,gene_symbol,".",$7}' > gtf/$spe.bed
	done

	for spe in rhesus;do
	zless gtf/$spe.gtf |awk -v OFS='\t' '$3=="transcript" {split($14,z,";");gene_symbol=z[1]; \
	gsub("\"","",gene_symbol);print $1,$4-1,$5,gene_symbol,".",$7}' > gtf/$spe.bed
	done

	grep -f gene_lst.txt -w gtf/human.bed > target.bed
	liftOver -minMatch=0.1 target.bed liftover/rbest/hg38.mm10.rbest.chain.gz hg38.mm10.Mapped hg38.mm10.unMapped
	bedtools intersect -a hg38.mm10.Mapped -b gtf/mouse.bed -wo | awk -v OFS='\t' '{$14=$13/($3-$2);$15=$13/($9-$8);print $0}' > inter.txt

	# intersect_region/orthologous_region maximum 
	# intersect_region/gene_length_in_orthologous_region >0.9
	# no anti-sense gene
	awk '{if($14>x[$4]){x[$4]=$14;y[$4]=$0}} END {for(i in x) {print y[i]}}' inter.txt|awk '$14>0.9&&$15>0.9' |grep -v '\-AS'| wc -l
#


# 2. liftover based on terminal exon ##
	grep -f gene_lst.txt -w ../evo_compare/h/terminal_exon.bed > target.te.bed

	liftOver -minMatch=0.1 target.te.bed liftover/rbest/hg38.mm10.rbest.chain.gz hg38.mm10.te.Mapped hg38.mm10.te.unMapped
	bedtools intersect -a hg38.mm10.te.Mapped -b ../evo_compare/m/terminal_exon.bed -wo| awk -v OFS='\t' '{$14=$13/($3-$2);$15=$13/($9-$8);print $0}' > inter.te.hm.txt
	awk '{if($14>x[$4]){x[$4]=$14;y[$4]=$0}} END {for(i in x) {print y[i]}}' inter.te.hm.txt |cut -f 4,10 |sort -k2,2 -> hm.lnc.txt
	join <(sort -k2,2 ../pro_dis_PA_compare/ortholog/human_mouse_ortholog.txt) hm.lnc.txt -1 2 -2 2 |cut -f 4 -d ' '> h.lnc.rm.txt
	grep -f h.lnc.rm.txt -w -v hm.lnc.txt |sort -k1,1 > hm.lnc.cl.txt

	liftOver -minMatch=0.1 target.te.bed liftover/rbest/hg38.rheMac8.rbest.chain.gz hg38.rheMac8.te.Mapped hg38.rheMac8.te.unMapped
	bedtools intersect -a hg38.rheMac8.te.Mapped -b ../evo_compare/r/terminal_exon.bed -wo| awk -v OFS='\t' '{$14=$13/($3-$2);$15=$13/($9-$8);print $0}' > inter.te.hr.txt
	awk '{if($14>x[$4]){x[$4]=$14;y[$4]=$0}} END {for(i in x) {print y[i]}}' inter.te.hr.txt |cut -f 4,10 |sort -k2,2 -> hr.lnc.txt
	join <(sort -k2,2 ../pro_dis_PA_compare/ortholog/human_macaque_ortholog.txt) hr.lnc.txt -1 2 -2 2 |cut -f 4 -d ' '> h.lnc.rm.txt
	grep -f h.lnc.rm.txt -w -v hr.lnc.txt |sort -k1,1 > hr.lnc.cl.txt

	join hr.lnc.cl.txt hm.lnc.cl.txt > h_lnc_rm.txt

	cat ../pro_dis_PA_compare/ortholog/hrm.ortholog.txt h_lnc_rm.txt > hrm.add_lnc.txt
#