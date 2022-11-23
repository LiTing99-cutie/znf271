#!/usr/bin/sh

################################################
#File Name: ../merge_sample.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 30 Sep 2022 01:54:33 PM CST
################################################

set -eou pipefail

pushd /home/user/data2/lit/project/ZNF271/02-APA-1/analysis

# generate pattern file for grep 

for stage in $(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt|sed 's/ /_/g'|cut -f7|tail -n +2|sort|uniq)
do
less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt | \
sed 's/ /_/g' | awk -v stage=$stage '$7==stage' | awk '{print $2}' |sort| uniq > pattern/$stage.txt
ls target_bam/*znf271.bam|egrep -f pattern/$stage.txt|egrep -i "brain|cerebellum" > bam_lst/$stage.bam.lst 
done 

for stage in $(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt|sed 's/ /_/g'|cut -f9|tail -n +2|sort|uniq)
do
less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt | \
sed 's/ /_/g' | awk -v stage=$stage '$9==stage' | awk '{print $2}' |sort| uniq > pattern/$stage.txt
ls target_bam/*znf271.bam|egrep -f pattern/$stage.txt|egrep -i "brain|cerebellum" > bam_lst/$stage.bam.lst 
done 


for bam_lst in $(ls bam_lst/*bam.lst);do
stage=$(basename $bam_lst | sed 's/.bam.lst//g')
samtools merge -f -o bam_lst/$stage.forebrain.bam -b $bam_lst 
done

for bam_lst in $(ls bam_lst/*bam.lst);do
stage=$(basename $bam_lst | sed 's/.bam.lst//g')
genomeCoverageBed -bg -trackline -ibam bam_lst/$stage.forebrain.bam -split > bedgraph/$stage.forebrain.bedgragh &
done
wait

for bam_lst in $(ls bam_lst/*bam.lst);do
stage=$(basename $bam_lst | sed 's/.bam.lst//g')
m_r=$(less /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapped_m_read_period.txt| sed 's/ /_/g'|\
grep $stage|cut -f 2)
less bedgraph/$stage.forebrain.bedgragh|tail -n +2|\
awk -v m_r=$m_r '{print $1,$2,$3,$4/m_r}' > bedgraph/$stage.forebrain.nm.bedgragh
done

for bam_lst in $(ls bam_lst/*bam.lst);do
stage=$(basename $bam_lst | sed 's/.bam.lst//g')
bedGraphToBigWig bedgraph/$stage.forebrain.nm.bedgragh /home/user/data/lit/database/public/genome/hg38/hg38_comChr.chrom.sizes \
bedgraph/$stage.forebrain.nm.bw
done