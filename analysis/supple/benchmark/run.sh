#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 22 Mar 2023 07:27:37 PM CST
################################################

set -eou pipefail


# 1. iso-seq pa annotation
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/h/PAusage.bed6+ ./PAusage.bed6+
## select only common chromosomes; convinced sites and sort
grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]' PAusage.bed6+| awk '$5>1'|sort -k1,1 -k2,2n  > PAusage.comChr.convinced.bed6+

# 2. public pa annotation
wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
gunzip -c atlas.clusters.2.0.GRCh38.96.bed.gz > atlas.clusters.2.0.GRCh38.96.bed
## 1-based to 0-based
less atlas.clusters.2.0.GRCh38.96.bed |cut -f 4| awk -v OFS='\t' -F ":" '{print $1,$2-1,$2,".",".",$3}' > atlas.pos.bed
## updown 30 bp
less atlas.pos.bed|awk -v OFS='\t' '{if($2-30>=0)print $1,$2-30,$3+30,$4,$5,$6;else print $1,1,$3,$4,$5,$6}'> atlas.upDown30.pos.bed
## select only common chromosomes and sort
grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]' atlas.upDown30.pos.bed |awk '{print "chr"$1,$2,$3,$4,$5,$6}'|sort -k1,1 -k2,2n > atlas.upDown30.pos.comChr.bed

# 3. dapars identification of PA sites using RNA-seq data
# export uniq_bam_path=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/human
# mkdir $uniq_bam_path/wig
# mkdir $uniq_bam_path/cov
# ls $uniq_bam_path/*bam| nohup parallel -j 30 --plus --joblog wig.log 'sample=`basename -s ".sorted.uniq.sorted.bam" {}`;bedtools genomecov -ibam {} -bga -split -trackline > $uniq_bam_path/wig/$sample.wig' &
# ls $uniq_bam_path/*bam | nohup parallel --plus -j 30 --joblog cov.log 'sample=`basename -s ".sorted.uniq.sorted.bam" {}`;number=$(samtools view -c {});echo -e $sample".wig""\t"$number >>$uniq_bam_path/cov/$sample.cov.txt' &
# cat *.cov.txt > mapping_wig_location_with_depth.txt
# ls *wig | tr '\n' ','

# 3.1 cp data in data dir to data3 dir
pushd /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/
mkdir human_ns
nohup find /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/ -type f |xargs -P 40 -I {} cp {} human_ns/ &
popd

# 3.2 prepare wig files
mkdir dapars && cd dapars
export uniq_bam_path=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/human_ns
mkdir $uniq_bam_path/wig
# mkdir $uniq_bam_path/cov
ls $uniq_bam_path/*bam| nohup parallel -j 30 --plus --joblog wig.log 'sample=`basename -s ".sorted.uniq.sorted.bam" {}`;bedtools genomecov -ibam {} -bga -split -trackline > $uniq_bam_path/wig/$sample.wig' &
ln -s $uniq_bam_path/wig/*wig ./

# 3.3 coverage files
cp $uniq_bam_path/library_size.txt library_size.txt
sed 's/.sorted/.wig/g' library_size.txt > mapping_wig_location_with_depth.txt

# 3.4 3'UTR annotation files
SCRIPT_PATH=/home/user/data2/lit/project/ZNF271/02-APA/src/DaPars2-master/src
Anno_path=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/human_dapars_anno/
gene_model=$Anno_path/*wholeGene_annotation.bed
mapping=$Anno_path/IDmapping.txt
python $SCRIPT_PATH/DaPars_Extract_Anno.py -b $gene_model -s $mapping -o 3UTR_annotation.bed

# 3.5 config file
cp /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/case_APA_devo_compare/dapars/Dapars2_configure_file ./
cp /home/user/data2/lit/project/ZNF271/02-APA/src/DaPars2-master/Dapars2_Test_Dataset/chrList.txt ./

# 3.6 run
## Aligned_Wig_files
ls *wig | tr '\n' ','
## run 
nohup bash -c "time python $SCRIPT_PATH/DaPars2_Multi_Sample_Multi_Chr.py Dapars2_configure_file chrList.txt" &> run.log &

# 3.7 merge pa sites of all chrs
mkdir Dapars2_test_chr
cat Dapars2_test_chr*/Dapars2_result_temp.chr*.txt > Dapars2_test_chr/Dapars2_result_temp.chr.txt
## remove duplicates
grep -v fit_value Dapars2_test_chr/Dapars2_result_temp.chr.txt|awk -v FS='|' '{print $3,$4,$5}'|awk '{print $1,$4-1,$4,".",".",$2,$1"_"$4"_"$2}' | \
awk '!a[$7]++{for (i=1;i<=6;i++)printf ("%s\t",$i);print ""}' > Dapars2_test_chr/Dapars2_result_temp.chr.bed6

# 4. updown15 & common chr include chr1-22 and chrXY & sort
cd ..
function updown15(){
awk '{if($2-15>=0)print $1,$2-15,$3+15,$4,$5,$6;else print $1,1,$3,$4,$5,$6}'|grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]'|sort -k1,1 -k2,2n
}
less PAusage.comChr.convinced.bed6+ | updown15 > PAusage.comChr.convinced.upDown15.bed6+
less dapars/Dapars2_test_chr/Dapars2_result_temp.chr.bed6 | updown15 > Dapars2_result_temp.chr.upDown15.bed6
## ensembl chromosome name to ucsc chromosome name
less atlas.pos.bed|awk '{print "chr"$1,$2,$3,$4,$5,$6}'|updown15 > atlas.upDown15.pos.comChr.bed

# 5. bedtools intersect
function inter(){
	bedtools intersect -a $1 -b $2 -wa -s| wc -l
}

# 6. write to output
inter_2_1=$(inter PAusage.comChr.convinced.upDown15.bed6+ atlas.upDown15.pos.comChr.bed)
# inter_3_1=$(inter Dapars2_result_temp.chr.upDown15.bed6 atlas.upDown15.pos.comChr.bed)
inter_2_3=$(inter PAusage.comChr.convinced.upDown15.bed6+ Dapars2_result_temp.chr.upDown15.bed6 )
# n_1=$(wc -l atlas.upDown15.pos.comChr.bed|cut -f1 -d ' ')
n_2=$(wc -l PAusage.comChr.convinced.upDown15.bed6+|cut -f1 -d ' ')
# n_3=$(wc -l Dapars2_result_temp.chr.upDown15.bed6|cut -f1 -d ' ')
echo -e "inter_2_1\tinter_2_3\tn_2" > intersect.txt
echo -e "$inter_2_1\t$inter_2_3\t$n_2" >> intersect.txt