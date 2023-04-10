#!/usr/bin/sh

################################################
#File Name: run.rat.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 16 Feb 2023 09:05:56 PM CST
################################################

set -eou pipefail
export gene_name=ENSRNOG00000049137
export spe=rat
export uniq_script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.NHi1.sh
[ -d /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe ] || mkdir /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe
export uniq_bam_path=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe
export bamCoverage=/home/user/BGM/lit/anaconda3/envs/devtools/bin/bamCoverage

echo -e "18\t15648791\t15664430" > target_bed/target.ENSRNOG00000049137.3c.bed

ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/$spe/bam/data/*sorted.bam > bam_lst/$spe.lst
echo "sample numbers of $spe tissue"
less bam_lst/$spe.lst | cut -d '.' -f 3 | sort | uniq -c
egrep 'Brain|Cerebellum' bam_lst/$spe.lst > bam_lst/$spe.brain.lst

# extract; uniq; sort
mkdir -p bam/$spe
nohup cat bam_lst/$spe.brain.lst | parallel \
-j 40 --plus --joblog log/target_bam.$spe.$gene_name.log \
'sample={/.};samtools view -h -L target_bed/target.$gene_name.3c.bed {} > bam/$spe/$sample.target.$gene_name.bam &&\
(samtools view -H bam/$spe/$sample.target.$gene_name.bam;samtools view bam/$spe/$sample.target.$gene_name.bam|grep -w 'NH:i:1') |\
 samtools sort -@ 60 -o bam/$spe/$sample.target.$gene_name.uniq.sorted.bam -' &

# uniq (take long time; maybe run ahead)
nohup cat bam_lst/$spe.brain.lst | parallel \
-j 30 --plus --joblog log/uniq_bam.wholeBam$spe.tag.log \
'sample={/.};(samtools view -H {};samtools view {}|grep -w 'NH:i:1') | samtools sort -@ 5 -o $uniq_bam_path/$sample.uniq.sorted.bam -' &

# libsize
mkdir -p bam_stat/$spe && pushd bam_stat/$spe
[ -d libsize ] || mkdir libsize
nohup ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe/*bam | parallel --plus -j 40 --joblog ../log/bam_stat.$spe.st.log \
'sample={/...};number=$(samtools view -c {});echo -e $sample"\t"$number >libsize/$sample.txt' &
popd
cat bam_stat/$spe/libsize/* > bam_stat/$spe/libsize.txt
sed -i 's/.sorted//g' bam_stat/$spe/libsize.txt

# bigwig
ln -s /home/user/data/lit/database/public/genome/rat/rn5.chrom.sizes chromSize/$spe.chrom.sizes
mv chromSize/$spe.chrom.sizes chromSize/$spe.chrom.sizes.before
less chromSize/$spe.chrom.sizes.before  | egrep -v 'chrUn|random' | sed 's/chr//g' > chromSize/$spe.chrom.sizes
for spe in $spe;do
for stage in embryo postnatal;do
SN=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/$spe.$stage.sampleName.txt
chromSize=chromSize/$spe.chrom.sizes
libsize=$(less bam_stat/$spe/libsize.txt | egrep -f $SN | awk 'BEGIN{sum=0} {sum+=$2}END {print sum/10^6}')
ls bam/$spe/*.target.*.uniq.* | egrep -f $SN >bam_lst/$spe.$stage.lst
samtools merge -f -o merge_bs/$spe.$stage.bam -b bam_lst/$spe.$stage.lst 
samtools index merge_bs/$spe.$stage.bam && $bamCoverage -b merge_bs/$spe.$stage.bam -o merge_bs/$spe.$stage.bedgragh -of bedgraph -bs 100 -p 10 2>/dev/null
sort -k1,1 -k2,2n merge_bs/$spe.$stage.bedgragh|egrep -v 'alt|random|chrUn|chrM|AAGW|GL|MT|AABR|JH'|awk -v ls=$libsize '{print $1,$2,$3,$4/ls/100*10^3}' > merge_bs/$spe.$stage.nm.bedgragh &&
bedGraphToBigWig merge_bs/$spe.$stage.nm.bedgragh $chromSize merge_bs/$spe.$stage.nm.bw
done
done