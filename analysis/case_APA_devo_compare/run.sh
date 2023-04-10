#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 08 Feb 2023 08:48:43 PM CST
################################################

set -eou pipefail

export gene_name_m=Zfp35
export gene_name_r=ENSMMUG00000049532

# target region
mkdir target_bed
zless /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.gtf.gz | grep -i $gene_name_m -w | awk '$3=="gene"' | cut -f 1-8  \
> target_bed/target.$gene_name_m.bed
# should be rhemac10plus not rhemac8 for rna-seq data
	less /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gtf | grep -i $gene_name_r -w | awk '$3=="transcript"' | cut -f 1-8 \
	> target_bed/target.$gene_name_r.tmp.bed
	paste -d '\t' <(head -n1 target_bed/target.ZNF271.tmp.bed|cut -f 1-4) <(tail -n1 target_bed/target.ZNF271.tmp.bed|cut -f 5-) > target_bed/target.ZNF271.bed
#
less /home/user/data/lit/database/in_house/rheMac10Plus/rheMac10Plus.addgeneName.gtf | grep -i $gene_name_r -w | awk '$3=="gene"' | cut -f 1-8 \
> target_bed/target.$gene_name_r.bed

cut -f 1,4,5 target_bed/target.$gene_name_m.bed > target_bed/target.$gene_name_m.3c.bed
cut -f 1,4,5 target_bed/target.$gene_name_r.bed > target_bed/target.$gene_name_r.3c.bed
echo -e "9\t76643533\t76655032" > target_bed/target.ENSOCUG00000029705.3c.bed

# select for tools
	# 3.4G
	test_bam=/home/user/data2/uplee/projects/1KMG/rna-seq/2019-nature-moreira/mapping/mouse/E-MTAB-6798.1740sTSm.Mouse.Brain.0dpb.Female.sorted.bam
	time samtools view -hb -L <(cut -f 1,4,5 target_bed/target.Zfp35.bed) $test_bam >test.bam && echo "samtools -L done" &
	time bedtools intersect -a $test_bam -b target_bed/target.Zfp35.bed >test.1.bam && echo "bedtools done" &
	# must provide index file
		samtools view $test_bam chr18:23989632-24005376 
		sambamba slice -L <(cut -f 1,4,5 target_bed/target.Zfp35.bed) $test_bam 
		sambamba view $test_bam chr18:23989632-24005376 
	#
#
# bam_lst
mkdir bam_lst
ls /home/user/data2/uplee/projects/1KMG/rna-seq/2019-nature-moreira/mapping/mouse/*sorted.bam > bam_lst/m.lst
echo "sample numbers of mouse tissue"
less bam_lst/m.lst | cut -d '.' -f 4 | sort | uniq -c
egrep 'Brain|Cerebellum' bam_lst/m.lst > bam_lst/m.brain.lst

ls /home/user/data2/uplee/projects/1KMG/rna-seq/2019-nature-moreira/mapping/rhesus/*sorted.bam > bam_lst/r.lst
echo "sample numbers of rhesus tissue"
less bam_lst/r.lst | cut -d '.' -f 4 | sort | uniq -c
egrep 'Brain|Cerebellum' bam_lst/r.lst > bam_lst/r.brain.lst

ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/rabbit/bam/data/*sorted.bam > bam_lst/rabbit.lst
echo "sample numbers of rabbit tissue"
less bam_lst/rabbit.lst | cut -d '.' -f 3 | sort | uniq -c
egrep 'Brain|Cerebellum' bam_lst/rabbit.lst > bam_lst/rabbit.brain.lst

# extract reads of target region
mkdir -p bam/m 
mkdir -p bam/r
mkdir log
nohup cat bam_lst/m.brain.lst | parallel \
-j 40 --plus --joblog log/target_bam.m.$gene_name_m.log \
'sample={/.};samtools view -hb -L target_bed/target.$gene_name_m.3c.bed {} > bam/m/$sample.target.$gene_name_m.bam' &

nohup cat bam_lst/rhesus.brain.lst | parallel \
-j 40 --plus --joblog log/target_bam.r.$gene_name_r.log \
'sample={/.};samtools view -hb -L target_bed/target.$gene_name_r.3c.bed {} > bam/rhesus/$sample.target.$gene_name_r.bam' &

mkdir -p bam/rabbit
nohup cat bam_lst/rabbit.brain.lst | parallel \
-j 40 --plus --joblog log/target_bam.rabbit.ENSOCUG00000029705.log \
'sample={/.};samtools view -hb -L target_bed/target.ENSOCUG00000029705.3c.bed {} > bam/rabbit/$sample.target.ENSOCUG00000029705.bam' &

# extract unique reads for only target region
for spe in m r;do
ls /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/case_APA_devo_compare/bam/$spe/*.target.*.bam > bam_lst/$spe.brain.target.lst
nohup bash /home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.sh bam_lst/$spe.brain.target.lst 60 bam/$spe \
&> log/uniq_bam.$spe.log &
done

for spe in rhesus;do
ls /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/case_APA_devo_compare/bam/$spe/*.target.ENSMMUG00000049532.bam > bam_lst/$spe.brain.target.lst
nohup bash /home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.sh bam_lst/$spe.brain.target.lst 60 bam/$spe \
&> log/uniq_bam.$spe.log &
done

rm -rf /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/case_APA_devo_compare/bam/$spe/*ZNF271*bam
rm -rf /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/case_APA_devo_compare/bam/$spe/*..bam

for spe in rabbit;do
uniq_script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.NHi1.sh
ls /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/case_APA_devo_compare/bam/$spe/*.target.*.bam > bam_lst/$spe.brain.target.lst
nohup bash $uniq_script bam_lst/$spe.brain.target.lst 60 bam/$spe \
&> log/uniq_bam.$spe.tag.log &
done

# extract unique reads for every species
mkdir /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq
mv bam_lst/m.brain.lst bam_lst/mouse.brain.lst
mv bam_lst/r.brain.lst bam_lst/rhesus.brain.lst
for spe in mouse rhesus rabbit;do
[ -d /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe ] || mkdir /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe
nohup bash /home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.sh bam_lst/$spe.brain.lst 60 /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe &> log/uniq_bam.wholeBam$spe.log &
done

# rabbit reads is mapped to GSNAP, so the cutoff of mapq 60 is not appropriate
# we use grep NH:i:1 to extract, which means number of hits is 1
# but this step is pretty slow, awaiting an improvement later
for spe in rabbit;do
uniq_script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.NHi1.sh
[ -d /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe ] || mkdir /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe
nohup bash $uniq_script bam_lst/$spe.brain.lst 60 /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe &> log/uniq_bam.wholeBam$spe.tag.log &
done

# library size
# bam_stat
	bam_stat=/home/user/BGM/lit/anaconda3/envs/py2/bin/bam_stat.py
	mkdir -p bam_stat/log
	for spe in mouse rhesus rabbit;do
	[ -d bam_stat/$spe ] || mkdir -p bam_stat/$spe
	pushd bam_stat/$spe
	nohup cat ../../bam_lst/$spe.brain.lst | parallel --plus -j 20 --joblog ../log/bam_stat.$spe.log \
	'sample={/.};$bam_stat -q 60 -i {} 1> $sample.txt 2>/dev/null' &
	popd
	done
#
# samtools
for spe in mouse rhesus rabbit;do
pushd bam_stat/$spe
export spe=$spe
[ -d libsize ] || mkdir libsize
nohup ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe/*bam | parallel --plus -j 20 --joblog ../log/bam_stat.$spe.st.log \
'sample={/...};number=$(samtools view -c {});echo -e $sample"\t"$number >libsize/$sample.txt' &
popd
done

## remove previous wrong uniquely mapped bams
rm -rf /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/rabbit/*uniq.bam
## count library size
for spe in rabbit;do
pushd bam_stat/$spe
export spe=$spe
[ -d libsize ] || mkdir libsize
nohup ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe/*bam | parallel --plus -j 20 --joblog ../log/bam_stat.$spe.st.log \
'sample={/...};number=$(samtools view -c {});echo -e $sample"\t"$number >libsize/$sample.txt' &
popd
done

ls bam_stat/rabbit/libsize/ | grep -v sorted | xargs rm -rf
for spe in mouse rhesus rabbit;do
	cat bam_stat/$spe/libsize/* > bam_stat/$spe/libsize.txt
done
sed -i 's/.sorted//g' bam_stat/rabbit/libsize.txt
#
# merge bams of embryos or postnatal stages
mkdir merge_bs
mv m mouse
mv r rhesus
mkdir chromSize
ln -s /home/user/data/lit/database/public/genome/rheMac10Plus/rheMac10Plus.comChr.chrom.sizes chromSize/rhesus.chrom.sizes
ln -s /home/user/data/lit/database/public/genome/mouse/mm10.chrom.sizes chromSize/mouse.chrom.sizes 
ln -s /home/user/data/lit/database/public/genome/rabbit/oryCun2.chrom.sizes chromSize/rabbit.chrom.sizes
mv chromSize/rabbit.chrom.sizes chromSize/rabbit.chrom.sizes.before
less chromSize/rabbit.chrom.sizes.before  | grep -v chrUn | sed 's/chr//g' > chromSize/rabbit.chrom.sizes
bamCoverage=/home/user/BGM/lit/anaconda3/envs/devtools/bin/bamCoverage
for spe in mouse rhesus rabbit;do
for stage in embryo postnatal;do
SN=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/$spe.$stage.sampleName.txt
chromSize=chromSize/$spe.chrom.sizes
libsize=$(less bam_stat/$spe/libsize.txt | egrep -f $SN | awk 'BEGIN{sum=0} {sum+=$2}END {print sum/10^6}')
ls bam/$spe/*.target.*.uniq.* | egrep -f $SN >bam_lst/$spe.$stage.lst
samtools merge -f -o merge_bs/$spe.$stage.bam -b bam_lst/$spe.$stage.lst 
samtools index merge_bs/$spe.$stage.bam && $bamCoverage -b merge_bs/$spe.$stage.bam -o merge_bs/$spe.$stage.bedgragh -of bedgraph -bs 100 -p 10 2>/dev/null
sort -k1,1 -k2,2n merge_bs/$spe.$stage.bedgragh|egrep -v 'alt|random|chrUn|chrM|AAGW|GL|MT'|awk -v ls=$libsize '{print $1,$2,$3,$4/ls/100*10^3}' > merge_bs/$spe.$stage.nm.bedgragh &&
bedGraphToBigWig merge_bs/$spe.$stage.nm.bedgragh $chromSize merge_bs/$spe.$stage.nm.bw
done
done
