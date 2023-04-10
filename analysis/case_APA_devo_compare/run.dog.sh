#!/usr/bin/sh

################################################
#File Name: run.dog.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 08 Mar 2023 07:09:11 PM CST
################################################

set -eou pipefail

echo -e "7\t54685661\t54698702" > target_bed/target.ENSCAFG00845006933.3c.bed

ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/dog/bam/*bam > bam_lst/dog.lst

mkdir -p bam/dog
nohup cat bam_lst/dog.lst | parallel \
-j 10 --plus --joblog log/target_bam.dog.ENSCAFG00845006933.log \
'sample={/.};samtools view -hb -L target_bed/target.ENSCAFG00845006933.3c.bed {} > bam/dog/$sample.target.ENSCAFG00845006933.bam' &


## count library size
spe=dog
mkdir -p bam_stat/$spe && pushd bam_stat/$spe
export spe=$spe
[ -d libsize ] || mkdir libsize
nohup ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/dog/bam/*bam | parallel --plus -j 10 --joblog ../log/bam_stat.$spe.st.log \
'sample={/...};number=$(samtools view -c {});echo -e $sample"\t"$number >libsize/$sample.txt' &
popd


cat bam_stat/$spe/libsize/* > bam_stat/$spe/libsize.txt

# sample names of embryo or postnatal
cd /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/dog
ls data | grep SRR | grep adult | cut -f 1 -d "_" | sort|uniq > postnatal.sampleName.txt
ls data | grep SRR | grep -v adult | cut -f 1 -d "_" | sort|uniq > embryo.sampleName.txt

# merge bams of embryos or postnatal stages
ln -s /home/user/data/lit/database/public/genome/dog/ROS_Cfam_1.0/Canis_lupus_familiaris.ROS_Cfam_1.0.chrom.sizes chromSize/dog.chrom.sizes
# mv chromSize/rabbit.chrom.sizes chromSize/rabbit.chrom.sizes.before
# less chromSize/rabbit.chrom.sizes.before  | grep -v chrUn | sed 's/chr//g' > chromSize/rabbit.chrom.sizes
bamCoverage=/home/user/BGM/lit/anaconda3/envs/devtools/bin/bamCoverage
for spe in dog;do
for stage in embryo postnatal;do
SN=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/dog/$stage.sampleName.txt
chromSize=chromSize/$spe.chrom.sizes
libsize=$(less bam_stat/$spe/libsize.txt | egrep -f $SN | awk 'BEGIN{sum=0} {sum+=$2}END {print sum/10^6/2}')
ls bam/$spe/*bam | egrep -f $SN >bam_lst/$spe.$stage.lst
samtools merge -f -o merge_bs/$spe.$stage.bam -b bam_lst/$spe.$stage.lst 
samtools index merge_bs/$spe.$stage.bam && $bamCoverage -b merge_bs/$spe.$stage.bam -o merge_bs/$spe.$stage.bedgragh -of bedgraph -bs 100 -p 10 2>/dev/null
sort -k1,1 -k2,2n merge_bs/$spe.$stage.bedgragh|egrep -v 'alt|random|chrUn|chrM|AAGW|GL|MT'|awk -v ls=$libsize '{print $1,$2,$3,$4/ls/100*10^3}' > merge_bs/$spe.$stage.nm.bedgragh &&
bedGraphToBigWig merge_bs/$spe.$stage.nm.bedgragh $chromSize merge_bs/$spe.$stage.nm.bw
done
done

# dapars2 annotate PA site
mkdir dapars && cd dapars
SCRIPT_PATH=/home/user/data2/lit/project/ZNF271/02-APA/src/DaPars2-master/src
gene_model=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/dog/Canis_lupus_familiaris.ROS_Cfam_1.0.wholeGene_annotation.bed
mapping=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/dog/IDmapping.txt
python $SCRIPT_PATH/DaPars_Extract_Anno.py -b $gene_model -s $mapping -o ROS_Cfam_1_3UTR_annotation.bed

echo "7" > chrList.txt
mv /home/user/data2/lit/project/ZNF271/02-APA/annotation/DAPARS2/Dapars2_configure_file ./

ls ../bam/dog/*bam| parallel -j 18 --plus --joblog wig.log 'sample=`basename -s ".uniq.sorted.target.ENSCAFG00845006933.bam" {}`;bedtools genomecov -ibam {} -bga -split -trackline > $sample.wig'
ls *wig | tr '\n' ','

ls ../bam/dog/*bam | parallel --plus -j 18 --joblog cov.log 'sample=`basename -s ".uniq.sorted.target.ENSCAFG00845006933.bam" {}`;number=$(samtools view -c {});echo -e $sample".wig""\t"$number >>$sample.cov.txt'

cat *.cov.txt > mapping_wig_location_with_depth.txt

python $SCRIPT_PATH/DaPars2_Multi_Sample_Multi_Chr.py Dapars2_configure_file chrList.txt

# whole genome
mkdir whole_genome
cd whole_genome

ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/dog/bam/*bam| parallel -j 18 --plus --joblog wig.log 'sample=`basename -s ".uniq.sorted.bam" {}`;bedtools genomecov -ibam {} -bga -split -trackline > $sample.wig' &
less ../../bam_stat/dog/libsize.txt | awk '{print $1".wig""\t"$2}' >mapping_wig_location_with_depth.txt
nohup python $SCRIPT_PATH/DaPars2_Multi_Sample_Multi_Chr.py Dapars2_configure_file_1 ../chrList.txt &>darpars.log &
