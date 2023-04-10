#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 24 Mar 2023 09:45:54 PM CST
################################################

set -eou pipefail

spe=human
# bam_lst
[ -d bam_lst ] || mkdir bam_lst
ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/$spe/*sorted.bam | egrep -v -i "\.(brain|forebrain|hindbrain|cerebellum)\.">bam_lst/$spe.lst

# extract uniquely mapped reads ~10h
[ -d log ] || mkdir log
[ -d /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe ] || mkdir /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe
nohup bash -c "time bash /home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.sh bam_lst/$spe.lst 60 /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/$spe" &> log/uniq_bam.wholeBam_$spe.log &

# stringtie ~5h
mkdir rpkm
output_dir=rpkm
annotation=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.gtf
ls /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/*bam|egrep -v -i "\.(brain|forebrain|hindbrain|cerebellum)\." >bam_lst/$spe.uniq_bam.lst
nohup bash stringtie.sh bam_lst/$spe.uniq_bam.lst $annotation $output_dir &> log/stringtie.log &
# merge
bash bin/merge_stringtie_tab.sh $output_dir
# merge with brain
cat rpkm/stringtie.rpkm.txt \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1_loose/stringtie.rpkm.txt \
> rpkm/human.all_tissue.stringtie.rpkm.txt

# output
Rscript run.apa_devo_mul_tissue.R