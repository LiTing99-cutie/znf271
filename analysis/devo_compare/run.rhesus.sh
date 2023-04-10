#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 14 Feb 2023 09:59:25 PM CST
################################################

set -eou pipefail

# mkdir 
mkdir -p rhesus/annotation/terminal_exon/
mkdir -p rhesus/iso_seq
# variable
gtf=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gtf
universal_script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/ref_based_develop_diff_pa_usage.universal.sh
terminal_exon_script=/home/user/data2/lit/project/ZNF271/02-APA-1/bin/terminal_exon_annotation.loose.py

# path
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/PAusage.r.bed6+ rhesus/iso_seq/PAusage.bed6+


[ -f rhesus/annotation/annotation.gpe ] && rm -rf rhesus/annotation/annotation.gpe && \
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gpe rhesus/annotation/annotation.gpe
less rhesus/annotation/annotation.gpe | \
awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12}
else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12}}' | \
sort -k4,4 > rhesus/annotation/terminal_exon.genePredExt

suffix=ref_based_all_1_loose
iso_anno=rhesus/iso_seq/PAusage.bed6+
work_path=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/devo_compare
# add in 2022/12/20
terminal_exon_annotation_script=$terminal_exon_script
# add in 2023/02/14
genePredExt=rhesus/annotation/terminal_exon.genePredExt
# add in 2023/02/15
script_path=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/devo_compare/script
# add in 2023/02/16
spe=rhesus

# 1. get processed terminal exon annoataion
    python $terminal_exon_annotation_script \
    --t_e $genePredExt \
    --iso_anno $iso_anno \
    --t_e_out $spe/annotation/terminal_exon/terminal_exon_annotation.$suffix.txt \
    --work_path $work_path

# 2. format to gtf 
    bash $script_path/proximal_distal_gtf.sh $work_path/$spe/annotation/terminal_exon/terminal_exon_annotation.$suffix.txt

# for rhesus, we need a liftover here
    # /rd1/brick/lixs/Data/General/liftOver/createLiftOver/rheMac10Plus/vsRheMac8/rheMac8TorheMac10Plus.2020-08-26/run.chain/rheMac8.rheMac10Plus.all.chain.gz
    mkdir liftover
    grep ZNF271 rhesus/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.gtf | \
    cut -f 1,4-5,7 > liftover/rheMac8.ZNF271.2Region.bed
    grep ZNF271 rhesus/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.gtf > liftover/rheMac8.ZNF271.2Region.gtf
    cd liftover
    liftOver rheMac8.ZNF271.2Region.bed chain/rheMac8.rheMac10Plus.all.chain.gz \
    rheMac8.rheMac10Plus.mapped rheMac8.rheMac10Plus.unmapped

    # modify in manual (modify strand) 2023/3/3
    touch rheMac8.ZNF271.2Region.rheMac10Plus.gtf
# 3. stringtie call rpkm
    nohup bash -c "time bash $script_path/stringtie.universal.sh liftover/rheMac8.ZNF271.2Region.rheMac10Plus.gtf \
    $work_path/$spe/analysis/stringtie_whole_genome_$suffix \
    frag_filter \
    $spe" &>log/stringtie_whole_genome_$suffix.$spe.log &

# 4. wilcox test
    # bash -c "time python bin/wilcox_test.py --rpkm analysis/stringtie_whole_genome_$suffix/stringtie.rpkm.txt \
    # --output output/final_list/final_res_$suffix.txt" &>log/wilcox.$suffix.log 
    echo "Wilcoxon test [Figure 1B]"
    Rscript $script_path/wilcox_test.R $work_path/$spe/analysis/stringtie_whole_genome_$suffix/stringtie.rpkm.txt