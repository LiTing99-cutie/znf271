#!/usr/bin/sh

################################################
#File Name: /home/user/data2/lit/project/ZNF271/02-APA/run/sh/$suffix_develop_diff_pa_usage.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 26 Oct 2022 09:29:44 PM CST
################################################

set -eou pipefail

suffix=$1
iso_anno=$2
work_path=$3
# add in 2022/12/20
terminal_exon_annotation_script=$4
# add in 2023/02/14
genePredExt=$5
# add in 2023/02/15
script_path=$6
# add in 2023/02/16
spe=$7

# 1. get processed terminal exon annoataion
    python $terminal_exon_annotation_script \
    --t_e $genePredExt \
    --iso_anno $iso_anno \
    --t_e_out $spe/annotation/terminal_exon/terminal_exon_annotation.$suffix.txt \
    --work_path $work_path

# 2. format to gtf 
    bash $script_path/proximal_distal_gtf.sh $work_path/$spe/annotation/terminal_exon/terminal_exon_annotation.$suffix.txt

# 3. stringtie call rpkm
    bash -c "time bash $script_path/stringtie.universal.sh $work_path/$spe/annotation/terminal_exon/terminal_exon_annotation.$suffix.gtf \
    $work_path/$spe/analysis/stringtie_whole_genome_$suffix \
    frag_filter \
    $spe" &>log/stringtie_whole_genome_$suffix.$spe.log 

# 4. wilcox test
    # bash -c "time python bin/wilcox_test.py --rpkm analysis/stringtie_whole_genome_$suffix/stringtie.rpkm.txt \
    # --output output/final_list/final_res_$suffix.txt" &>log/wilcox.$suffix.log 
    echo "Wilcoxon test [Figure 1B]"
    Rscript $script_path/wilcox_test.R $work_path/$spe/analysis/stringtie_whole_genome_$suffix/stringtie.rpkm.txt

# # 5. add gene name
#     python bin/add_gene_name.py \
#     --vs annotation/map/ensembl_gene_id_type_symbol.txt \
#     --file output/final_list/final_res_$suffix.txt