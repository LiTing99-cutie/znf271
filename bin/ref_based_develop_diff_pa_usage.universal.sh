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

# 1. get processed terminal exon annoataion
    python bin/terminal_exon_annotation.py \
    --t_e annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt \
    --iso_anno $iso_anno \
    --t_e_out annotation/terminal_exon/terminal_exon_annotation.$suffix.txt \
    --work_path $work_path

# 2. format to gtf 
    bash bin/proximal_distal_gtf.sh annotation/terminal_exon/terminal_exon_annotation.$suffix.txt

# 3. stringtie call rpkm
    bash -c "time bash bin/stringtie.universal.sh annotation/terminal_exon/terminal_exon_annotation.$suffix.gtf \
    analysis/stringtie_whole_genome_$suffix \
    frag_filter" &>log/stringtie_whole_genome_$suffix.log 

# 4. wilcox test
    bash -c "time python bin/wilcox_test.py --rpkm analysis/stringtie_whole_genome_$suffix/stringtie.rpkm.txt \
    --output output/final_res_$suffix.txt" &>log/wilcox.$suffix.log 

# 5. add gene name
    python bin/add_gene_name.py \
    --vs annotation/ensembl_gene_id_vs_symbol.txt \
    --file output/final_res_$suffix.txt