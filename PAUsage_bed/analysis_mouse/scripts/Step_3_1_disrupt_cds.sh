#!/usr/bin/sh

################################################
#File Name: scripts/Step_3_1_disrupt_cds.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 11:30:35 PM CST
################################################

set -eo pipefail

gtf=$1
gpe=$2
archive=$3
species=$4
compare=$5
# less $gtf|awk -v OFS='\t' '$3=="gene"{print $10,$12,$14}'|sed 's/"//g;s/;//g' > output/ensembl_gene_id_type_symbol.txt
less $gpe|awk -v OFS='\t' '{print $1,$2,$3,$12}' > output/transcript_id_chr_strand_symbol.txt
## Predict lncRNA cds region
bash scripts/lncRNA.sh $archive $species
## Annotate proximal PAS
Rscript ../scripts/cds_anno.R "$gpe" "output/do_te_per_fil.bed" "output/pro_dis_bin.txt" "output/do_te_pas.txt" "output/do_te_fil.bed" "output/lncRNA_cds_s_e.txt" "output"
## Proximal PAS supported
if [ $compare =="yes" ];then
bash scripts/proximal_pas_supported.sh
fi