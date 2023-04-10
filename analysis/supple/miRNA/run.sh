#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 21 Mar 2023 09:00:20 PM CST
################################################

set -eou pipefail

# predicted miRNA-binding sites downloaded from TargetScan
# 23-3-21
wget https://www.targetscan.org/vert_72/vert_72_data_download/Predicted_Target_Locations.default_predictions.hg19.bed.zip
unzip Predicted_Target_Locations.default_predictions.hg19.bed.zip

# liftover
lo_path=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/lncRNA_ortholog/liftover
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz -P $lo_path/hg19/
liftOver Predicted_Target_Locations.default_predictions.hg19.bed $lo_path/hg19/hg19ToHg38.over.chain.gz Predicted_Target_Locations.default_predictions.hg38.bed unMapped

# get refined terminal exon coordinate
# also get proximal region and distal region
# /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/liftover_inter.sh
te=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.txt
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.pro.bed hg38.pro.bed
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.dis.bed hg38.dis.bed
awk -v OFS='\t' '{if($3=="+"){print $2,$4,$7,$1,".",$3}else{print $2,$6-1,$5,$1,".",$3}}' $te > hg38.te.bed

# APA and miRNA-binding sites correlation
for type in te pro dis;do
bedtools intersect -a hg38.$type.bed -b Predicted_Target_Locations.default_predictions.hg38.bed -wo | awk -v OFS='\t' '$20=$19/($9-$8)' | awk '$20>0.5'> hg38.$type.miRNA.txt
done

# Transcriptome-wide Discovery of microRNA Binding Sites in Human Brain
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52082/suppl/GSE52082_Ago2_binding_clusters.BED.gz
gunzip -c GSE52082_Ago2_binding_clusters.BED.gz > GSE52082_Ago2_binding_clusters.BED
liftOver GSE52082_Ago2_binding_clusters.BED $lo_path/hg19/hg19ToHg38.over.chain.gz GSE52082_Ago2_binding_clusters.hg38.BED GSE52082_Ago2_binding_clusters.unMapped

mkdir AGO2
cd AGO2
for type in te pro dis;do
bedtools intersect -a ../hg38.$type.bed -b ../GSE52082_Ago2_binding_clusters.hg38.BED -wo | awk -v OFS='\t' '$14=$13/($9-$8)' | awk '$14>0.5'> hg38.$type.miRNA.txt
done

# intersect between AGO2 and TargetScan
cd ..
bedtools intersect -a GSE52082_Ago2_binding_clusters.hg38.BED -b Predicted_Target_Locations.default_predictions.hg38.bed