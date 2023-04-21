#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 11 Apr 2023 07:54:32 PM CST
################################################

set -eou pipefail

# scripts

## cp scripts
    cp /home/user/data2/lit/project/ZNF271/02-APA-1/bin/terminal_exon_annotation.loose.py scripts/t_e_anno.py
    cp $project_path/analysis/orf_predict/loose_n_l/run.sh scripts/lncRNA.sh
    cp $project_path/analysis/orf_predict/loose_n_l/run/get_fasta.ipynb scripts
    cp $project_path/analysis/orf_predict/loose_n_l/run/cds_p_to_genomic_p.ipynb scripts
    cp $project_path/run/R/lncRNA.cds.R scripts
    cp $script_path/run/R/ScatterPlot_c.R scripts/ScatterPlot.R
    cp $script_path/run/R/pie_multiple_circles_c.R scripts/pie_multiple_circles.R
    cp $script_path/run/R/GO_c.R scripts/GO.R
    cp /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/miRNA/run.R scripts/miRNA.R


## prepare some files
    less $project_path/annotation/gencode.v41.basic.annotation.gpe |\
    cut -f 1,2,3,12 > $project_path/annotation/map/transcript_id_chr_strand_symbol.txt

    # fragmentation score of samples
    nohup bash /home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/scripts/fragmention.sh \
    /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/mouse 50 mouse_brain \
    /home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output_m \
    /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.genePredExt &>log/frag.mouse_brain.time.log &

    nohup bash /home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/scripts/fragmention.sh \
    /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/rhesus 25 rhesus_brain \
    /home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output_r \
    /home/user/data/lit/database/in_house/rheMac10Plus/rheMac10Plus.addgeneName.gpe &>log/frag.rhesus_brain.time.log &

## define variables
    script_path=/home/user/data2/lit/project/ZNF271/02-APA-1/
    project_path=/home/user/data2/lit/project/ZNF271/02-APA-1/

# 1. get human terminal exons
    # extract chr,strand,terminal exon start,terminal exon end,ensembl gene id
    cd $project_path/
    gtfToGenePred -genePredExt -geneNameAsName2 annotation/gencode.v41.basic.annotation.gtf.gz annotation/gencode.v41.basic.annotation.gpe
    less annotation/gencode.v41.basic.annotation.gpe | \
    awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12,$1}
    else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12,$1}}' | \
    sort -k4,4 > annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt

    Rscript scripts/dominant_te.R
    Rscript scripts/all.exon.R
    Rscript scripts/dominant_te_inter_other_trans.R &> log/dominant_te_inter_other_trans.log
    Rscript scripts/dominant_te_fil.R

    Rscript scripts/dominant_te.R
    Rscript scripts/all_exon.sh
    Rscript scripts/dominant_te_inter_other_trans_other_gene.R

# 2. human PAS

## motif (3 species)
# Galaxy: /rdq/user/lit/project/271/PA/analysis/motif

## compare with 3'-seq and dapars2
bash scripts/compare_3_seq_dapars2.sh
scripts/compare_3_seq_dapars2.R

# 3. PAS map to backbone

## Select PAS cutoff
scripts/cutoff_pas.R
scripts/fil.PAS.R
## Map
scripts/pas_map.R

# 4. Proximal PAS disrupt cds
## Predict lncRNA cds region
scripts/lncRNA.sh
scripts/get_fasta.ipynb
scripts/cds_p_to_genomic_p.ipynb
scripts/lncRNA.cds.R

## Annotate proximal PAS
Rscript $project_path/PAUsage_bed/scripts/cds_anno.R
## Proximal PAS supported
bash $project_path/PAUsage_bed/scripts/proximal_pas_supported.sh

# 5. RNA-seq quantification
## GTF
tail -n +2 output/pro_dis_bin.txt >  output/pro_dis_bin.nh.txt
bash $script_path/bin/proximal_distal_gtf.sh output/pro_dis_bin.nh.txt
## stringtie
nohup bash -c "time bash $script_path/bin/stringtie.universal.sh \
output/pro_dis_bin.nh.gtf \
output/stringtie \
frag_filter" &>log/stringtie.log &
## Wilcoxon test 
cp $script_path/bin/wilcox_test.R scripts/wilcox_test.R
nohup Rscript scripts/wilcox_test.R &

# 6. Description analysis and functional set
## Scatter plot
scripts/ScatterPlot.R


## Pie plot
scripts/pie_multiple_circles.R

## GO plot for functional set
scripts/GO.R

## miRNA bingding for functional set  
### generate proximal or distal bed
bash scripts/te_pro_dis_bed.sh
### intersect with miRNA-binding site
scripts/miRNA.sh
# compare with different functional gene set
scripts/miRNA.R

## disrupt gene set three species

# 7. case

## case APA (three species)

## case fetal expression

## case multiple tissue PDUI
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/znf271_apa_mul_tissue