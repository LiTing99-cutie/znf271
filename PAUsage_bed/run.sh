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
    cp $script_path/analysis/supple/miRNA/run.R scripts/miRNA.R
    cp $script_path/bin/wilcox_test.R scripts/wilcox_test.R
    cp $script_path/bin/develop_case.R scripts/develop_case.R
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
    cp output_r/* /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/rhesus && rm -rf output_r 
    cp output_m/* /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/mouse && rm -rf output_m

    nohup bash /home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/scripts/fragmention.sh \
    /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/r10/ 48 r10_brain \
    /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/r10/ \
    /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/rheMac10/Macaca_mulatta.Mmul_10.107.ucsc.gpe &>log/frag.r10_brain.time.log &

## define variables
    script_path=/home/user/data2/lit/project/ZNF271/02-APA-1/
    project_path=/home/user/data2/lit/project/ZNF271/02-APA-1/

# 1. get human terminal exons
bash scripts/Step_1_te_backbone.sh
    # extract chr,strand,terminal exon start,terminal exon end,ensembl gene id
    cd $project_path/
    gtfToGenePred -genePredExt -geneNameAsName2 annotation/gencode.v41.basic.annotation.gtf.gz annotation/gencode.v41.basic.annotation.gpe
    less annotation/gencode.v41.basic.annotation.gpe | \
    awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12,$1}
    else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12,$1}}' | \
    sort -k4,4 > annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt

    # Rscript scripts/dominant_te.R
    # Rscript scripts/all.exon.R
    # Rscript scripts/dominant_te_inter_other_trans.R &> log/dominant_te_inter_other_trans.log
    # Rscript scripts/dominant_te_fil.R
    Rscript scripts/dominant_te.R
    Rscript scripts/all_exon.sh
    Rscript scripts/dominant_te_inter_other_trans_other_gene.R
#
# 1.1. human PAS
bash scripts/Step_1_1_cutoff_filter_compare.sh
    ## motif (3 species)

    # Galaxy: /rdq/user/lit/project/271/PA/analysis/motif

    ## select cutoff 
    Rscript scripts/cutoff_pas.R

    ## filter
    Rscript scripts/fil.PAS.R

    ## compare with 3'-seq and dapars2
    bash scripts/compare_3_seq_dapars2.sh
    Rscript scripts/compare_3_seq_dapars2.R
#
# 2. PAS map to backbone
## Map
bash scripts/Step_2_pas_map_to_backbone.sh
    Rscript scripts/pas_map.R
#
# 3. RNA-seq quantification [take time]
nohup bash scripts/Step_3_RNA_seq.sh &> log/Step_3_RNA_seq.log &
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
#
# 3.1. Proximal PAS disrupt cds [take time]
nohup bash -c "time bash scripts/Step_3_1_disrupt_cds.sh" &> log/Step_3_1_disrupt_cds.log &
    ## Predict lncRNA cds region
    bash scripts/lncRNA.sh
    Rscript scripts/get_fasta.ipynb
    Rscript scripts/cds_p_to_genomic_p.ipynb
    Rscript scripts/lncRNA.cds.R

    ## Annotate proximal PAS
    Rscript $project_path/PAUsage_bed/scripts/cds_anno.R
    ## Proximal PAS supported
    bash $project_path/PAUsage_bed/scripts/proximal_pas_supported.sh
#
# 4. Description analysis and functional set
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

## case APA and coverage
Rscript scripts/develop_case.R "$PWD/output/stringtie/stringtie.rpkm.txt" \
"/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/human.metadata.clean.txt" \
"/home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.txt" \
"ZNF271P" \
"$PWD/output/"

nohup bash scripts/visual_normalized_bam_h_m.sh &>log/visual_normalized_bam_h_m.log &
nohup bash scripts/visual_normalized_bam_m_r.sh &>log/visual_normalized_bam_m_r.log &
## case fetal expression

## case multiple tissue PDUI
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/znf271_apa_mul_tissue

# 8.adjust the directory structure
cp -r mouse analysis_mouse && rm -rf mouse/
cp -r macaque analysis_rhesus && rm -rf macaque/