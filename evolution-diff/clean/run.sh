#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 23 Nov 2022 10:53:16 PM CST
################################################

set -eou pipefail

# gene_lst
    # new list
        less /home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt | cut -f 2 -d ','|tail -n +2>cds_type_II_III.gene_name.h.lst
        echo "Gene type of type_II_III"
        less /home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt | cut -f 3 -d ','|tail -n +2|sort | uniq -c
        # 2 lncRNA
        # 16 protein_coding
        # 1 transcribed_unitary_pseudogene
    #
#

# homolog in other species
    # ensembl 107 biomart
    # human_macaque_ortholog.txt
    # macaque_mouse_ortholog.txt
    echo -e "ZNF271P ZNF271 Zfp35" >homolog_identified.txt
    grep -f cds_type_II_III.gene_name.h.lst -w human_macaque_ortholog.txt | awk '$3~/ortholog_one2one/'| sort -k1,1 | cut -f2|tr -s '\n'> cds_type_II_III.gene_name.r.lst
    # human macaque mouse
    join -1 2 -2 1 <(grep -f cds_type_II_III.gene_name.h.lst -w human_macaque_ortholog.txt | sort -k2,2) <(grep -f cds_type_II_III.gene_name.r.lst -w macaque_mouse_ortholog.txt| sort -k1,1) | awk -F' ' '$3~/ortholog_one2one/ && $5~/ortholog_one2one/ {print $2,$1,$4}' > homolog.txt
    cat homolog.txt homolog_identified.txt > homolog.add.txt

    cat homolog.add.txt | cut -f 3 -d ' '> cds_type_II_III.gene_name.m.lst
    cat homolog.add.txt | cut -f 2 -d ' '> cds_type_II_III.gene_name.r.lst

#

# softlink
    # human
    [ -f h/cds_type_II_III.gene_name.lst ] && rm -rf h/cds_type_II_III.gene_name.lst && \
    ln -s ../cds_type_II_III.gene_name.h.lst h/cds_type_II_III.gene_name.lst
    # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/PAusage.h.c_anno.bed6+ h/PAusage.bed6+
    # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gpe h/annotation.gpe
    # mouse 
    # [ -f m/annotation.gpe ] && rm -rf m/annotation.gpe && \
    # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.genePredExt m/annotation.gpe
    [ -f m/cds_type_II_III.gene_name.lst ] && rm -rf m/cds_type_II_III.gene_name.lst && \
    ln -s ../cds_type_II_III.gene_name.m.lst m/cds_type_II_III.gene_name.lst
    # [ -f m/PAusage.bed6+ ] && rm -rf m/PAusage.bed6+ && \
    # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/PAusage.m.c_anno.bed6+ m/PAusage.bed6+

    # rhesus
    # [ -f r/annotation.gpe ] && rm -rf r/annotation.gpe && \
    # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gpe r/annotation.gpe
    [ -f r/cds_type_II_III.gene_name.lst ] && rm -rf r/cds_type_II_III.gene_name.lst && \
    ln -s ../cds_type_II_III.gene_name.r.lst r/cds_type_II_III.gene_name.lst
    # [ -f r/PAusage.bed6+ ] && rm -rf r/PAusage.bed6+ && \
    # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/PAusage.r.bed6+ r/PAusage.bed6+

#

# [stable] terminal exon annotaion
    # for spe in m r h;do
    #     less $spe/annotation.gpe | \
    #     awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12}
    #     else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12}}' | \
    #     sort -k4,4 > $spe/terminal_exon.genePredExt

    #     python /home/user/data2/lit/project/ZNF271/02-APA-1/bin/terminal_exon_annotation.loose_for_m_r_h.py \
    #     --t_e terminal_exon.genePredExt \
    #     --iso_anno PAusage.bed6+ \
    #     --t_e_out terminal_exon_annotation.txt \
    #     --work_path /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/$spe
    # done
#

# fill for PA site


    for spe in m r h;do
    pushd $spe
    less annotation.gpe | cut -f 2,3,4,5,6,7,12,1 | awk '{print $0}' > t_id_chr_str_t_s_e_cds_s_e_gene_name.txt
    grep -f cds_type_II_III.gene_name.lst -w PAusage.bed6+ | awk '$5>=2'> PAusage.type_II_III.bed6+
    less terminal_exon_annotation.txt | \
    # awk -v OFS='\t' '{if($3=="+")print $2,$4,$7,$1,".",$3;else print $2,$6-1,$5,$1,".",$3}' | awk '$2>=0'> terminal_exon.bed
    awk -v OFS='\t' '{if($3=="+")print $2,$4,$7,$1,".",$3;else print $2,$6,$5,$1,".",$3}'> terminal_exon.bed
    bedtools intersect -a terminal_exon.bed -b PAusage.type_II_III.bed6+ -wb -s| cut -f 7-12 > PAusage.type_II_III.fil.bed6+
    popd
    done
#

# cds length 

    ## mouse
        # download from biomart
        # cds_l_transcript_l_ensembl_98_m.txt
    ##

    ## rhesus c_anno
        # download from biomart
        # cds_l_transcript_l_ensembl_97_r.txt
    ##    

    ## rhesus
        pushd r

        # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/rheMac8.fa rheMac8.fa
        # ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gtf annotation.gtf

        grep -f cds_type_II_III.gene_name.lst -w annotation.gtf > cds_type_II_III.gtf

        # gff to gtf 
        gffread cds_type_II_III.gtf -o- > cds_type_II_III.gff3

        # 
        gffread -x - -g  rheMac8.fa cds_type_II_III.gff3 > cds_type_II_III.cds.nucl.fa
        seqkit fx2tab -l cds_type_II_III.cds.nucl.fa | sed 's/gene=//g' | tr ' ' '\t'|awk -v OFS='\t' '{print $2,$1,$4,$4/3}' > cds_type_II_III.cds_l_transcript_l.txt
        popd
    ##
#

# clean dir
    # rm -rf homolog.v2.txt
    # rm -rf multi_region.*
    # rm -rf ensembl_107_human_mouse_macaque_homology.clean.txt
    # rm -rf cds_type_II_III.gene_id.r.*
    # rm -rf cds_type_II_III.gene_name.gene_id.r.lst
    # rm -rf t_b2g_b.R r.vs.sorted.txt test_1.ipynb UCSC.h.coordinate.txt ensembl_107_human_mouse_macaque_homology.clean.add.txt cds_type_II_III.gene_name.r.sorted.lst
#

# R 
    Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/run/R/cds_h.R
    Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/run/R/cds_m.R
    Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/run/R/cds_r.R