#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 23 Nov 2022 10:53:16 PM CST
################################################

set -eou pipefail

# gene_lst
    less /home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/cds_type.txt | egrep -i 'another|non_coding' \
    | cut -f 2 -d ' ' > cds_type_II_III.gene_name.h.lst
#

# softlink
    # human
    ln -s ../cds_type_II_III.gene_name.h.lst h/cds_type_II_III.gene_name.lst
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/PAusage.h.c_anno.bed6+ h/PAusage.bed6+
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gpe h/annotation.gpe
    # mouse 
    [ -f m/annotation.gpe ] && rm -rf m/annotation.gpe && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.genePredExt m/annotation.gpe
    [ -f m/cds_type_II_III.gene_name.lst ] && rm -rf m/cds_type_II_III.gene_name.lst && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/cds_type_II_III.gene_name.m.lst \
    m/cds_type_II_III.gene_name.lst
    [ -f m/PAusage.bed6+ ] && rm -rf m/PAusage.bed6+ && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/PAusage.m.c_anno.bed6+ m/PAusage.bed6+

    # rhesus
    [ -f r/annotation.gpe ] && rm -rf r/annotation.gpe && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gpe r/annotation.gpe
    [ -f r/cds_type_II_III.gene_name.lst ] && rm -rf r/cds_type_II_III.gene_name.lst && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/cds_type_II_III.gene_name.r.lst r/cds_type_II_III.gene_name.lst
    [ -f r/PAusage.bed6+ ] && rm -rf r/PAusage.bed6+ && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/PAusage.r.bed6+ r/PAusage.bed6+

    # r_c_anno
    [ -f r_c_anno/annotation.gpe ] && rm -rf r_c_anno/annotation.gpe && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Macaca_mulatta.Mmul_8.0.1.97.add.gpe r_c_anno/annotation.gpe
    [ -f r_c_anno/cds_type_II_III.gene_name.lst ] && rm -rf r_c_anno/cds_type_II_III.gene_name.lst && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/cds_type_II_III.gene_id.r.add.lst r_c_anno/cds_type_II_III.gene_name.lst
    [ -f r_c_anno/PAusage.bed6+ ] && rm -rf r_c_anno/PAusage.bed6+ && \
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/PAusage.r.c_anno.bed6+ r_c_anno/PAusage.bed6+

#

# homolog in other species
    # -w match whole word only
    echo -e "ZNF271P\tZfp35\tZNF271" >homolog_identified.txt
    # echo -e "ZNF271P\tZfp35\tZNF271" >homolog_identified.txt
    cat ensembl_107_human_mouse_macaque_homology.clean.txt homolog_identified.txt > ensembl_107_human_mouse_macaque_homology.clean.add.txt

    grep -f cds_type_II_III.gene_name.h.lst -w ensembl_107_human_mouse_macaque_homology.clean.add.txt > homolog.txt

    # 17
    grep -f cds_type_II_III.gene_name.h.lst -w ensembl_107_human_mouse_macaque_homology.clean.add.txt | \
    cut -f 2 > cds_type_II_III.gene_name.m.lst

    # 
    grep -f cds_type_II_III.gene_name.h.lst -w ensembl_107_human_mouse_macaque_homology.clean.add.txt | \
    cut -f 3 > cds_type_II_III.gene_name.r.lst

    sort -k1,1 cds_type_II_III.gene_name.r.lst > cds_type_II_III.gene_name.r.sorted.lst
    sort -k2,2 /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/ensembl_gene_id_vs_symbol.txt > r.vs.sorted.txt

    join -1 1 -2 2 cds_type_II_III.gene_name.r.sorted.lst r.vs.sorted.txt > cds_type_II_III.gene_name.gene_id.r.lst
    less cds_type_II_III.gene_name.gene_id.r.lst | cut -f 2 -d ' ' > cds_type_II_III.gene_id.r.lst
    sed '1i ENSMMUG00000049532' cds_type_II_III.gene_id.r.lst > cds_type_II_III.gene_id.r.add.lst
    
#

# terminal exon annotaion
    for spe in m r r_c_anno h;do
        less $spe/annotation.gpe | \
        awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12}
        else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12}}' | \
        sort -k4,4 > $spe/terminal_exon.genePredExt

        python /home/user/data2/lit/project/ZNF271/02-APA-1/bin/terminal_exon_annotation.loose_for_m_r_h.py \
        --t_e terminal_exon.genePredExt \
        --iso_anno PAusage.bed6+ \
        --t_e_out terminal_exon_annotation.txt \
        --work_path /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/$spe
    done
#

# fill for PA site


    for spe in m r r_c_anno h;do
    pushd $spe
    less annotation.gpe | cut -f 2,3,4,5,6,7,12,1 | awk '{print $0}' > t_id_chr_str_t_s_e_cds_s_e_gene_name.txt
    grep -f cds_type_II_III.gene_name.lst -w PAusage.bed6+ | awk '$5>=2'> PAusage.type_II_III.bed6+
    less terminal_exon_annotation.txt | \
    awk -v OFS='\t' '{if($3=="+")print $2,$4,$7,$1;else print $2,$6-1,$5,$1}' | awk '$2>=0'> terminal_exon.bed
    bedtools intersect -a terminal_exon.bed -b PAusage.type_II_III.bed6+ -wb | cut -f 5-10 > PAusage.type_II_III.fil.bed6+
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
        cd r

        ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/rheMac8.fa rheMac8.fa
        ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gtf annotation.gtf

        grep -f cds_type_II_III.gene_name.lst \
        -w annotation.gtf > cds_type_II_III.gtf

        # gff to gtf 
        gffread cds_type_II_III.gtf -o- > cds_type_II_III.gff3

        # 
        gffread -x - -g  rheMac8.fa cds_type_II_III.gff3 > cds_type_II_III.cds.nucl.fa
        seqkit fx2tab -l cds_type_II_III.cds.nucl.fa | sed 's/gene=//g' | tr ' ' '\t'|awk -v OFS='\t' '{print $2,$1,$4,$4/3}' > cds_type_II_III.cds_l_transcript_l.txt
    ##
#