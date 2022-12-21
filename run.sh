#!/usr/bin/sh
set -eou pipefail

# 1. get terminal exons
    # extract chr,strand,terminal exon start,terminal exon end,ensembl gene id
    gtfToGenePred -genePredExt -geneNameAsName2 annotation/gencode.v41.basic.annotation.gtf.gz annotation/gencode.v41.basic.annotation.gpe
    less annotation/gencode.v41.basic.annotation.gpe | \
    awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12}
    else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12}}' | \
    sort -k4,4 > annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt

# 2. Five steps
    # get processed terminal exon annoataion | annotation/terminal_exon/terminal_exon_annotation.$suffix.txt
    # format to gtf | annotation/terminal_exon/terminal_exon_annotation.$suffix.gtf
    # stringtie call rpkm | analysis/stringtie_whole_genome_$suffix/stringtie.rpkm.txt
    # wilcox test | output/final_res_$suffix.txt
    # add gene name | output/final_res_$suffix.with_geneName.txt
    nohup bash -c "time bash bin/ref_based_develop_diff_pa_usage.universal.sh \
    ref_based_all_1 \
    annotation/terminal_exon/PAusage.h.c_anno.bed6+ \
    /home/user/data2/lit/project/ZNF271/02-APA-1" &> log/ref_based_all_1.log &

    # 3'-UTR numbers
    wc -l /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1.txt
    # 3607

# 3. Plot
    # 3.1 [deprecated] pheatmap output/final_list/final_res_ref_based_all_1.with_geneName.txt | output/R/pheatmap.pdf
        run/R/pheatmap.R
    # 3.2 [deprecated] hist output/final_list/final_res_ref_based_all_1.with_geneName.txt | output/R/hist.pdf;output/R/final_gene_type.pdf;output/final_list/diff.lst
        run/R/tab_and_hist.R
    # 3.3 [deprecated] pie plot 
        # input output/final_list/diff.lst;annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1.txt;
        # annotation/cds_l_transcript_l_ensembl_107.txt;annotation/gencode.v41.annotation.cds.s_e.txt;annotation/cds_lncRNA.rds
        less annotation/gencode.v41.basic.annotation.gpe | cut -f 1,3,6,7,12 > annotation/gencode.v41.basic.annotation.cds.s_e.txt
        run/R/cds.R
    # 3.4 case developmental stages 
        # ZDHHC8
            less /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.txt | \
            grep -w ZDHHC8> \
            /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.ZDHHC8.txt
            
            mkdir /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ZDHHC8/
            Rscript bin/develop_case.R ZDHHC8
        #

        # ZNF271P
                less /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.txt | \
                grep -w ZNF271P> \
                /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.ZNF271P.txt
                mkdir /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ZNF271P/
                Rscript bin/develop_case.R ZNF271P
        #
        # CXCL1
            less /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.txt | \
            grep -w CX3CL1> \
            /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.CX3CL1.txt
            mkdir /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/CX3CL1/
            Rscript bin/develop_case.R CX3CL1
        #
    #
    # 3.6 cds effect
        run/R/cds.R
    #
    # 3.5 ScatterPlot
        run/R/ScatterPlot.R
    #
    # 3.7 multiple pies
        run/R/pie_multiple_circles.R
    #
    # 3.6 evo pheatmap
        evolution-diff/clean/run.sh
        run/R/cds_h.R
        run/R/cds_m.R
        run/R/cds_r.R
        run/R/cds_compara.R
    #


# 4. [stable] frag score
    nohup ls /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/*uniq.sorted.bam | parallel --joblog log/index.brain.log -j 20 samtools index {} &

    nohup bash -c "time bash fragmentation_parallel.sh <(ls /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/*uniq.sorted.bam) 57 \
    log/frag.brain.log /home/user/data2/lit/project/ZNF271/02-APA/analysis/fragmentation/brain" &>log/frag.brain.time.log &

    for frag in $(ls /home/user/data2/lit/project/ZNF271/02-APA/analysis/fragmentation/brain/*frag.txt);do
    perl ../src/calRnaSeqRbScore.modi.pl --fragmentation $frag
    done | cat > /home/user/data2/lit/project/ZNF271/02-APA/analysis/fragmentation/brain/fragmentation.score.txt

# 5. [change in need] lncRNA orf prediction

    run/py/get_fasta.ipynb
    run/sh/cds.annotation.lncRNA.sh
    run/py/cds_p_to_genomic_p.ipynb
    run/R/lncRNA.cds.R
    # input: annotation/map/gene_id_vs_transcript_id.txt;annotation/lncRNA_orf_prot_co_filter_genomic.txt

# 6. [stable] case coverage
    # occupy memory
    # 16min
    nohup bash -c "time bash run/sh/target_bam.sh znf271 /home/user/data2/lit/project/ZNF271/02-APA-1/data/lst/frag_fil_develop_bam.lst" & 
    run/sh/merge_sample.sh

# 7. [In Galaxy] iso-seq PA identification
    @ Galaxy

# 8. [stable] ribo-seq coverage
    ls /home/user/data/lit/project/ZNF271/data/ribo_seq/E-MTAB-7247/human_brain_ribo_*.bam > /home/user/data2/lit/project/ZNF271/02-APA-1/data/lst/ribo.list
    nohup bash -c \
    "time bash /home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.tophat2.v2.1.1.sh \
    /home/user/data2/lit/project/ZNF271/02-APA-1/data/lst/ribo.list\
    40\
    /home/user/data/lit/project/ZNF271/data/ribo_seq/E-MTAB-7247" &>log/E-MTAB-7247.ribo.uniq.log &
# 9. [deprecated] species compare in UCSC genome browser
    cd /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/
    lst=/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/cds_type_II_III.gene_name.h.lst
    gtf=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gtf.gz
    homolog=/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/homolog.txt
    zless $gtf | grep -f $lst -w| awk -v OFS='\t' '$3=="gene" {split($14,z,";");gene_symbol=z[1];gsub("\"","",gene_symbol);print $1,$4,$5,gene_symbol}' > multi_region.bed
    echo -e "#database hg38\n#shortDesc typeII_III_gene\n#padding 6" > multi_region.header.txt
    cat multi_region.header.txt multi_region.bed >multi_region.withHeader.bed
    ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/multi_region.withHeader.bed /home/user/data2/rbase/ucsc2/htdocs/data/lit/
    join -1 4 -2 1 -t $'\t' <(sort -k4,4 multi_region.bed) <(sort -k1,1 $homolog|cut -f1) | awk -v OFS='\t' '{print $1,$2":"$3"-"$4}'

    lst=/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/cds_type_II_III.gene_name.m.lst
    gtf=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.gtf.gz
    zless $gtf | grep -f $lst -w| awk -v OFS='\t' '$3=="gene" {split($14,z,";");gene_symbol=z[1];gsub("\"","",gene_symbol);print $1,$4,$5,gene_symbol}' > multi_region.m.bed
    less multi_region.m.bed | awk -v OFS='\t' '{print $4,$1":"$2"-"$3}'

    lst=/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/cds_type_II_III.gene_name.r.lst
    gpe=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/Final.426909.withName.gpe
    zless $gpe | grep -f $lst -w| cut -f 2,4,5,12 > multi_region.r.bed
    # get gene bed from transcript bed
    /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/t_b2g_b.R
    less multi_region.r.gene.bed | awk -v OFS='\t' '{print $1,$2":"$3"-"$5}'
#

# 10. loose for terminal exon annotation
    # 5116 to 7529 terminal exons of genes
    python bin/terminal_exon_annotation.loose.py \
    --t_e annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt \
    --iso_anno annotation/terminal_exon/PAusage.h.c_anno.bed6+ \
    --t_e_out annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.txt \
    --work_path /home/user/data2/lit/project/ZNF271/02-APA-1

    # plan to run in 22/12/20, change suffix
        nohup bash -c "time bash bin/ref_based_develop_diff_pa_usage.universal.sh \
        ref_based_all_1_loose \
        annotation/terminal_exon/PAusage.h.c_anno.bed6+ \
        /home/user/data2/lit/project/ZNF271/02-APA-1 \
        bin/terminal_exon_annotation.loose.py" &> log/ref_based_all_1_loose.log &
#