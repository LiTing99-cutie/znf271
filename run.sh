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
    annotation/terminal_exon/PAusage.bed6+ \
    /home/user/data2/lit/project/ZNF271/02-APA-1" &> log/ref_based_all_1.log &

    # 3'-UTR numbers
    wc -l /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1.txt
    # 3607

# 3. Plot
    # 3.1 pheatmap output/final_list/final_res_ref_based_all_1.with_geneName.txt | output/R/pheatmap.pdf
        run/R/pheatmap.R
    # 3.2 hist output/final_list/final_res_ref_based_all_1.with_geneName.txt | output/R/hist.pdf;output/R/final_gene_type.pdf;output/final_list/diff.lst
        run/R/tab_and_hist.R
    # 3.3 pie plot 
        # input output/final_list/diff.lst;annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1.txt;
        # annotation/cds_l_transcript_l_ensembl_107.txt;annotation/gencode.v41.annotation.cds.s_e.txt;annotation/cds_lncRNA.rds
        run/R/cds.R
    # 3.4 case developmental stages 
        less /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.txt | \
        grep ENSG00000257267 > \
        /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.ENSG00000257267.txt

        less /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.txt | \
        grep ENSG00000240694 > \
        /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.ENSG00000240694.txt

        run/R/develop_case.R
    

# 4. frag score
    nohup ls /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/*uniq.sorted.bam | parallel --joblog log/index.brain.log -j 20 samtools index {} &

    nohup bash -c "time bash fragmentation_parallel.sh <(ls /home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/*uniq.sorted.bam) 57 \
    log/frag.brain.log /home/user/data2/lit/project/ZNF271/02-APA/analysis/fragmentation/brain" &>log/frag.brain.time.log &

    for frag in $(ls /home/user/data2/lit/project/ZNF271/02-APA/analysis/fragmentation/brain/*frag.txt);do
    perl ../src/calRnaSeqRbScore.modi.pl --fragmentation $frag
    done | cat > /home/user/data2/lit/project/ZNF271/02-APA/analysis/fragmentation/brain/fragmentation.score.txt

# 5. lncRNA orf prediction
    run/sh/cds.annotation.lncRNA.sh
    run/py/get_fasta.ipynb
    run/py/cds_p_to_genomic_p.ipynb
    run/R/lncRNA.cds.R
    # input: annotation/map/gene_id_vs_transcript_id.txt;annotation/lncRNA_orf_prot_co_filter_genomic.txt

# 6. case coverage
    # occupy memory
    # 16min
    nohup bash -c "time bash run/sh/target_bam.sh znf271 /home/user/data2/lit/project/ZNF271/02-APA-1/data/lst/frag_fil_develop_bam.lst" & 
    run/sh/merge_sample.sh

# 7. iso-seq PA identification
    @ Galaxy

# 8. ribo-seq coverage
    ls /home/user/data/lit/project/ZNF271/data/ribo_seq/E-MTAB-7247/human_brain_ribo_*.bam > /home/user/data2/lit/project/ZNF271/02-APA-1/data/lst/ribo.list
    nohup bash -c \
    "time bash /home/user/data2/lit/project/ZNF271/02-APA-1/bin/uniq.tophat2.v2.1.1.sh \
    /home/user/data2/lit/project/ZNF271/02-APA-1/data/lst/ribo.list\
    40\
    /home/user/data/lit/project/ZNF271/data/ribo_seq/E-MTAB-7247" &>log/E-MTAB-7247.ribo.uniq.log &