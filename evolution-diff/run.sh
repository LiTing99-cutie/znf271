#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 21 Nov 2022 05:07:14 PM CST
################################################

set -eou pipefail

ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/PAusage.bed6+ PAusage.bed6+

grep ENSG00000240694.9 PAusage.bed6+ | awk -v OFS='\t' '{print $1":"$2"-"$3,$4,$5,$6,$7}' > PAusage.ENSG00000240694.9.bed6+
grep ENSG00000171817.17 PAusage.bed6+ | awk -v OFS='\t' '{print $1":"$2"-"$3,$4,$5,$6,$7}' > PAusage.ENSG00000171817.17.bed6+

less /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.txt | \
grep ENSG00000171817 > \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.ENSG00000171817.txt
mkdir /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ENSG00000171817/
Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/run/R/develop_case.R ENSG00000171817

# draw plot for visual
less /home/user/data2/lit/project/ZNF271/02-APA/output/final_list/cds_type.txt | grep -i another | cut -f 1 -d ' ' > cds_type_II.gene.lst

for gene in $(cat cds_type_II.gene.lst);do
gene_no_suffix=$(echo $gene | cut -d '.' -f 1)
less /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.txt | \
grep $gene_no_suffix > \
/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1/stringtie.rpkm.$gene_no_suffix.txt
mkdir /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/$gene_no_suffix/
Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/run/R/develop_case.R $gene_no_suffix
done

less ../annotation/gencode.v41.annotation.gpe | cut -f 2,3,4,5,6,7,12,1 | awk '{print $0}' > t_id_chr_str_t_s_e_cds_s_e_gene_name.txt

grep -f cds_type_II.gene.lst PAusage.bed6+ | awk '$5>=2'> PAusage.type_II.bed6+

less /home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1.txt | \
awk -v OFS='\t' '{if($3=="+")print $2,$4,$7,$1;else print $2,$6-1,$5,$1}' > terminal_exon.bed

bedtools intersect -a terminal_exon.bed -b PAusage.type_II.bed6+ -wb | cut -f 5-10 > PAusage.type_II.fil.bed6+

# cp from galaxy
# ensembl_107_human_mouse_macaque_homology.clean.txt

# -w match whole word only
grep -f cds_type_II.gene.gene_name.lst -w ensembl_107_human_mouse_macaque_homology.clean.txt | cut -f 2 > cds_type_II.gene.gene_name.m.lst
grep -f cds_type_II.gene.gene_name.lst -w ensembl_107_human_mouse_macaque_homology.clean.txt | cut -f 3 > cds_type_II.gene.gene_name.r.lst

less gencode.vM23.annotation.genePredExt | \
awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12}
else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12}}' | \
sort -k4,4 > m.terminal_exon.genePredExt

less Final.426909.withName.gpe | \
awk -v OFS='\t' '{if($3=="+"){split($9,x,",");split($10,y,",");x_l=length(x);y_l=length(y);print $2,$3,x[x_l-1],y[y_l-1],$12}
else {split($9,x,",");split($10,y,",");print $2,$3,x[1],y[1],$12}}' | \
sort -k4,4 > r.terminal_exon.genePredExt

for spe in m;do
python /home/user/data2/lit/project/ZNF271/02-APA-1/run/py/terminal_exon_annotation.py \
--t_e ${spe}.terminal_exon.genePredExt \
--iso_anno PAusage.${spe}.bed6+ \
--t_e_out terminal_exon_annotation.${spe}.txt \
--work_path /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff
done

for spe in m r;do
python /home/user/data2/lit/project/ZNF271/02-APA-1/run/py/terminal_exon_annotation.loose_for_m_r.py \
--t_e ${spe}.terminal_exon.genePredExt \
--iso_anno PAusage.${spe}.bed6+ \
--t_e_out terminal_exon_annotation.${spe}.txt \
--work_path /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff
done

less cds_type_II.gene.gene_name.r.lst | sed '1i ZNF271' > cds_type_II.gene.gene_name.r.add_nc.lst
less cds_type_II.gene.gene_name.m.lst | sed '1i Zfp35' > cds_type_II.gene.gene_name.m.add_nc.lst

grep -f cds_type_II.gene.gene_name.m.add_nc.lst -w PAusage.m.bed6+ | awk '$5>=2'> PAusage.type_II.m.bed6+

less terminal_exon_annotation.m.txt | \
awk -v OFS='\t' '{if($3=="+")print $2,$4,$7,$1;else print $2,$6-1,$5,$1}' > terminal_exon.m.bed

bedtools intersect -a terminal_exon.m.bed -b PAusage.type_II.m.bed6+ -wb | cut -f 5-10 > PAusage.type_II.fil.m.bed6+

less gencode.vM23.annotation.genePredExt | cut -f 2,3,4,5,6,7,12,1 | awk '{print $0}' > t_id_chr_str_t_s_e_cds_s_e_gene_name.m.txt

# download from biomart
# cds_l_transcript_l_ensembl_98_m.txt

# gff to gtf 
gffread type_II.r.gtf -o- > type_II.r.gff3

# rheMac8 fa
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/rheMac8.fa.gz

# 
gunzip -c rheMac8.fa.gz > rheMac8.fa
gffread -x - -g  rheMac8.fa type_II.r.gff3 > type_II.cds.nucl.fa
seqkit fx2tab -l type_II.cds.nucl.fa | sed 's/gene=//g' | tr ' ' '\t'|awk -v OFS='\t' '{print $2,$1,$4,$4/3}' > type_II.cds_l_transcript_l.r.txt


grep -f cds_type_II.gene.gene_name.r.add_nc.lst -w PAusage.r.bed6+ | awk '$5>=2'> PAusage.type_II.r.bed6+

less terminal_exon_annotation.r.txt | \
awk -v OFS='\t' '{if($3=="+")print $2,$4,$7,$1;else print $2,$6-1,$5,$1}' > terminal_exon.r.bed

bedtools intersect -a terminal_exon.r.bed -b PAusage.type_II.r.bed6+ -wb | cut -f 5-10 > PAusage.type_II.fil.r.bed6+

less Final.426909.withName.gpe | cut -f 2,3,4,5,6,7,12,1 | awk '{print $0}' > t_id_chr_str_t_s_e_cds_s_e_gene_name.r.txt


