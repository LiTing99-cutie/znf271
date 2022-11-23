#!/usr/bin/sh

################################################
#File Name: /home/user/data2/lit/project/ZNF271/02-APA/run/sh/ref_based_develop_diff_pa_usage.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 26 Oct 2022 09:29:44 PM CST
################################################

set -eou pipefail

# get processed terminal exon annoataion
python /home/user/data2/lit/project/ZNF271/02-APA/run/py/terminal_exon_annotation.py \
--t_e annotation/terminal_exon/gencode.v41.annotation.terminal_exon.genePredExt \
--iso_anno annotation/terminal_exon/ENCFF678TEN.PAusage.bed6+ \
--t_e_out annotation/terminal_exon/terminal_exon_annotation.ref_based.txt

# format to gtf 
bash /home/user/data2/lit/project/ZNF271/02-APA/bin/proximal_distal_gtf.sh annotation/terminal_exon/terminal_exon_annotation.ref_based.txt

# stringtie call rpkm
nohup bash -c "time bash bin/stringtie.universal.sh annotation/terminal_exon/terminal_exon_annotation.ref_based.gtf \
/home/user/data2/lit/project/ZNF271/02-APA/analysis/stringtie_whole_genome_ref_based \
frag_filter" &>log/stringtie_whole_genome_ref_based.log &

wait

# wilcox test
nohup bash -c "time python run/py/wilcox_test.py --rpkm /home/user/data2/lit/project/ZNF271/02-APA/analysis/stringtie_whole_genome_ref_based/stringtie.rpkm.txt \
--output /home/user/data2/lit/project/ZNF271/02-APA/output/final_res_ref_based.txt" &>log/wilcox.ref_based.log &

wait

python /home/user/data2/lit/project/ZNF271/02-APA/bin/add_gene_name.py \
--vs annotation/ensembl_gene_id_vs_symbol.txt \
--file /home/user/data2/lit/project/ZNF271/02-APA/output/final_res_ref_based.txt

cat /home/user/data2/lit/project/ZNF271/02-APA/output/final_res_ref_based.with_geneName.txt | \
awk '$3<0.001 && $5-$4>0.1' > /home/user/data2/lit/project/ZNF271/02-APA/output/diff_develop/sig.ref_based.ba.lst

cat /home/user/data2/lit/project/ZNF271/02-APA/output/final_res_ref_based.with_geneName.txt | \
awk '$3<0.001 && $5-$4>0.1' > /home/user/data2/lit/project/ZNF271/02-APA/output/diff_develop/sig.ref_based.ab.lst