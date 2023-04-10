#!/usr/bin/sh
set -eou pipefail


echo "Generate annotation and call rpkm of common and extended region [Figure 1A left and middle]" 
bash -c "time bash bin/ref_based_develop_diff_pa_usage.universal.sh \
ref_based_all_1_loose \
annotation/terminal_exon/PAusage.h.c_anno.bed6+ \
/home/user/data2/lit/project/ZNF271/02-APA-1 \
bin/terminal_exon_annotation.loose.py" &> log/ref_based_all_1_loose.log

echo "Plot ScatterPlot and remove outliers [Figure 1B]"
Rscript run/R/ScatterPlot_c.R

echo "Add ensembl gene id for further analysis"
python bin/add_gene_name.py \
--vs annotation/map/ensembl_gene_id_type_symbol.txt \
--file output/final_list/res.rm_outlier.txt

echo "GO analysis [Figure 1D and E]"
Rscript run/R/GO_c.R

echo "lncRNA orf prefict"
pushd /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/orf_predict/loose_n_l/
bash /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/orf_predict/loose_n_l/run.sh 
popd

echo "Proximal proximal PA effect [Figure 1A right]" 
Rscript run/R/cds_c.R

echo "Plot pie plot [Figure 1C]"
Rscript run/R/pie_multiple_circles_c.R

echo "Prepare for heatmap"
pushd /home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean
bash run.sh 
popd

# Picture different when run in command line or in R studio
echo "Plot heatmap [Figure 1F]"
Rscript run/R/cds_compara.R


python bin/terminal_exon_annotation.loose.py \
--t_e annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt \
--iso_anno annotation/terminal_exon/PAusage.h.c_anno.bed6+ \
--t_e_out annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.txt \
--work_path /home/user/data2/lit/project/ZNF271/02-APA-1