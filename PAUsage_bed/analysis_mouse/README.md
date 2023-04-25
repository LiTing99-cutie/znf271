
# cp scripts
# 2023/4/25
mkdir scripts
cp ../scripts/Step* scripts
cp ../scripts/*sh scripts
cp /home/user/data2/lit/project/ZNF271/02-APA-1/bin/stringtie.universal.sh scripts
cp ../scripts/get_fasta.ipynb scripts
cp ../scripts/cds_p_to_genomic_p.ipynb scripts
jupyter-nbconvert --to python scripts/get_fasta.ipynb scripts/get_fasta.py
jupyter-nbconvert --to python scripts/cds_p_to_genomic_p.ipynb scripts/cds_p_to_genomic_p.py

# run
gtf=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.gtf
gpe=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/mouse/gencode.vM23.basic.annotation.genePredExt
md=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/mouse.metadata.clean.txt
frag=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/mouse/fragmentation.score.clean.txt
uniq_bam_path=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/mouse
archive=sep2019
species=mus_musculus
case=Zfp35

bash scripts/Step_1_te_backbone.sh $gpe $gtf

bash scripts/Step_1_1_cutoff_filter_compare.sh PAusage_bed/PAusage.bed6+ output PAusage_bed no

bash scripts/Step_2_pas_map_to_backbone.sh PAusage_bed/PAusage.fil.bed6+ output/do_te_fil.bed output

nohup bash scripts/Step_3_RNA_seq.sh \
$md \
$frag \
$uniq_bam_path &>log/Step_3_RNA_seq.log &

## wrong so rerun
	nohup Rscript ../scripts/wilcox_test.R "output/stringtie/stringtie.rpkm.txt" \
	"/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/mouse.metadata.clean.txt" \
	"/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/mouse/fragmentation.score.clean.txt" \
	"output" &

nohup bash -c "time bash scripts/Step_3_1_disrupt_cds.sh $gtf $gpe $archive $species no" &> log/Step_3_1_disrupt_cds.log &

bash scripts/Step_4_case.sh $md $frag $case