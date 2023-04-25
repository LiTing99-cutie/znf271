#! /bin/sh

script_path=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/analysis_mouse/scripts
pro_dis_bin=output/pro_dis_bin.txt
map=output/transcript_id_chr_strand_symbol.txt
output_path=output/predict
archive=$1
species=$2

source activate base

# get gene_id to predict
Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/scripts/get_lncRNA_gene_id.R \
$pro_dis_bin \
$map \
$output_path

pushd $output_path

### get fasta ###
python3 $script_path/get_fasta.py --input to_predict_gene_id.txt --output cdna --archive $archive --species $species
#

### get orf ###
cat cdna/*fa > lncRNA_cdna.fa
# get orf 
getorf -find 1 -minsize 300 -sequence lncRNA_cdna.fa -outseq lncRNA_orf_prot.fa
# get orf length and for each transcript, the longest orf is kept
paste -d '\t' <(seqkit fx2tab -l lncRNA_orf_prot.fa | cut -f 1 | cut -d ' ' -f 1) <(seqkit fx2tab -l lncRNA_orf_prot.fa | cut -f 4) > lncRNA_orf_prot_len.txt
# output lncRNA_orf_prot_len.filter.txt transcript_id orf_id orf_length
Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/bin/filter_longest_orf.R lncRNA_orf_prot_len.txt
# get cdna coordinate
less lncRNA_orf_prot.fa | grep '>' | sed 's/>//g;s/\[//g;s/\]//g;s/-//g' | sed 's/(REVERSE SENSE)/-/g'|awk -v OFS='\t' '{if ($4=="-")print $1,$3,$2,$4 ;else print $1,$2,$3,"+"}'|sort -k1,1>lncRNA_orf_prot_co.txt
join -1 2 -2 1 <(sort -k1,1 -u lncRNA_orf_prot_len_filter.txt|sort -k2,2) lncRNA_orf_prot_co.txt | awk -F " " -v OFS='\t' '{split($1,x,".");print x[1],$2,$4,$5,$3,$6}'>lncRNA_orf_prot_co_filter.txt
#	
# get genomic position	
python3 $script_path/cds_p_to_genomic_p.py --archive $archive --species $species
#

# output
popd
Rscript ../scripts/lncRNA.cds.R "output/predict/lncRNA_orf_prot_co_filter_genomic.txt" "$map" "output"
#