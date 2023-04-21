

PAS_iso=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/human/PAusage.fil.bed6+
PAS_DaPars2=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/benchmark/dapars/Dapars2_test_chr/Dapars2_result_temp.chr.bed6
PAS_3_atlas=/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/benchmark/atlas.pos.bed
output_path=/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output
# function
select_com_chr(){
    grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]'
}
updown15(){
awk -v OFS='\t' '{if($2-15>=0)print $1,$2-15,$3+15,$4,$5,$6;else print $1,1,$3,$4,$5,$6}'|sort -k1,1 -k2,2n
}
inter(){
	bedtools intersect -a $1 -b $2 -wa -s| wc -l
}


# output
pushd $output_path
less $PAS_iso|select_com_chr|updown15 > PAS_iso.bed
less $PAS_DaPars2|select_com_chr|updown15 > PAS_DaPars2.bed
less $PAS_3_atlas|awk -v OFS='\t' '{print "chr"$1,$2,$3,$4,$5,$6}'|select_com_chr|updown15 > PAS_3_atlas.bed

inter_2_1=$(inter PAS_iso.bed PAS_3_atlas.bed)
inter_2_3=$(inter PAS_iso.bed PAS_DaPars2.bed)
n_2=$(wc -l PAS_iso.bed|cut -f1 -d ' ')
echo -e "inter_2_1\tinter_2_3\tn_2" > intersect.txt
echo -e "$inter_2_1\t$inter_2_3\t$n_2" >> intersect.txt

Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/scripts/compare_3_seq_dapars2.R intersect.txt $PWD

popd