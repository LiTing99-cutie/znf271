

# pas supported samples
cd /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/typeII-III/pas_supported_sample
mkdir pas_supported_sample
cd pas_supported_sample
mkdir nine_samples

# proximal PAS location
less ../p_e_r_c.txt |awk '{if($2-30>=0)print $1,$2-30,$3+30,$4,$5,$6,$7;else print $1,1,$3,$4,$5,$6,$7}' > p_e_r_c_updown30.txt
[ -f map_2_sample.txt ] && rm -rf map_2_sample.txt

# proximal PAS location map to nine samples
for file in ../nine_samples/ENCFF*.PAusage.bed6+;do
dir=$(dirname $file)
sample=$(basename -s .PAusage.bed6+ $file)
while read gene;do grep "[[:space:]]$gene[[:space:]]" $file;done < ../gene.lst > nine_samples/$sample.type_II_III.txt
bedtools intersect -a p_e_r_c_updown30.txt -b nine_samples/$sample.type_II_III.txt -wa -wb -s | awk '$4==$11' | awk '{print $0"\t""'$sample'"}' >> map_2_sample.txt
done

# count supported samples for each gene
awk '{count[$4]++;} END {for(i in count) {print i"\t"count[i]}}' map_2_sample.sorted.txt > supported_sample.txt
# 35
awk '$12>1' map_2_sample.sorted.txt| awk '{count[$4]++;} END {for(i in count) {print i"\t"count[i]}}' > supported_high_confidence_sample.txt
join -t $'\t' -1 1 -2 1 -a1 -a2 <(sort -k1,1 supported_sample.txt) <(sort -k1,1 supported_high_confidence_sample.txt) > supported_sample.final.txt

# validated 
# 61
awk '{x[$4] += $12} END {for(i in x){print i, x[i]}}' map_2_sample.sorted.txt | sort -k1,1 > map.g_cnt.txt
join -t $'\t' -1 1 -2 4 map.g_cnt.txt <(sort -k4,4 ../p_e_r_c.txt) | awk '$2==$6'|wc -l
