#!/usr/bin/sh

# terminal_exon_annotation.txt
terminal_exon_annotation=$1
dir=$(dirname $terminal_exon_annotation)
prefix=$(basename $terminal_exon_annotation | sed 's/.txt//g')
less $terminal_exon_annotation| awk -v OFS='\t' '{print $2,"custom","transcript",$4+1,$5,".",$3,".","gene_id \""$1"_proximal\"; transcript_id \""$1"_proximal\"; gene_name \""$1"\""}' \
> $dir/$prefix.proximal.gtf
less $terminal_exon_annotation|awk -v OFS='\t' '{print $2,"custom","transcript",$6+1,$7,".",$3,".","gene_id \""$1"_distal\"; transcript_id \""$1"_distal\"; gene_name \""$1"\""}' \
> $dir/$prefix.distal.gtf
cat $dir/$prefix.proximal.gtf $dir/$prefix.distal.gtf > $dir/$prefix.gtf
rm $dir/$prefix.proximal.gtf $dir/$prefix.distal.gtf