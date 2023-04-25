#!/usr/bin/sh

gtf=$1
output_path=$2

less $gtf | awk -v OFS='\t' '$3=="exon" {print $1,$4,$5,$12,$16,$7}'|sed 's/"//g;s/;//g'|awk -v OFS='\t' '{print $1,$2-1,$3,$4,$5,$6}' > $output_path/exon.bed6
