#!/usr/bin/sh

gtf=/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gtf
output_path=output

less $gtf | awk -v OFS='\t' '$3=="exon" {print $1,$4,$5,$12,$16,$7}'|sed 's/"//g;s/;//g'|awk -v OFS='\t' '{print $1,$2-1,$3,$4,$5,$6}' > $output_path/exon.bed6
