#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 28 Nov 2022 07:41:18 PM CST
################################################

set -eou pipefail

chrom_size_h=/home/user/data/lit/database/public/genome/hg38/hg38.chrom.sizes
chrom_size_m=/home/user/data/lit/database/public/genome/mouse/mm10.chrom.sizes

for spe in h;do
cat processed.$spe.bed12+| cut -f 1-12 | awk 'length($1)<6' | LC_COLLATE=C sort -k1,1 -k2,2n > processed.$spe.sorted.comChr.bed12

bedToBigBed -type=bed12 processed.$spe.sorted.comChr.bed12 $chrom_size_h processed.$spe.sorted.bb

[ -f /home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/hg38/iso_seq/processed.$spe.sorted.bb ] && rm -rf /home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/hg38/iso_seq/processed.$spe.sorted.bb
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/data/iso_seq_reads/processed.$spe.sorted.bb \
/home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/hg38/iso_seq/processed.$spe.sorted.bb    
done

for spe in m;do
cat processed.$spe.bed12+| cut -f 1-12 | awk 'length($1)<6' | LC_COLLATE=C sort -k1,1 -k2,2n > processed.$spe.sorted.comChr.bed12

bedToBigBed -type=bed12 processed.$spe.sorted.comChr.bed12 $chrom_size_m processed.$spe.sorted.bb

[ -f /home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/hg38/iso_seq/processed.$spe.sorted.bb ] && rm -rf /home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/hg38/iso_seq/processed.$spe.sorted.bb
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/data/iso_seq_reads/processed.$spe.sorted.bb \
/home/user/data2/rbase/ucsc2/htdocs/data/lit/myHub/mm10/iso_seq/processed.$spe.sorted.bb
done