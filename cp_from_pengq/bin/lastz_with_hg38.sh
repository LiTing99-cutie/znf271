################################################
#File Name: lastz_with_BGMv1.sh
#Author: Wanqiu Ding    
#Mail: wanqiuding@163.com
#Created Time: Tue Aug  7 18:35:34 2018
################################################

#!/bin/sh 

BGMDir=/rd1/user/dingwq/rhesus_genome/newGenome/liftOver_chain/rheMac8_chr
queryDir=/rd1/user/dingwq/rhesus_genome/newGenome/liftOver_chain/hg38
genome="hg38"
cat hg38_vs_rheMac8_chrPair.txt | while read queryNum BGMNum
do
lastz $BGMDir/chr${BGMNum}.fa $queryDir/chr${queryNum}.fa --chain --gapped --gfextend --step=20 --format=maf --rdotplot=rheMac8_chr${BGMNum}_${genome}_chr${queryNum}.rdotplot > rheMac8_chr${BGMNum}_${genome}_chr${queryNum}.maf 2>lastz_hg38_chr${queryNum}.err &

grep -v 'NA' rheMac8_chr${BGMNum}_${genome}_chr${queryNum}.rdotplot|sed '1d' |sed "1i rheMac8_chr${BGMNum}\t${genome}_chr${queryNum}" >rheMac8_chr${BGMNum}_${genome}_chr${queryNum}.remNA.rdotplot
Rscript /rd1/brick/dingwq/bin/rdotplot.R -p="rheMac8_chr${BGMNum}_${genome}_chr${queryNum}.remNA.rdotplot.pdf" -m="rheMac8_chr${BGMNum} vs ${genome}_chr${queryNum}" -x="rheMac8_chr${BGMNum}" -y="${genome}_chr${queryNum}" <rheMac8_chr${BGMNum}_${genome}_chr${queryNum}.remNA.rdotplot

done
