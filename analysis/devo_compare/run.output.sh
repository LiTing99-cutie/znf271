#!/usr/bin/sh

################################################
#File Name: run.output.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 17 Feb 2023 08:36:05 PM CST
################################################

set -eou pipefail

Rscript script/wilcox_test.R mouse Zfp35 > output/mouse/stat.txt;Rscript script/wilcox_test.R rhesus ZNF271 > output/rhesus/stat.txt

Rscript script/develop_case.R mouse Zfp35 ; Rscript script/develop_case.R rhesus ZNF271

Rscript /home/user/data2/lit/project/ZNF271/02-APA-1/bin/develop_case.R ZNF271P
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ZNF271P/rpkm.pro.2.pdf output/human/rpkm.pro.2.pdf
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ZNF271P/rpkm.dis.2.pdf output/human/rpkm.dis.2.pdf
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ZNF271P/pau.2.pdf output/human/pau.2.pdf

file=/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.txt
(head -n1 $file; grep ZNF271P $file) > output/human/stat.txt

# 2023/3/27 rerun; previous script use wrong rpkm.txt for human
rm -rf human/analysis/stringtie_whole_genome_ref_based_all_1_loose/stringtie.rpkm.txt
ln -s /home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1_loose/stringtie.rpkm.txt human/analysis/stringtie_whole_genome_ref_based_all_1_loose/
Rscript script/wilcox_test.R human ZNF271P > output/human/stat.txt
Rscript script/develop_case.R human ZNF271P

# 2023/3/27 dog to dog/stat.txt
Rscript script/wilcox_test_dog.R dog ENSCAFG00845006933 >output/dog/stat.txt