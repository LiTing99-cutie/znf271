#!/usr/bin/sh

################################################
#File Name: scripts/visual_normalized_bam_m_s.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 24 Apr 2023 09:15:30 PM CST
################################################

set -eo pipefail



pushd analysis_mouse
[ -d log ] || mkdir log
bash -c "time bash ../scripts/visual_normalized_bam.sh \
/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/mouse \
output/visual_n_bam \
/home/user/data/lit/project/ZNF271/visual_n_bam/mouse \
/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/mouse" &>log/visual_normalized_bam_mouse.log
popd

pushd analysis_rhesus
[ -d log ] || mkdir log
bash -c "time bash ../scripts/visual_normalized_bam.sh \
/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/uniq/r10/ \
output/visual_n_bam \
/home/user/data/lit/project/ZNF271/visual_n_bam/r10 \
/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/r10" &>log/visual_normalized_bam_r10.log
popd