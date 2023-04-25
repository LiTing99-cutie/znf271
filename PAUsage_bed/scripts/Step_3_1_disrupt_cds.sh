#!/usr/bin/sh

################################################
#File Name: scripts/Step_3_1_disrupt_cds.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 21 Apr 2023 11:30:35 PM CST
################################################

set -eo pipefail


## Predict lncRNA cds region
bash scripts/lncRNA.sh
## Annotate proximal PAS
Rscript scripts/cds_anno.R
## Proximal PAS supported
bash scripts/proximal_pas_supported.sh