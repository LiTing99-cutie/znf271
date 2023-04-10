#!/usr/bin/sh

################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 13 Mar 2023 08:33:43 PM CST
################################################

set -eou pipefail

nohup bash run.dog.sh &>log/dog.log &
