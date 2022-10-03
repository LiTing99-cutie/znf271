#!/bin/sh

if [ $1 == "-h" ] || [ $# -eq 0 ];then
    echo "This script will output the union peaks of two replications of macs2 peaks.
Usage: sh $0 [broad|narrow] rep1.Peak rep2.Peak outFile"
exit 0
fi

if [ $1 == "broad" ];then
    bedtools intersect -wo -a $2 -b $3 |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$11){print $2"\t";}else{print $11"\t";}if($3<=$12){print $12"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge -i stdin >$4
fi
if [ $1 == "narrow" ];then
    bedtools intersect -wo -a $2 -b $3 |awk -v OFS="\t" -v ORS="" '{print $1"\t";if($2<=$12){print $2"\t";}else{print $12"\t";}if($3<=$13){print $13"\n"}else{print $3"\n"}}'|sort|uniq |sort -k1,1 -k2,2n|bedtools merge -i stdin >$4
fi
