#!/bin/bash

if [ $1 == "-h" ] || [ $# -eq 0 ];then
    echo "This script will map the input fastq with bwa mem and call nucleosome occupancy with DANPOS.
    Please make sure the following programs are in your PATH:
    bwa        mapping
    DANPOS2    nucleosome calling
Usage: sh $0 shreads readLength bwa.db.prefix <read1.fq> <read2.fq>"
exit 0
fi

threads=$1
readL=$2
dbPrefix=$3
chrSize=$4
read1=$5
read2=$6

if [ -d log ];then
    [ -d log_old ] && rm -r log_old
    mv log log_old
fi
mkdir log

if [ $readL -lt 70 ];then
    echo "Running bwa aln------"
    bwa aln -t $threads $dbPrefix $read1 >out1.sai 2>log/aln_1.log &
    bwa aln -t $threads $dbPrefix $read2 >out2.sai 2>log/aln_2.log
    echo "Running bwa sampe------"
    bwa sampe out1.sai out2.sai $read1 $read2 2>log/sampe.log |samtools view -bSu - 2>log/samtoolsView.log| samtools sort -@ $threads -n - out.name.sorted
else
    echo "Running bwa mem------"
    bwa mem -t $threads $dbPrefix $read1 $read2 |samtools view -bSu - 2>log/samtoolsView.log|samtools sort - -@ $threads -n out.name.sorted ;
fi
echo "Run bamstats for all reads-----"
bamtools stats -in out.name.sorted.bam >out.name.sorted.bamStats ;
echo "Run filter(properPair,insertSize[50,300],mapQuality>=20) for pair-end reads-----"
mismatch=$(($readL/10));
perl /home/pengq/bin/bam_insert_filter.pl -i out.name.sorted.bam -q 20 -n $mismatch | samtools view -bSu - | samtools sort - -@ $threads myFiltered.sorted ;
echo "Running bamstats for filtered reads-----";
bamtools stats -in myFiltered.sorted.bam >myFiltered.sorted.bamStats ;

echo "Call nucleosomes by danpos.py dpos----"
danpos.py dpos myFiltered.sorted.bam -jd 147 -u 1 -m 1 --extend 74 -c 10000000 -o danpos_rst >log/danpos.log 2>log/danpos.err;
cd ./danpos_rst/pooled
ls *.wig |sed 's/.wig//g' | while read file;do wigToBigWig -clip ${file}.wig $chrSize ${file}.bw 2>wigTobw.log;done
