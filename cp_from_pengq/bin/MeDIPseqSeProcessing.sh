#!/bin/sh
#######################
###Date: 2017/5/7######
###Name: Ding Wanqiu###
#######################

<<!
if [ $1 = -h ] || [ $# -eq 0 ];then
    echo "This script will map the single-end fastq with bwa mem to genome and call 6mA peaks with MACS2.
    Please make sure the following programs are in your PATH:
    bwa		for mapping
    MACS2	for peak calling
    bamtools & samtools		for bam/sam file analysis
Usage: sh $0 -p threads -i bwa.db.prefix -g chrSize -d data.dir -c sample.list"
exit 0
fi
!
usage(){
    cat << EOF
Description: Perform paired-end reads mapping with bwa mem to genome and peak calling with MACS2.
Dependencies: bwa, samtools, picard, bamtools and MACS2.
Options:
    p	INT	number of threads
    i	FILE	index prefix of BWA [which can be built through 'bwa index -a is genome.fa' when the genome is short]
    g	FILE	genome size [chrom.sizes is a two-column file/URL: <chromosome name> <size in bases> which can be built though ]
    s   INT     effective genome size for macs2 callpeak
    d	DIR	the directory of the fastq data files
    c	FILE	the list of samples which will be analyzed in the -d directory
Usage: sh $0 -p threads -i bwa.db.prefix -g chrSize -d data.dir -c sample.list"
sample.list is a list where every line is a sample prefix such as 
[the directory name is the same as sample prefix which fastq file must include this prefix xx.1.clean.fq.gz]:
PC30_LY_NT-1-input
PC30_LY_NT-1-IP
PC30_LY_NT-2-input
PC30_LY_NT-2-IP
EOF
    exit 0
}
[ $1 ] || usage

sample=sample.conf
threads=5
while getopts "hp:i:g:s:d:c:" OPTION
do
    case $OPTION in
        h) usage;;
        p) threads=$OPTARG;;
        i) bwaPrefix=$OPTARG;;
        g) chrSize=$OPTARG;;
	s) genomeSize=$OPTARG;;
        d) dataDir=$OPTARG;;
        c) sample=$OPTARG;;
        ?) usage
	   exit 0
	   ;;
    esac
done

dataDir=realpath $dataDir

if [ -d log ];then
    [ -d log_old ] && rm -r log_old
    mv log log_old
fi
mkdir -p log

while read line
do
    echo "Started on $(date)" >>log/run.log
    #echo "-----$(date):Quality control (cutadaptor, filter short length reads) in $line-----" >>log/run.log
#	cutadapt -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $dataDir/$line/trimmed.${line}.1.clean.fq.gz -p $dataDir/$line/trimmed.${line}.2.clean.fq.gz $dataDir/$line/${line}.1.clean.fq.gz $dataDir/$line/${line}.2.clean.fq.gz >$dataDir/$line/trimmed.${line}.cutadp || exit $?
#	cutadapt --minimum-length=36 -o $dataDir/$line/trimmedL36.${line}.1.clean.fq.gz -p $dataDir/$line/trimmedL36.${line}.2.clean.fq.gz $dataDir/$line/trimmed.${line}.1.clean.fq.gz $dataDir/$line/trimmed.${line}.2.clean.fq.gz >$dataDir/$line/trimmedL36.${line}.cutadpL36 || exit $?

	echo "-----$(date):Running bwa mem in single end mode in $line-----" >>log/run.log
    bwa mem -t $threads $bwaPrefix $dataDir/${line}.fastq 2>log/${line}.bwa.mem.log|samtools view -h - |samtools view -bSh >${line}.se.bam || exit $?
	
    echo "-----$(date):Extract unique mapped reads, sorted reads and remDuplicated reads in $line-----" >>log/run.log
	samtools view -h -q 60 -F 0x800 ${line}.se.bam | samtools view -bSuh - | samtools sort -@ $threads - -o ${line}.uniq.sorted.bam || exit $? 
	#samtools view -h ${line}.aln.pe.bam|grep -v 'XA'|samtools view -bSuh - 2>log/samtoolsView.log|samtools sort -@ threads - -o ${line}.uniq.sorted.bam || exit $?
    #replace the "-v XA" unique mapping parameter with -q 60 and -F 0x800 to get high quality mapped reads with the highest mapping quality and to remove the supplymentary mapped reads
	samtools index ${line}.uniq.sorted.bam || exit $?
    echo "$(date):Finish mapping step and successfully index bam" >>log/run.log	

    java -Xmx16g -XX:ParallelGCThreads=$threads -jar /rd1/user/liym/tools/picard-2.17.6/picard.jar MarkDuplicates I=${line}.uniq.sorted.bam O=${line}.uniq.sorted.markDup.bam M=${line}.uniq.sorted.markDup.metrix.txt REMOVE_DUPLICATES=true 2>log/${line}.uniq.sorted.markDup.log || exit $?
    samtools index ${line}.uniq.sorted.markDup.bam || exit $?
    echo "$(date):Successfully finish the markDuplicate and index step" >>log/run.log

    echo "-----$(date):Run bamStats for reads in $line-----" >>log/run.log
    bamtools stats -in ${line}.se.bam -insert >${line}.se.bamStats.txt || exit $? &
    samtools flagstat ${line}.se.bam >${line}.se.flagstat.txt || exit $? &
    bamtools filter -isMapped false -in ${line}.se.bam >${line}.unmapped.bam || exit $? &
    bamtools stats -in ${line}.uniq.sorted.bam -insert > ${line}.uniq.sorted.bamStats.txt || exit $? &
    bamtools stats -in ${line}.uniq.sorted.markDup.bam -insert > ${line}.uniq.sorted.markDup.bamStats.txt || exit $? &
        
    echo "-----$(date):Convert bam to bigwig file in $line-----" >>log/run.log
    bamCoverage --bam ${line}.uniq.sorted.markDup.bam --outFileFormat bedgraph --outFileName ${line}.uniq.sorted.markDup.normRPKM.bdg --normalizeUsing RPKM --extendReads -e 100 --binSize 10 2>log/${line}.bamCoverage.log || exit $?
    bedGraphToBigWig ${line}.uniq.sorted.markDup.normRPKM.bdg $chrSize ${line}.uniq.sorted.markDup.normRPKM.bw || exit $? &
done <$sample

resultPath=$(pwd)
mkdir -p macs2_results
peakDir=$resultPath/macs2_results
ls *input.uniq.sorted.markDup.bam | sed 's/input.uniq.sorted.markDup.bam//' | while read label
    do
        echo "-----$(date):Call peaks by MACS2 in $label-----" >>log/run.log
	macs2 callpeak -t ${label}IP.uniq.sorted.markDup.bam -c ${label}input.uniq.sorted.markDup.bam -n $label --outdir $peakDir -f BAM --bdg -g $genomeSize --fix-bimodal --extsize 100 2> log/macs2_${label}.log || exit $? &

	echo "-----$(date):bamCompare IP to input in read coverage in $label-----" >>log/run.log
	bamCompare -b1 ${label}IP.uniq.sorted.markDup.bam -b2 ${label}input.uniq.sorted.markDup.bam --operation subtract --extendReads 100 --binSize 10 --scaleFactorsMethod None --normalizeUsing RPKM --outFileFormat bedgraph -o ${label}subtract.bdg 2>log/${label}bamCompare.log || exit $?
	bedGraphToBigWig ${label}subtract.bdg $chrSize ${label}subtract.bw || exit $?

	echo "Finished on $(date)" >>log/run.log
    done
