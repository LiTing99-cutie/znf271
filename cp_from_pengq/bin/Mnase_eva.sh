#!/bin/sh

usage(){
  cat << EOF
Description:
    The tool wraps steps of reads mapping,(Illumina sequence reads up to 100bp) a set of post-mapping evaluations and nucleosome position calling:
    
    Step		        Program
    mapping:			bwa
    post-mapping evaluations:	samtools, bamtools,region_dinucleotideFreq.pl,nucleosome_coverage.pl,TSS_nucleosomeCov.R,RSeQC
    nucleosome calling:         DANPOS2
    
    So, please make sure above programs are in your PATH.
    Before running this tool, you must export the arguments needed by bwa (bwaAlnArgs,bwaSamseArgs/bwaSampeArgs) and DANPOS2(danpos2Args).

Usage:
    export bwaAlnArgs='All bwa aln arguments' 
    export bwaSamseArgs/bwaSampeArgs='All bwa samse/sampe arguments'
    export bwaMemArgs='All bwa mem arguments' [You should set one of bwaAlnArgs,bwaSampeArgs/bwaSampeArgs or bwaMemArgs based on your need.]
    export danpos2Args='All danpos2 arguments' [Not include the -o argument]
    Note: above arguments must be exported into $(basename $0) using export command
    $(basename $0) Options >$(basename $0 .sh).log 2>$(basename $0 .sh).err
	
Example:
    export bwaAlnArgs='-l 25 -k 2 -t 5 -I /rd1/user/data/genome/bwa/hg19/bwtsw' 
    export bwaSamseArgs/bwaSampeArgs='-n 3 /rd1/user/data/genome/bwa/hg19/bwtsw'
    export bwaMemArgs='-t 5 /rd1/user/data/genome/bwa/hg19/bwtsw'
    export danpos2Args='-jd -147 -u 1 --extend 74' 
    # for single end reads
    $(basename $0) -f /data/genome/fna/rheMac2/all.fa read.fastq >$(basename $0 .sh).log 2>$(basename $0 .sh).err
    
    # for pair end reads
    $(basename $0) -f /data/genome/fna/rheMac2/all.fa read_1.fastq read_2.fastq >$(basename $0 .sh).log 2>$(basename $0 .sh).err

Options:
    -o  DIR     Output directory for all results[Current Directory]
    -e  LOGIC   Run post-mapping evaluation and its following steps
    -f  FILE    Reference sequence in fasta format
    -y  LOGIC   The header of reference genome is ">chr*" format.[default:">*"]
    -g  FILE    A gene structure file in gpe format(no bin column)
    -c  FILE    fastqc result directory for reads(If the reads are pair-end,seperate the two directories by comma)
    -l  INT     Read length
    -1  FILE    read1 in fastq format
    -2  FILE    read2 in fastq format if your data is paired end
    -h          Show this help information
EOF
    exit 0
}
[ $1 ] || usage

inputDir=$PWD
outputDir=$PWD
perlScript=~/bin

while getopts "ho:ef:sg:yl:1:2:cd:" OPTION
do
    case $OPTION in
        h) usage;;
        o) outputDir=$OPTARG;;
        e) evaluation=1;;
        f) refFa=$OPTARG;;
        s) skip=$OPTARG;;
        g) gpeStructure=$OPTARG;;
        y) chrom=$OPTARG;;
        c) fastqcDirs=$OPTARG;;
        d) sample=$OPTARG;;      
        l) readLen=$OPTARG;;
        1) read1=$OPTARG;;
        2) read2=$OPTARG;;
        ?) usage;;
    esac
done

echo "--------------- Begin at $(date) ---------------"
shift $((OPTIND - 1))

if [ -d log ];then
    [ -d log_old ] && rm -r log_old
    mv log log_old
fi
mkdir log

#echo "Running bwa-----"
#if [ $read2 ];then
 #   if [ $bwaAlnArgs ];then
  #      echo "Running bwa aln $bwaSampeArgs $read1 >out1.sai 2>log/aln_1.log ..."
   #     bwa aln $bwaAlnArgs $read1 >out1.sai 2>log/aln_1.log &
    #    echo "Running bwa aln $bwaAlnArgs $read2 >out2.sai 2>log/aln_2.log ..."
     #   bwa aln $bwaAlnArgs $read2 >out2.sai 2>log/aln_2.log
      #  echo "Running bwa sampe $bwaSampeArgs out1.sai out2.sai $read1 $read2 2>log/sampe.log ..."
       # bwa sampe $bwaSampeArgs out1.sai out2.sai $read1 $read2 2>log/sampe.log |samtools view -bSu -| samtools sort -n - out.name.sorted
  #  else
   #     echo "Running bwa mem $bwaMemArgs $read1 $read2 2>log/mem.log ..."
    #    bwa mem $bwaMemArgs $read1 $read2 |samtools view -bSu -|samtools sort -n - out.name.sorted
   # fi    
#else
 #   if [ $bwaAlnArgs ];then
  #      echo "Running bwa aln $bwaAlnArgs $read1 $read1 >out.sai 2>log/aln.log ..."
   #     bwa aln $bwaAlnArgs $read1 >out.sai 2>log/aln.log
    #    echo "Running bwa samse $bwaSamseArgs out.sai $read1 2>log/samse.log ..."
     #   bwa samse $bwaSamseArgs out.sai $read1 2>log/samse.log |samtools view -bSu -| samtools sort -n - out.name.sorted
   # else
    #    echo "Running bwa mem $bwaMemArgs $read1 2>log/mem.log ..."
     #   bwa mem $bwaMemArgs $read1 |samtools view -bSu -|samtools sort -n - out.name.sorted
#    fi
 #fi

#if [ $evaluation ];then
  echo "Run bamstats for all reads-----"
   bamtools stats -in ../${sample}/out.name.sorted.bam >bamStats.txt ;
  #  echo "Run filter(properPair,insertSize[50,300],mapQuality>=30) for pair-end reads or (mapQuality>=30,unique) for single-end reads-----"
  #  if [ $read2 ];then
   #     bamtools filter -mapQuality ">=30" -tag "XT:U" -isProperPair true -in sorted.name.filtered.bam -out uniq.sorted.filtered.bam 2>log/bamtoolsFilter.log  
  #  else
   #     perl $perlScript/bam_insert_filter.pl -i sorted.name.filtered.bam | samtools view -bSu - | samtools sort -@ 3 - uniq.sorted.filtered
   # fi 
  #ORIGIN#  #samtools view -h -q 10 out.sorted.bam | awk '/^@/ || $12=="XT:A:U"' | samtools view -bSu - 2>log/samtoolsView.log | samtools sort - uniq.sorted 2> log/samtoolsSort.log
    bamtools stats -in out.sorted.uniq.bam >filteredbamStats.txt;
    
   # echo "Call nucleosomes by danpos.py dpos uniq.sorted.bam $danpos2Args -o danpos_rst >log/danpos.log 2>log/danpos.err ..."
   # danpos.py dpos uniq.sorted.filtered.bam $danpos2Args -o danpos_rst >log/danpos.log 2>log/danpos.err 
   # awk -v OFS="\t" '{if($1!="chr"){print $1,$4-1,$4}}' danpos_rst/uniq.sorted.filtered.smooth.positions.xls >danpos_rst/pooled/nucleosome_dyad.bed3 

    echo "Run post-mapping evaluation-----"
    if [ -d evaluation ];then
        [ -d evaluation_old ] && rm -r evaluation_old
        mv evaluation evaluation_old
    fi
    mkdir evaluation
   # echo "Calculate duplication level(fastqc)-----"
   # if [ $read2 ];then
    #    fastqc1=$(echo $fastqcDirs |cut -d "," -f1)
     #   fastqc2=$(echo $fastqcDirs |cut -d "," -f2)
        #TotalReads=$(grep "^Total Sequences" $fastqc1/fastqc_data.txt |awk '{print $3}')
        #TotalReads=$((${TotalReads}*2))
    #    echo -n "Read1 " >evaluation/DupLevel.txt 
     #   grep "^#Total Duplicate Percentage" $fastqc1/fastqc_data.txt |sed 's/#//' >>evaluation/DupLevel_fastqc.txt
      #  echo -n "Read2 " >>evaluation/DupLevel.txt
    #    grep "^#Total Duplicate Percentage" $fastqc2/fastqc_data.txt |sed 's/#//' >>evaluation/DupLevel_fastqc.txt
   # else
        #TotalReads=$(grep "^Total Sequences" $fastqcDirs/fastqc_data.txt |awk '{print $3}')
  #      grep "^#Total Duplicate Percentage" $fastqcDirs/fastqc_data.txt >evaluation/DupLevel_fastqc.txt
  #  fi
    echo "Calculate mapping rate-----"
    TotalReads=$(grep "^Total reads" bamStats.txt |awk '{print $3}')
    MappedR=$(grep "^Total reads" bamStats.txt |awk '{print $3}')
    uniqMappedR=$(grep "^Total reads" filteredbamStats.txt |awk '{print $3}')
    cd evaluation
    echo "Total reads: $TotalReads " >mappingRate.txt
    echo -n "Mapping Rate: " >>mappingRate.txt
    echo "scale=6;$MappedR/$TotalReads"|bc >>mappingRate.txt
    echo -n "Unique properly paired mapping rate: " >>mappingRate.txt
    echo "scale=6;$uniqMappedR/$TotalReads"|bc >>mappingRate.txt
    
    if [ -d log ];then
        [ -d log_old ] && rm -r log_old
        mv log log_old
    fi
    mkdir log
    echo "Fragment size distribution-----"
    java -jar /home/pengq/anaconda3/envs/python2/share/picard-2.26.10-0/picard.jar CollectInsertSizeMetrics -I ../out.sorted.uniq.rmdup.bam -O insert_size_metrics.txt -H insert_size_histogram.pdf > log/collectInsertS.log 2> log/collectInsertS.err
    #if [ $read2 ];then
    #    samtools rmdup ../uniq.sorted.filtered.bam uniq.sorted.rmdup.bam 2>log/uniq.sorted.rmdup.log
    #else
    #    samtools rmdup -s ../uniq.sorted.filtered.bam uniq.sorted.rmdup.bam 2>log/uniq.sorted.rmdup.log
    #fi
    #/rd1/user/liym/ToolKit/distogram_phasogram.py uniq.sorted.rmdup.bam -d -o distogram.pdf -x 0 300 2>log/distogram.log &
    
    echo "Calculate mismatch profile-----"
    mismatch_profile.py -i ../out.sorted.uniq.rmdup.bam -l 100 -o mismatchProfile 2>log/mismatchProfile.log &
    
    echo "Calculate dinucleotide frequency across nucleosomes-----"
    perl $perlScript/peReadsFragmentMid.pl -i ../out.sorted.uniq.rmdup.bam -f 147,147 >n147_mid.txt 
    perl $perlScript/region_dinucleotideFreq.pl -b n147_mid.txt -n AA,AT,TA,TT -u -150 -d 150 -f $refFa -c $chrom >AT_freq.tsv 2>log/AT_freq.log;
    perl $perlScript/region_dinucleotideFreq.pl -b n147_mid.txt -n GG,GC,CG,CC -u -150 -d 150 -f $refFa -c $chrom >GC_freq.tsv 2>log/GC_freq.log;
    if [[ -s AT_freq.tsv ]] && [[ -s GC_freq.tsv ]];then
        Rscript $perlScript/dinucleotide_plot.R AT_freq.tsv GC_freq.tsv dinucleotide.pdf 2>log/dinucleotide_R.log   
 fi
    
    echo "Calculate nucleosome profile across TSS----"
 #   danpos.py profile ../danpos_rst/pooled/uniq.sorted.filtered.smooth.wig --genefile_paths $gpeStructure --genomic_sites TSS --flank_up 1000 --flank_dn 1000 >log/profile.log 2>log/profile.err
    awk '{if($3 == "+" && $4>1000 ){print $2"\t"$4-1000"\t"$4+1001"\t"$2"_"$4"\t0\t"$3;}else if($3 == "-" && $4>1000 ){print $2"\t"$5-1001"\t"$5+1000"\t"$2"_"$5-1"\t0\t"$3}}' $gpeStructure |sort|uniq >TSS_updn1k.bed6
    perl $perlScript/nucleosome_coverage.pl -b TSS_updn1k.bed6 -i 6 -n ../danpos_rst/pooled/nucleosome_dyad.bed3 |sort -k1,1n >TSS_nucleosomeCov.tsv
    Rscript $perlScript/TSS_nucleosomeCov.R TSS_nucleosomeCov.tsv TSS_nucleosomeCov.pdf  2>log/TSS_nucleosomeCov.log
    cd ..
#fi

echo "--------------- End at $(date) ---------------"
