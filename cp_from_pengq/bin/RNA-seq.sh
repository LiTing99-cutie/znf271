#!/bin/bash
usage(){
  cat<<EOF
Description:
  This is for the RNA-Seq evaluation~
           Step		                                   Program
  post-mapping evaluations:   samtools, bamtools, Picard,mutationRate_tophat.pl, bedtools,
                              samSsEva.pl,samFlag.pl,coverage.R, mutRate_tophat.R,
                              gpeFeature.pl, gpeMerge.pl, gpe2bed.pl, fragmentation.py,
                              fragmentation.R and samContam.pl
  calculate RB Score:         calRnaSeqRbScore.pl
  calculate expression level: geneRPKM.pl, transRPKM.pl, geneFPKM.pl
  calculate coverage          wiggles
  Options:
    -a  FILE    Input fastqc_data.txt files (reads1)
    -b  FILE    Input fastqc_data.txt files (reads2)
    -i  FILE    Input uniq-sorted-bam files
    -o  DIR     Output directory for all results[Current Directory]
    -f  FILE    Reference sequence in fasta format
    -g  FILE    A gene structure file in gpe format
    -s  INT     The slopping length from the exon-intron joint to intron[Half of minAnchor]
    -l  STR     Sequencing library type,fr-unstranded,fr-firststrand,fr-secondstrand
    -r  DOU     The RNA Integrity Number[0]
    -h          Show this help information
EOF
    exit 0
}
[ $1 ] || usage
echo "--------------- Begin at $(date) ---------------"

outputDir=$PWD
RIN=0
ssRate=0
perlPath=/home/pengq/bin

while getopts "hi:o:a:b:f:g:r:s:l:" OPTION
do
    case $OPTION in
        h) usage;;
	a) fqInput1=$OPTARG;;
	b) fqInput2=$OPTARG;;
        i) InputBam=$OPTARG;;
        o) outputDir=$OPTARG;;
	f) refFa=$OPTARG;;
        g) gpeStructure=$OPTARG;;
        s) slop=$OPTARG;;
        l) libType=$OPTARG;;
        r) RIN=$OPTARG;;
        ?) usage;;
    esac
done
shift $((OPTIND - 1))
fqInput1=$(realpath $fqInput1)
fqInput2=$(realpath $fqInput2)
if  [ "X$libType" != "Xfr-unstranded" ] && [ "X$libType" != "Xfr-firststrand" ] && [ "X$libType" != "Xfr-secondstrand" ];then
  echo "Please specify correct libType by -l" >&2
  exit 1
fi
InputBam=$(readlink -f $InputBam)
outputDir=$(readlink -f $outputDir)
if [ ! -e $refFa ];then
  echo "The Reference sequences file $refFa doesn't exist" >&2
    exit 1
    fi
refFa=$(readlink -f $refFa)
if [ ! -e $gpeStructure ];then
  echo "The gene structure file file $gpeStructure doesn't exist" >&2
  exit 1
fi
gpeStructure=$(readlink -f $gpeStructure)

mkdir -p $outputDir
if [ $? -ne 0 ];then
  echo "Can't create directory $outputDir" >&2
  exit 1
fi
cd $outputDir
if [ -d log ];then
  [ -d log_old ] && rm -rf log_old
    mv log log_old
    fi
mkdir log
mkdir -p evaluation && cd evaluation
    
    echo 'Mapping Statistics (Background Running)'
    bamtools stats -in $InputBam -insert >bamStats.txt 
    
    if [ $(head -n 1 $gpeStructure | awk '{print NF}') -eq 15 ];then
      ln -s $gpeStructure myGpe.gpe
    else
      if [ $(head -n 1 $gpeStructure | awk '{print NF}') -eq 16 ];then
        cut -f2- $gpeStructure >myGpe.gpe
      else
        echo "Please offer the gene structure file in correct gpe format" >&2
        exit 1
      fi
    fi
    gpeStructure=$outputDir/evaluation/myGpe.gpe
 
    echo 'Fragmentation Evaluation (Background Running)'
    $perlPath/fragmentation.py --gpe $gpeStructure --bam $InputBam >fragmentation.tsv && $perlPath/fragmentation.R <fragmentation.tsv 

    #echo 'Mutation Rate Evaluation (Background Running)'
    #$perlPath/mutationRate_tophat.pl --sam $InputBam --ref $refFa >mutRate.tsv 2>$outputDir/log/mutRate.log && $perlPath/mutRate_tophat.R <mutRate.tsv 

    echo 'Coverage Evaluation (Background Running)'
    if [ "X$libType" == "Xfr-unstranded" ];then
        bedtools coverage -b $InputBam -a <($perlPath/gpeFeature.pl -e $gpeStructure | cut -f1-6) -hist -split >cov.tsv
    else
      if [ "X$libType" == "Xfr-secondstrand" ];then
        $perlPath/samFlag.pl -2 --rev $InputBam | samtools view -buS - 2>/dev/null | bedtools coverage -s -b /dev/stdin -a <($perlPath/gpeFeature.pl -e $gpeStructure | cut -f1-6) -hist -split >cov.tsv

     else
        $perlPath/samFlag.pl -1 --rev $Input | samtools view -buS - 2>/dev/null | bedtools coverage -s -b /dev/stdin -a <($perlPath/gpeFeature.pl -e $gpeStructure | cut -f1-6) -hist -split >cov.tsv
      fi
    fi &&grep ^all cov.tsv | cut -f 2,3,5 >cov_all.tsv && $perlPath/coverage.R <cov_all.tsv >cov_all_cum.tsv
       
    echo 'DNA Contamination Evaluation (Background Running)'
    if [ "X$libType" == "Xfr-unstranded" ];then
      $perlPath/samContam.pl -l $libType -g $gpeStructure -s $slop $InputBam >contam.tsv 2>$outputDir/log/contam.log
    else  
      $perlPath/samContam.pl -l $libType -g $gpeStructure -s $slop $InputBam >contam.tsv 2>$outputDir/log/contam.log
      echo 'Strand-specific Evaluation (Background Running)'
      $perlPath/gpeMerge.pl -n $gpeStructure | $perlPath/gpe2bed.pl >structure.bed
      bedtools subtract -a <(awk '$6=="+"' structure.bed) -b <(awk '$6=="-"' structure.bed) >noAntisenseRegion.bed
      bedtools subtract -a <(awk '$6=="-"' structure.bed) -b <(awk '$6=="+"' structure.bed) >>noAntisenseRegion.bed
      $perlPath/samSsEva.pl -b noAntisenseRegion.bed -l $libType $Input >ssEva.tsv 2>$outputDir/log/ssEva.log && ssRate=$(grep "^Correct Strand-specific Reads:" ssEva.tsv | grep -oP "[0-9.]+%" | sed 's/%$//')
    fi && echo "Run calRnaSeqRbScore.pl -r $RIN -s $ssRate -u $(grep "^Total reads:" bamStats.txt | grep -oP "\d+") --rpkmExon2Intron $(sed -n '5p;6q' contam.tsv | cut -f 1) --rpkmExon2Intergenic $(sed -n '6p;7q' contam.tsv | cut -f 1) -d $(awk '$1==10{print $2}' cov_all_cum.tsv) --fragmentation fragmentation.tsv -m $(awk '$5<=0.01' mutRate.tsv | wc -l),$(awk '$9<=0.01' mutRate.tsv | wc -l) -f $fqInput2 $fqInput1 >rbScore.tsv 2>$outputDir/log/rbScore.err"

   $perlPath/calRnaSeqRbScore.pl -r $RIN -s $ssRate -u $(grep "^Total reads:" bamStats.txt | grep -oP "\d+") --rpkmExon2Intron $(sed -n '5p;6q' contam.tsv | cut -f 1) --rpkmExon2Intergenic $(sed -n '6p;7q' contam.tsv | cut -f 1) -d $(awk '$1==10{print $2}' cov_all_cum.tsv) --fragmentation fragmentation.tsv -m $(awk '$5<=0.01' mutRate.tsv | wc -l),$(awk '$9<=0.01' mutRate.tsv | wc -l) -f $fqInput2 $fqInput1 >rbScore.tsv 2>$outputDir/log/rbScore.err

echo "--------------- End at $(date) ---------------"
