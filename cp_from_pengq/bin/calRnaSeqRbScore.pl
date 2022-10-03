#!/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use Getopt::Long;
use File::Basename;

my ($fastqc, $lowMismatchBase, $uniqReadCount, $highDepthFraction, $fragmentationFile, $rpkmExon2Intron, $rpkmExon2Intergenic);
my ($RIN, $minQual, $intervalN, $ssRate) = (0, 20, 25, 0);
GetOptions(
            'r|RIN=i'               => \$RIN,
            'f|fastqc=s'            => \$fastqc,
            'q|minQual=i'           => \$minQual,
            's|ssRate=s'            => \$ssRate,
            'u|uniqRead=i'          => \$uniqReadCount,
            'm|lowMismatchBase=s'   => \$lowMismatchBase,
            'rpkmExon2Intron=s'     => \$rpkmExon2Intron,
            'rpkmExon2Intergenic=s' => \$rpkmExon2Intergenic,
            'd|highDepthFraction=s' => \$highDepthFraction,
            'fragmentation=s'       => \$fragmentationFile,
            'i|intervalN=i'         => \$intervalN,
            'h|help'                => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";
open FRAG, "$fragmentationFile" or die "Can't open $fragmentationFile: $!";

my ($totalReadCount, $read1Len);
my $read1HighQualBase = 0;
while(<IN>){
    chomp;
    $totalReadCount = $1 if /Total Sequences\t(\d+)/;
    $read1Len = $1 if /Sequence length\t(\d+)/;
    if(/^>>Per base sequence quality/){
        <IN>; # skip title line
        while(my $line = <IN>){
            chomp $line;
            last if $line =~ /^>>END_MODULE/;
            my ($pos, $mean) = split "\t", $line;
            next if $mean < $minQual;
            if($pos =~ /-/){
                my ($start, $end) = split "-", $pos;
                $read1HighQualBase += $end - $start + 1;
            }else{
                $read1HighQualBase++;
            }     
        }
    }
}

my ($read2Len);
my $read2HighQualBase = 0;
if(defined $fastqc){
    open F, "$fastqc" or die "Can't open $fastqc: $!";
    while(<F>){
        chomp;
        $totalReadCount += $1 if /Total Sequences\t(\d+)/;
        $read2Len = $1 if /Sequence length\t(\d+)/;
        if(/^>>Per base sequence quality/){
            <F>;
            while(my $line = <F>){
                chomp $line;
                last if $line =~ /^>>END_MODULE/;
                my ($pos, $mean) = split "\t", $line;
                next if $mean < $minQual;
                if($pos =~ /-/){
                    my ($start, $end) = split "-", $pos;
                    $read2HighQualBase += $end - $start + 1;
                }else{
                    $read2HighQualBase++;
                } 
            }
        }
    }
    
}

my $sFragmentation = 0;
while(<FRAG>){
    chomp;
    my $fraction = (split "\t")[1];
    $sFragmentation += abs($fraction - 1/$intervalN);
}
$sFragmentation = 1 - $sFragmentation/2/(1-1/$intervalN);

my $sRIN = ($RIN >= 10 ? 1 : $RIN/10);

my ($sReadLen, $sHighQual, $sMuta);
my @lowMismatchBases = split ',', $lowMismatchBase;
if(defined $fastqc){
    $sReadLen = ( ($read1Len + $read2Len) >= 100 ? 1 : ($read1Len + $read2Len)/ 100);
    $sHighQual = ($read1HighQualBase / $read1Len + $read2HighQualBase / $read2Len)/2;
    $sMuta = ($lowMismatchBases[0] / $read1Len + $lowMismatchBases[1] / $read2Len) / 2;
}else{
    $sReadLen = ($read1Len >= 100 ? 1 : $read1Len / 100);
    $sHighQual = $read1HighQualBase / $read1Len;
    $sMuta = $lowMismatchBases[0] / $read1Len;
}

my $sMapEff = $uniqReadCount / $totalReadCount;

my $sContam = ($rpkmExon2Intron > 80 ? 1 : $rpkmExon2Intron / 80) + ($rpkmExon2Intergenic > 100 ? 1 : $rpkmExon2Intergenic / 100);

my $sSsRate = $ssRate / 100;

my $total = $sRIN + $sReadLen + $sHighQual + $sSsRate + $sMapEff + $sMuta + $sContam + $highDepthFraction + $sFragmentation;

say "RIN\t$sRIN";
say "ReadLength\t$sReadLen";
say "BaseQuality\t$sHighQual";
say "StrandSpecificity\t$sSsRate";
say "MappingEfficiency\t$sMapEff";
say "MutationFraction\t$sMuta";
say "Contamination\t$sContam";
say "TranscriptomeCoverage\t$highDepthFraction";
say "Fragmentation\t$sFragmentation";
say "Total\t$total";

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName fastqc_data.txt >score.tsv
    If fastqc_data.txt isn't specified, input from STDIN
Options:
    -r --RIN                    INT         The RNA Integrity Number[0]
    -f --fastqc                 FILE        (Optional) The fastqc_data.txt file for read2
    -q --minQual                INT         The minimal quality to be taken as high-quality base[20]
    -s --ssRate                 DOU         The ratio (0-100) of reads that are consistent with strand-specific principle[0]
    -u --uniqRead               INT         The uniq read number
    -m --lowMismatchBase        INT1[,INT2] The number(s) of bases with low mismatch rate
       --rpkmExon2Intron        DOU         The ration of exon RPKM to intron RPKM
       --rpkmExon2Intergenic    DOU         The ration of exon RPKM to intergenic RPKM
    -d --highDepthFraction      DOU         The fraction of high depth bases in whole transcriptome
       --fragmentation          FILE        The fragmentation evaluation result file
    -i --intervalN              INT         The interval number used in fragmentation evaluation[25]
    -h --help                               Print this help information
HELP
    exit(-1);
}