#!/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::gpeParser;
use pm::samParser;

my ($gpeFile, $bin, $genomeLen, $chrSizeFile);
my ($slop, $libType) = (0, 'fr-unstranded');
GetOptions(
            'g|gpe=s'       => \$gpeFile,
            'b|bin'         => \$bin,
            'l|libType=s'   => \$libType,
            'gLen=s'        => \$genomeLen,
            'c|chrSize=s'   => \$chrSizeFile,
            's|slop=i'      => \$slop,
            'h|help'        => sub{&usage()}
        ) || usage();

if(defined $ARGV[0]){
    if (-B $ARGV[0]){
        open IN, "samtools view -h $ARGV[0]|" or die "Can't open $ARGV[0]: $!";
    }else{
        open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";
    }
}else{
    open IN, "-";
}

open ANO, "$gpeFile" or die "Can't open $gpeFile: $!";

chomp(my $line = <IN>);
if( defined $genomeLen ){
    if($genomeLen eq 'hg19'){
        $genomeLen = 3_137_161_264;
    }elsif($genomeLen eq 'hg18'){
        $genomeLen = 3_107_677_273;
    }elsif($genomeLen eq 'rheMac2'){
        $genomeLen = 2_864_122_635;
    }else{
        die "Please specify the correct --gLen\n" unless $_ =~ /^\d+$/
    }
}elsif( defined $chrSizeFile ){
    $genomeLen = 0;
    open CHR,"$chrSizeFile" or die "Can't open $chrSizeFile: $!";
    while(<CHR>){
        chomp;
        $genomeLen += (split "\t")[1];
    }
}else{
    $genomeLen = 0;
    while(defined $line){
        last if $line !~ /^\@/;
        $genomeLen += $1 if $line =~ /^\@SQ\tSN:.+\tLN:(\d+)/;
        $line = <IN>;
    }
}
die "Total length of genome is 0\n" if $genomeLen == 0;

my %gpeHash;
my $exomeLen = 0;
my $intromeLen = 0;
while(<ANO>){
    chomp;
    my @fields = split "\t";
    shift @fields if defined $bin;
    my ($name, $chr, $strand, $start, $end) = @fields[0..4];
    my @blockStarts = split ",", $fields[8];
    my @blockEnds = split ",", $fields[9];
    my $exonsLen = gpeParser::getExonsLength(\@blockStarts, \@blockEnds);
    $exomeLen += $exonsLen;
    $intromeLen += $end - $start - $exonsLen;
    push @{$gpeHash{$chr}{$strand}}, [ $start, $end, \@blockStarts, \@blockEnds ];
}

while(my ($chr, $chrV) = each %gpeHash){
    for my $strand (keys %$chrV){
        my @sortedStrandV = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]}@{$chrV->{$strand}};
        $chrV->{$strand} = \@sortedStrandV;
    }
}

my ($exonicReadN, $intronicReadN, $intergenicReadN, $properMapReadN, $totalReadN) = (0, 0, 0, 0, 0);
READ:
while(defined $line){
    $totalReadN++;
    chomp $line;
    my @fields = split "\t", $line;
    my ($name, $flag, $ref, $start, $cigar) = @fields[0, 1, 2, 3, 5];
    my $tags = join "\t", @fields[11..$#fields];
    say STDERR "#$totalReadN reads ($name) have been processed" if $totalReadN % 100000 == 0;
    if($cigar eq '*'){
        $line = <IN>;
        next READ;
    }
    my $codingStrand = samParser::determineCodingStrand($libType, $flag);
    my $XS = $1 if $tags =~ /XS:A:([+-])/;
    if(!defined $codingStrand){
        die "Please specify correct library type by --libType\n";
    }elsif($codingStrand eq '.'){# fr-unstranded
        $codingStrand = $1 if defined $XS;
    }elsif($codingStrand eq ''){
        say STDERR join "\t", ("mappedAgainstLibType", $line);
        $line = <IN>;
        next;
    }elsif(defined $XS && $codingStrand ne $XS){
        say STDERR join "\t", ("strandContradict", $line);
        $line = <IN>;
        next;
    }
    $properMapReadN++;
    $start--;
    my $end = $start;
    $end += $_ for( $cigar =~ /(\d+)[MDN=X]/g );
    my $isIntergenicRead = 1;
    my @strands = $codingStrand eq '.' ? ('+', '-') : ($codingStrand);
    for my $strand(@strands){
        for my $trans (@{$gpeHash{$ref}{$strand}}){
            my ($transStart, $transEnd, $blockStarts, $blockEnds) = @$trans;
            last if $transStart - $slop >= $end;
            if( $start >= $transStart - $slop && $end <= $transEnd + $slop){#read embeded in transcript
                $isIntergenicRead = 0;
                my ($startI, $endI,);
                my ($isStartInExon, $isEndInExon) = (0, 0);
                for( $startI = 0; $startI < @$blockStarts; $startI++){
                    if ($start >= $blockStarts->[$startI] - $slop && $start < $blockEnds->[$startI] + $slop){
                        $isStartInExon = 1;
                        last;
                    }
                    last if $startI < @$blockStarts -1 && $start < $blockStarts->[$startI +1] - $slop && $start >= $blockEnds->[$startI] + $slop;
                }
                for($endI = $startI; $endI < @$blockStarts; $endI++){
                    if( $end >= $blockStarts->[$endI] - $slop && $end < $blockEnds->[$endI] + $slop ){
                        $isEndInExon = 1;
                        last;
                    }
                    last if $endI <@$blockStarts-1 && $end < $blockStarts->[$endI + 1] - $slop && $end >= $blockEnds->[$endI] + $slop;
                }
                if( $isStartInExon == 1 && $isEndInExon == 1 ){
                    $line = <IN>;
                    $exonicReadN++;
                    next READ;
                }
            }
        }
    }
    if($isIntergenicRead == 1){
        $intergenicReadN++;
    }else{
        $intronicReadN++;
    }
    $line = <IN>;
}

my $intergenomeLen = $genomeLen*2 - $exomeLen - $intromeLen;
my $exomeRPKM = $exonicReadN / ($exomeLen/1e3) / ($properMapReadN/1e6);
my $intromeRPKM = $intronicReadN / ($intromeLen/1e3) / ($properMapReadN/1e6);
my $intergenomeRPKM = $intergenicReadN / ($intergenomeLen/1e3) / ($properMapReadN/1e6);

say "#Total reads=\t$totalReadN";
say "#Properly mapped reads=\t$properMapReadN";

say join "\t", ("#RPKM Ratio of ExonmeRPKM to this", "RPKM", "Length", "Read NO");
say join "\t", (1, $exomeRPKM, $exomeLen, $exonicReadN);
say join "\t", ($exomeRPKM/$intromeRPKM, $intromeRPKM, $intromeLen, $intronicReadN);
say join "\t", ($exomeRPKM/$intergenomeRPKM, $intergenomeRPKM, $intergenomeLen, $intergenicReadN);

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName -g gene_structure.gpe -s 4 INPUT.SAM/BAM >result.tsv 2>report.log
    If INPUT.SAM/BAM isn't specified, input from STDIN
    Output result to STDOUT
    This tool can be used to calculate RPKM of exome, introme and intergeome and their ratios
    
    -g --gpe        FILE    A gpe file of transcript structure
    -b --bin                With bin column in --gpe
    -l --libType    STR     The library type, it can be
                                fr-unstranded: for Standard Illumina (default)
                                fr-firststrand: for dUTP, NSR, NNSR
                                fr-secondstrand: for Ligation, Standard SOLiD and Illumina Directional Protocol
        --gLen      INT/STR (Optional)The haploid genome length of organsim assembly. It can be:
                            'hg19': length = 3,137,161,264
                            'hg18': length = 3,107,677,273
                            'rheMac2': length = 2,864,122,635
                            INT: specify a number                            
    -c --chrSize    FILE    (Optional)If --gLen isn't specified, a tab-separated chromosome size file
                            with 'chr name' and 'size' columns can be specified. If both --gLen and
                            --chrSizes aren't specified, get the chromosome size from sam SQ header and
                            summarize genome length
    -s --slop               Specify the slopping length from the exon-intron joint to intron[0]
    -h --help               Print this help information
HELP
    exit(-1);
}
