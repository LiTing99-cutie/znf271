#!/bin/usr/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inFile,$chromSizeFile,$bedFile,$bigWig);
GetOptions(
            'i|bam=s'          => \$inFile,
            'c|chrom=s'        => \$chromSizeFile,
	    'b|bed=s'          => \$bedFile,
	    'w|bw=s'           => \$bigWig,
            'h|help'           => sub{usage()}
	  ) || usage();
my ($IN,$step,$chromSize,$totalCov,$meanCov);
open $IN,"$inFile" or die "Can't open file $inFile:$!";
#Calculate total coverage
while(<$IN>){
    chomp;
    next if /^#/;
    next if /^tract/;
    if(/^fixedStep/){
        $_=~/.* step=(\d+)/;
        $step=$1;
        next;
    }else{
        $totalCov+=$_*$step;
    }
}

if(defined $chromSizeFile){
    open CH,"$chromSizeFile" or die "Can't open file $chromSizeFile:$!";
    while(<CH>){
	chomp;
	my($chr,$size)=split /\t/;
	$chromSize+=$size;
    }
    $meanCov=$totalCov/$chromSize;
}
#Need bwtool 
if(defined $bedFile){
    `cut -f1-3 $bedFile|bwtool summary stdin $bigWig region_sum.txt -total -header -with-sum`;
    open SUM,"region_sum.txt" or die "Can't open file region_sum.txt:$!";
    my $line1=<SUM>;
    chomp(my $line2=<SUM>);
    my @split=split /\t/,$line2;
    $meanCov=$split[-1]/$split[3];
}

seek $IN,0,0;
while(<$IN>){
    chomp;
    next if /^#/;
    next if /^track/;
    if(/^fixedStep/){
        say $_;
    }else{
	my $res=sprintf "%.2f",$_/$meanCov;
        say $res;
    }
}
          
sub usage{
print <<HELP;
Usage:perl $0 -i *.wig >out.wig
Author:Yumei Li,2016-4-14
       Revised by Yumei Li,2016-6-8: Adding normalizing by the average base pair coverage of specific input region.Requiring bwtool in your PATH.
Description:This script will normalize the input wig file by dividing the average base pair coverage of the whole genome or specific input region.
Options:
    -i|--wig    The input wig format file(fixed step size)
    -c|--chrom  The genomic chromsome size file.[Two columns:chr    size]
    -b|--bed    The specific region used to normalize data in bed format.(Incompatible with -c,force -w)
    -w|--bw     The same input file in bigWig format.
    -h|--help   Print this help information.
HELP
    exit(-1);
}
