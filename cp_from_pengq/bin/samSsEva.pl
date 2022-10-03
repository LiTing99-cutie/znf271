#!/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;
use pm::samParser;

my $bedFile;
my ($slop, $libType) = (0, 'fr-firststrand');
GetOptions(
            'b|bed=s'           => \$bedFile,
            's|slop=i'	        => \$slop,
            'l|libraryType=s'   => \$libType,
            'h|help'	        => sub{&usage}
          ) || usage();

if(defined $ARGV[0]){
    if (-B $ARGV[0]){
        open IN, "samtools view $ARGV[0]|" or die "Can't open $ARGV[0]: $!";
    }else{
        open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";
    }
}else{
    open IN, "-";
}

open BED, $bedFile or die "Can't open $bedFile: $!";

my %structure;
while(<BED>){
    chomp;
    my ($chr, $start, $end, $strand) = (split "\t")[0..2,5];    
    push @{$structure{$chr}}, [$start, $end, $strand];
}
for my $chr(keys %structure){
    my @sorted = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]}@{$structure{$chr}};
    $structure{$chr} = \@sorted;
}

my $line;
while($line = <IN>){ #skip header line
    last if $line !~ /^@/;
}

my ($totalNo, $readNo, $properReadNo, $strandErrorNo) = (0, 0, 0, 0);

while(defined $line){
    chomp $line;
    $totalNo++;
    my @fields = split "\t", $line;
    say STDERR common::getFormatedTime() . " $totalNo reads have been processed" if ($totalNo%100000 == 0);
    my ($chr, $CIGAR) = @fields[2, 5];
    if($CIGAR eq '*'){
	$line = <IN>;
	next;
    }
    my $flag = $fields[1];
    my $tags = join "\t", @fields[11..$#fields];
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

    my ($readS, $readE);
    $readE = $readS = $fields[3] - 1;
    $readE += $_ for($CIGAR=~/(\d+)[MDN=X]/g);
    if(!exists $structure{$chr}){ #this ref isn't in structure
	$line=<IN>;
	next;
    }
    for my $region(@{$structure{$chr}}){
	my ($start, $end, $strand) = @$region;
	last if $start - $slop >= $readE;
	if($start - $slop <= $readS && $readE <= $end + $slop){
	    if($strand eq $codingStrand){
		$properReadNo++;
	    }else{
		$strandErrorNo++;
	    }
	    $readNo++;
	    last;
	}
    }
    $line=<IN>;
}

say "Total Reads falling in the providing regions :\t$readNo";
say "Correct Strand-specific Reads:\t$properReadNo (". $properReadNo/$readNo*100 . "%)";
say "Incorrect Strand-specific Reads:\t$strandErrorNo (" . $strandErrorNo/$readNo*100 . "%)";

sub usage{
    my $scriptName = basename $0;    
print STDERR <<HELP;
Usage: perl $scriptName -b gene_structure.bed [-s 0] SAMPLE.SAM/BAM >result.log 2>error.log
    if SAMPLE.SAM/BAM isn't specified, input is from STDIN
    output to STDOUT
    -b --bed	        FILE	A bed6 or bed6+ file with non-overlapped gene body regions between two strands
                                A merged bed is suggested for speeding up assignment of read to gene body regions
    -s --slop	        INT     Specify the slop length[0]
    -l --libraryType    STR	The library type, it can be
				    fr-firststrand: for dUTP, NSR, NNSR (default)
				    fr-secondstrand: for Ligation, Standard SOLiD and Illumina Directional Protocol
    -h --help	                Print this help information    
HELP
    exit(-1);
}