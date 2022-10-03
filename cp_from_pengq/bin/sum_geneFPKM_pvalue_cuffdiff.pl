#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($geneExp,$genesRPKM,$sampleNum,%hash);
GetOptions(
            'ip|pvalue=s'       => \$geneExp,
	    'ir|rpkm=s'         => \$genesRPKM,
            'n|sample=i'        => \$sampleNum,
	    'h|help'            => sub{usage()}
	  ) || usage();
open IP,"$geneExp" or die "Can't open file $geneExp:$!";
open IR,"$genesRPKM" or die "Can't open file $genesRPKM:$!";
while(<IP>){
    chomp;
    next if /^test_id/;
    my @split=split;
    $hash{$split[0]}="$split[11]"."_"."$split[12]"; #store pvalue for each gene
}
#print file header
print "#Gene";
my $tag=0;
while(<IR>){
    chomp;
    next if /^tracking_id/;
    $tag++;
    my @split=split;
    my $symbol="$split[1]"."_"."$split[2]";
    print "\t$symbol";
    last if $tag==$sampleNum;
}
print "\tp_value\tq_value\n";
#print FPKM and pvalue
seek IR,0,0;
$tag=0;
while(<IR>){
    chomp;
    next if /^tracking_id/;
    $tag++;
    my @split=split;
    my $fpkm=$split[6];
    if($tag==1){
        print "$split[0]\t$fpkm";
    }elsif($tag==$sampleNum){
        my @value=split /_/,$hash{$split[0]};
        print "\t$fpkm\t$value[0]\t$value[1]\n";
        $tag=0;
    }else{
        print "\t$fpkm";
    }
}
sub usage{
print <<HELP;
Description: Output each gene's FPKM in each sample and the pvalue from cuffdiff results.
Usage: perl $0 -ip gene_exp.diff -ir genes.read_group_tracking >result
Output file:gene  condition1_rep1 condition1_rep2..  condition2_rep1  condition2_rep2.. pvalue  FDR
Options:
    -ip  FILE          The file contained genes and pvalue.
    -ir  FILE          The file contained genes' FPKM in each sample.
    -n   INT           The total sample number used to run cuffdiff.
    -h --help          Print this help information
HELP
    exit(-1);
}