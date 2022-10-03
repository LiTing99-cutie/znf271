#!/bin/env perl

use warnings;
use strict;
#use lib "/rd1/user/chenjy/perl_lib/lib/perl5/";
use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case bundling);
use List::Util qw (sum min max shuffle);
use vars qw($sam $ref  $uniq $length $title $version);
GetOptions(
	"sam|i:s"	=>	\$sam,
	"ref|r:s"	=>	\$ref,
	"uniq|u"	=>	\$uniq,
	"title|t"	=>	\$title,
	"version|v"	=>	sub{&version;exit(-1);},
	"help|h"	=>	sub{&usage;exit(-1);}
);
my $fh;
unless($ref){&usage;exit(0);}
unless($sam){$fh = \*STDIN;}
else{unless (-B $sam){open $fh,$sam;} else{open $fh,"samtools view $sam|";}}


my (%mis,%insertion,%deletion);
my $line = 0;
my $line1 = 0;
my $line2 = 0;
while(<$fh>){
	next if /^@/;
# only uniquely-mapped reads are taken into consideration	
	if(defined $uniq){
		next unless ($_ =~ /XT:A:U/ || $_ =~ /NH:i:1\D/ || $_ =~ /X0:i:1\D/);
	}
# chomp after pattern match
	chomp;
	my @field = split("\t",$_);
	my ($flag,$chr,$start,$cigar,$seq) = @field[1,2,3,5,9];
	next if $cigar =~ /[SHP=X]/;
	
# only mapped reads are taken into consideration
	next if ($flag&0x4 || $cigar eq "*");
	
	$line ++;
	if ($flag & 0x80){
		$line2 ++;
	}else{
		$line1 ++;
	}
	showLog("$line lines have been processed...") unless ($line%10000);
	next if (/NM:i:0/);
	
	if (!defined $length){$length = length($field[9]);}
# record deletions
	if ($cigar =~ /D/){
		my $copy = $cigar;
		my $num;
		if ($flag&0x10){
			my @new;
			foreach(;$copy=~s/(\d+[MIDN])//;){
				unshift @new,$1;
			}
			$copy = join('',@new);
		}
		foreach(;$copy =~ s/(\d+)([MID])//;){
			my ($tmp1,$tmp2) = ($1,$2);
			if ($tmp2 ne 'D'){$num += $tmp1;}
			if ($tmp2 eq 'D'){
				unless ($flag&0x80){
					$deletion{'1st'}->{$num} ++;
				}else{
					$deletion{'2nd'}->{$num} ++;
				}
			}
		}
	}
# record insertion & mismatch
	my ($MIS,$INSERTION) = &mismatch_position(&get_genome_seq_base_on_start_and_cigar($ref,$chr,$start,$cigar),$seq);
	foreach(@{$MIS}){
		unless ($flag&0x80){ 
			unless($flag&0x10){$mis{'1st'}->{$_}++;}
			else{$mis{'1st'}->{length($seq)-$_+1}++;}
		}else{
			unless($flag&0x10){$mis{'2nd'}->{$_}++;}
			else{$mis{'2nd'}->{length($seq)-$_+1}++;}
		}
	}
	foreach(@{$INSERTION}){
		unless ($flag&0x80){
			unless($flag&0x10){$insertion{'1st'}->{$_}++;}
			else{$insertion{'1st'}->{length($seq)-$_+1}++;}
		}else{
			unless($flag&0x10){$insertion{'2nd'}->{$_}++;}
			else{$insertion{'2nd'}->{length($seq)-$_+1}++;}
		}
	}
}

if (defined $title){
	if ($line2){
		print "#total reads under check:$line($line1,$line2)\n";
		print "#pos\tmis1st\tins1st\tdel1st\tper_read1\tmis2nd\tins2nd\tdel2nd\tper_read2\tper_total\n";
	}else{
		print "#total reads under check:$line($line1)\n";
		print "#pos\tmis1st\tins1st\tdel1st\tper\n";
	}
	
}

foreach(1..$length){
	print "$_\t";

	unless (exists $mis{'1st'}->{$_}){$mis{'1st'}->{$_} = 0;}
	unless (exists $insertion{'1st'}->{$_}){$insertion{'1st'}->{$_} = 0;}
	unless (exists $deletion{'1st'}->{$_}){$deletion{'1st'}->{$_} = 0;}
	if($line2 != 0){
		unless (exists $mis{'2nd'}->{$_}){$mis{'2nd'}->{$_} = 0;}
		unless (exists $insertion{'2nd'}->{$_}){$insertion{'2nd'}->{$_} = 0;}
		unless (exists $deletion{'2nd'}->{$_}){$deletion{'2nd'}->{$_} = 0;}
	}	
	print $mis{'1st'}->{$_},"\t",$insertion{'1st'}->{$_},"\t",$deletion{'1st'}->{$_},"\t";
	print sum($mis{'1st'}->{$_},$insertion{'1st'}->{$_},$deletion{'1st'}->{$_})/$line1;
	if ($line2 != 0){
		print "\t";
		print $mis{'2nd'}->{$_},"\t",$insertion{'2nd'}->{$_},"\t",$deletion{'2nd'}->{$_},"\t";
		print sum($mis{'2nd'}->{$_},$insertion{'2nd'}->{$_},$deletion{'2nd'}->{$_})/$line2;
		print "\t";
		print sum($mis{'1st'}->{$_},$insertion{'1st'}->{$_},$deletion{'1st'}->{$_},$mis{'2nd'}->{$_},$insertion{'2nd'}->{$_},$deletion{'2nd'}->{$_})/$line;
	}
	print "\n";
}


# subrountines

# fetch genome sequence based on start and cigar
# insertion and clipping were mark with '-'
sub get_genome_seq_base_on_start_and_cigar{
        my ($ref,$chr,$start,$cigar) = @_;
        my $db = Bio::DB::Fasta->new($ref);
        my $seq;
	my $geno_pos = &genome_location_base_on_cigar($chr,$start,$cigar);
        my @starts = split(',',$geno_pos->[0]);        
        my @ends = split(',',$geno_pos->[1]);
        for(my $i = 0; $i <= $#starts ;$i++){
                if ($starts[$i] =~ /(\d+)_bp/){
                        $seq .= '-'x$1;
                }else{
                        $seq .= $db->seq($chr,$starts[$i],$ends[$i]);
                }
        }
        return (uc($seq));
}

# genome location based on CIGAR
sub genome_location_base_on_cigar{
        my ($chr,$start,$cigar) = @_;
        my (@starts,@ends);
        my $start_before = $start;

        for(;$cigar=~s/(\d+)([DNIM])//;){
                my ($num,$type) = ($1,$2);
                if ($type eq 'M'){
                        push @starts,$start_before;
                        push @ends,$start_before + $num - 1;
                        $start_before = $start_before + $num - 1;
                }
                if ($type eq 'D' || $type eq 'N'){
                        $start_before += ($num + 1);
                }
                if ($type eq 'I'){
                        push @starts,"${num}_bp_insertion";
                        push @ends,"${num}_bp_insertion";
                        $start_before ++;
                }
        }
        my @return = (join(',',@starts),join(',',@ends));
        return(\@return);
}

# record the postion of insertion and mismatch
sub mismatch_position{
        my ($seq1,$seq2) = @_;
        my @mis;
        my @insertion;
        for(my $i=0;$i<length($seq1);$i++){
                my $tmp = substr($seq2,$i,1);
                if(substr($seq1,$i,1) eq '-'){
                        push @insertion,$i+1;
                }elsif(substr($seq1,$i,1) !~ /$tmp/i){
                        push @mis,$i+1;
                }
        }
        return(\@mis,\@insertion);
}

# trace the process
sub showLog{
        my ($info) = @_;
        my @times = localtime; # sec, min, hour, day, month, year
        print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}

# usage
sub usage{
print STDERR <<HELP 
Usage:	perl $0 --sam [input] --ref [reference.fa] --uniq|-u  --title|-t  --help|-h 
	--sam|-i:s		sam or bam file or STDIN
	--ref|-r:s		reference genome, fasta format 
	--uniq|-u		only unique reads defined by XT:A:U, NH:i:1 or X0:i:1 tag are considered
	--title|-t		get result with title
	--version|-v		release note
	--help|-h		print this help message
HELP
}

sub version{
print STDERR<<VERSION

Author:
	Jerry Chen
Last Modified:
	07/02/2014
Function:
	mismatch fraction for strand-specific RNA-Seq data (1st first cDNA, 2nd second cDNA strand). For non-strand-specific & single-end RNA-seq, however, it still works. MD tag-independent.
Limitation:
	1). only work for SAM/BAM with CIGAR, either one of XT,NH,X0 tags.
	2). all reads should have the same length
Release Note
	\@version 1.0
		beta version
	\@version 1.1
		add --version|-v option
		add percentage in output
	\@version 2.0
		single-end RNA-Seq is surpported
	\@version 2.1
		chomp after pattern match

VERSION
}
