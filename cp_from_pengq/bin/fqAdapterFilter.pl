#!/usr/bin/perl -w
use strict;
use warnings;
use 5.010;
use Getopt::Long;
use File::Basename;
my $scriptName=(fileparse($0))[0];
#my ($source,$target)=('fastq-illumina','fastq-illumina');
my ($adaptor,$adaptor2)=('ACACTCTTTCCCTACACGACGCTCTTCCGATCT','CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT');
my ($input,$output,$input2,$output2,$report,$help);
my $opt=GetOptions(
                        'a=s'=>\$adaptor,
                        'i=s' => \$input,
                        'o=s' => \$output,
                        'a2=s'=>\$adaptor2,
                        'i2=s' => \$input2,
                        'o2=s' => \$output2,
                        'r|report=s' =>\$report,
                        'h|help'=> \$help
                        );
usage() if(defined $help);
if(defined $input){
   #Add processing gz file by Yumei Li
   if($input=~/.gz$/){
      open STDIN,"gunzip -c $input|" or die "Can't open file $input:$!";   
   }else{
      open STDIN,"<$input" or die "$input:$!";
   }
}
unless( (defined $input2 && defined $output2) || (!defined $input2 && !defined $output2) ){
    say "Please specify all the argument values of --i2 and --o2"
}
if(defined $output){
   open OUT,">$output" or die "$output:$!";
   select OUT;
}
my ($id,$seq,$plus,$quality,@adaptors);
my ($counterOfDiscard,$counterOfAll)=(0,0);
if(-f $adaptor){
   open AD,"$adaptor" or die "Can't open file $adaptor:$!";
   while(<AD>){
      push @adaptors,$_;
   }
}else{
   @adaptors=split /,/,$adaptor;
}
#SE reads
if (!defined $input2){
   while($id=<STDIN>){
      $seq=<STDIN>;
      $plus=<STDIN>;
      $quality=<STDIN>;
      my $tag=0;
      for(my $i=0;$i<=$#adaptors;$i++){
         if($seq=~/^$adaptors[$i]/){
            $tag++;
            last;
         }
      }
      if($tag==0){
      	print $id,$seq,$plus,$quality;
      }else{
      	$counterOfDiscard++;
      }
      $counterOfAll++;
   }
}else{#PE reads 
    my ($id2,$seq2,$plus2,$quality2,@adaptor2s);
    if($input2=~/.gz$/){
      open IN2,"gunzip -c $input2|" or die "Can't open file:$!";  
    }else{
      open IN2,"<$input2" or die "$input2:$!";
    }
    open OUT2,">$output2" or die "$output2:$!";
    if(-f $adaptor2){
      open AD2,"$adaptor2" or die "Can't open file $adaptor2:$!";
      while(<AD2>){
         push @adaptor2s,$_;
      }
    }else{
      @adaptor2s=split /,/,$adaptor2;
    }
    while($id=<STDIN>){
       $seq=<STDIN>;
       $plus=<STDIN>;
       $quality=<STDIN>;
       $id2=<IN2>;
       $seq2=<IN2>;
       $plus2=<IN2>;
       $quality2=<IN2>;
       my($tag,$tag2)=(0,0);
       for(my $i=0;$i<=$#adaptors;$i++){
         if($seq=~/^$adaptors[$i]/){
            $tag++;
            last;
         }
       }
       for(my $j=0;$j<=$#adaptor2s;$j++){
         if($seq2=~/^$adaptor2s[$j]/){
            $tag2++;
            last;
         }
       }
       if($tag==0 && $tag2==0){
          print $id,$seq,$plus,$quality;
          print OUT2 $id2,$seq2,$plus2,$quality2;
       }
       else{
            $counterOfDiscard++;
       }
       $counterOfAll++;
    }

}
#revised by Yumei Li
#if(defined $report){
#    select STDOUT;
#    say "No.Reads of All=\t$counterOfAll";
#    say "No.Reads discarded by $scriptName=\t$counterOfDiscard";
#}
open REPORT,">$report" or die $!;
say REPORT "No.Reads of All=\t$counterOfAll";
say REPORT "No.Reads discarded by $scriptName=\t$counterOfDiscard";
sub usage{
print <<HELP;
Usage:perl $scriptName [-a ACACTCTTTCCCTACACGACGCTCTTCCGATCT] [-i input.fq(.gz)] [-o output.fq(.gz)] [--i2 input2.fq(.gz) --o2 output2.fq(.gz)]
Revised by Yumei Li: Add the processing of mutiple adaptors.(2017-5-15)
Revised by Yumei Li: Add processing .gz input files;
Revised by Yumei Li: Add processing mutiple adaptors stored in a file (2018-12-26);
    -a         adapter sequence,default is ACACTCTTTCCCTACACGACGCTCTTCCGATCT(Illumina Paired End Adapter 1)(Mutiple adapters can be provided separated by comma or stored in a file with each adapter in a line)
    -i         input file in fq or fq.gz format,default is STDIN
    -o         output file,default is STDOUT
    --a2       adaptor2 sequence,default is CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT(Illumina Paired End Adapter 2)(Mutiple adapters can be provided separated by comma or stored in a file with each adapter in a line)
    --i2       (optional)input file 2 if the reads are pair in fq or fq.gz format
    --o2       (optional)output file 2 if the reads are pair
    -r --report Print result report to this file.
    -h --help   This help information screen
HELP
    exit(-1);
}
