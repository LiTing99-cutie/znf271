#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use List::Util qw(min max);
my ($bed,$peakFile,$heatmap);
my ($upstream,$downstream,$point)=(1000,1000,"c");
my $opt=GetOptions(
    'i|peak=s'      => \$peakFile,
    'b|bed=s'       => \$bed,
    'p|point=s'     => \$point,
    'u|up=i'        => \$upstream,
    'd|down=i'      => \$downstream,
	'm|hmp=s'       => \$heatmap,
    'h|help'	=> sub{&usage;exit(-1);}
                  );
my $IN;
chomp(my $bedColumn=`head -n1 $bed|awk '{print NF}'`);
if($bedColumn>=6){
	if($point eq "c"){
		open $IN,"awk -v OFS='' -v ORS='' -v up=$upstream -v dn=$downstream '{print \$1\"\\t\";if(\$6==\"+\"){print \$2+int((\$3-\$2)/2)-up\"\\t\"\$2+int((\$3-\$2)/2)+1+dn;}else{print \$3-int((\$3-\$2)/2)-1-dn\"\\t\"\$3-int((\$3-\$2)/2)+up};for(i=4;i<=NF;i++){print \"\\t\"\$i}print \"\\n\"}' $bed|awk '\$2>0'|bedtools intersect -wao -a stdin -b $peakFile|" or die "Can't open file:$!";
	}elsif($point eq "s"){
		open $IN,"awk -v OFS='' -v ORS='' -v up=$upstream -v dn=$downstream '{print \$1\"\\t\";if(\$6==\"+\"){print \$2-up\"\\t\"\$2+1+dn;}else{print \$3-1-dn\"\\t\"\$3+up;};for(i=4;i<=NF;i++){print \"\\t\"\$i}print \"\\n\"}' $bed|awk '\$2>0'|bedtools intersect -wao -a stdin -b $peakFile|" or die "Can't open file:$!";
	}elsif($point eq "e"){
		open $IN,"awk -v OFS='' -v ORS='' -v up=$upstream -v dn=$downstream '{print \$1\"\\t\";if(\$6==\"+\"){print \$3-1-up\"\\t\"\$3+dn;}else{print \$2-dn\"\\t\"\$2+1+up;};for(i=4;i<=NF;i++){print \"\\t\"\$i}print \"\\n\"}' $bed|awk '\$2>0'|bedtools intersect -wao -a stdin -b $peakFile|" or die "Can't open file:$!";
	}else{
		say "Please set correct parameter for -p";
	}
}else{
	if($point eq "c"){
		open $IN,"awk -v OFS='' -v ORS='' -v up=$upstream -v dn=$downstream '{print \$1\"\\t\";print \$2+int((\$3-\$2)/2)-up\"\\t\"\$2+int((\$3-\$2)/2)+1+dn;for(i=4;i<=NF;i++){print \"\\t\"\$i}print \"\\n\"}' $bed|awk '\$2>0'|bedtools intersect -wao -a stdin -b $peakFile|" or die "Can't open file:$!";
	}elsif($point eq "s"){
		open $IN,"awk -v OFS='' -v ORS='' -v up=$upstream -v dn=$downstream '{print \$1\"\\t\";print \$2-up\"\\t\"\$2+1+dn;for(i=4;i<=NF;i++){print \"\\t\"\$i}print \"\\n\"}' $bed|awk '\$2>0'|bedtools intersect -wao -a stdin -b $peakFile|" or die "Can't open file:$!";
	}elsif($point eq "e"){
		open $IN,"awk -v OFS='' -v ORS='' -v up=$upstream -v dn=$downstream '{print \$1\"\\t\";print \$3-1-up\"\\t\"\$3+dn;for(i=4;i<=NF;i++){print \"\\t\"\$i}print \"\\n\"}' $bed|awk '\$2>0'|bedtools intersect -wao -a stdin -b $peakFile|" or die "Can't open file:$!";
	}else{
		say "Please set correct parameter for -p";
	}
}
my $outStart=-1*$upstream;
my %hash;
for(my $i=$outStart;$i<=$downstream;$i++){
	$hash{$i}=0;
}
my $lineNum=0;
if(defined $heatmap){
	open OUT,">$heatmap" or die "Can't open file $heatmap:$!";
}
while(<$IN>){
	chomp;
	$lineNum+=1;
	my @split=split /\t/;
	if($split[$#split]==0){
		if(defined $heatmap){
			for(my $i=$outStart;$i<$downstream;$i++){
				print OUT "0\t";
			}
			say OUT "0";
		}
	}else{
		my $overlap_s=($split[1]>$split[1+$bedColumn])?$split[1]:$split[1+$bedColumn];
		my $overlap_e=($split[2]>$split[2+$bedColumn])?$split[2+$bedColumn]:$split[2];
		if($bedColumn>=6){
			if($split[5] eq "+"){
				if(defined $heatmap){
					for(my $j=$outStart;$j<$outStart+$overlap_s-$split[1];$j++){
						print OUT "0\t";
					}
				}
				for(my $i=$overlap_s;$i<$overlap_e;$i++){
					$hash{$outStart+$i-$split[1]}+=1;
					if(defined $heatmap){								
						print OUT "1\t";
					}
				}
				if(defined $heatmap){
					for(my $k=$outStart+$overlap_e-$split[1];$k<=$downstream;$k++){
						print OUT "0\t";
					}
					print OUT "\n";
				}
			}else{
				if(defined $heatmap){
					for(my $j=$outStart;$j<=$outStart+$split[2]-$overlap_e-1;$j++){
						print OUT "0\t";
					}
				}
				for(my $i=$overlap_s;$i<$overlap_e;$i++){
					$hash{$outStart+$split[2]-$i-1}+=1;
					if(defined $heatmap){
						print OUT "1\t";
					}
				}
				if(defined $heatmap){
					for(my $k=$outStart+$split[2]-$overlap_s;$k<=$downstream;$k++){
						print OUT "0\t";
					}
					print OUT "\n";
				}
			}
		}else{
			if(defined $heatmap){
				for(my $j=$outStart;$j<$outStart+$overlap_s-$split[1];$j++){
					print OUT "0\t";
				}
			}
			for(my $i=$overlap_s;$i<$overlap_e;$i++){
				$hash{$outStart+$i-$split[1]}+=1;
				if(defined $heatmap){
					print OUT "1\t";
				}
			}
			if(defined $heatmap){
				for(my $k=$outStart+$overlap_e-$split[1];$k<=$downstream;$k++){
						print OUT "0\t";
				}
				print OUT "\n";
			}
		}
	}
}
foreach my $pos(sort{$hash{$a}<=>$hash{$b}} keys %hash){
	say join "\t",$pos,$hash{$pos},$hash{$pos}/$lineNum;
}
sub usage{
print STDERR <<HELP 
Usage:	perl $0 -i <peak.bed> -b <*.bed> -p [c|s|e] -u <INT> -d <INT> >result.file 2>log
        Statistics the peak coverage across the input bed region. 
Author: Yumei Li, 2018-2-28
Revision: Yumei Li, 2018-3-6, Adding output the data to plot heatmap.
Revision: Yumei Li, 2018-8-16. Set -m as optional parameter.
Output: relative_position count frequency
    'i|peak'   FILE     The peak regions in bed format (Only the first three columns are used).
    'b|bed'    FILE     The targeted regions in bed format.( If it is in bed6 or bed6+ format, the strand will be taken into consideration)
    'p|point'  STRING   The reference point for the output coverage.[c for the region's center(default); s for start; e for end]
    'u|up'     INT      The upstream length from the reference point[default: 1000]
    'd|down'   INT      The downstream length from the reference point[default: 1000]
    'm|hmp'     FILE     The output file name for data to plot heatmap of the peak distribution (Optional).
    'help|h'            Print this help message    
HELP
}