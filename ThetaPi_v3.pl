##This program is for ThetaPi calculation, the principle equation could refer to libsequence, under the tag of PolySNP.c/ThetaPi
#(http://molpopgen.github.io/libsequence/doc/html/classSequence_1_1PolySNP.html#a9c26d27c7b0aaf1faf859272871b3cf5)	.
##NOTICE: this program will only for homozygous SNP calculation, heterozygous situations have not yet been taken care of.

#!/usr/bin/perl -w
use strict;

my $usage="This program is written for thetaPi calculation for homozygous SNP data
		Usage: \$ perl ThetaPi.ori.pl infile missing_rate(o*) window_size(o) start_site(o) > outfile
			* o: optional
			- missing_rate: [0, 1], default: 0.50, sequence with missing rate higher than this value will be discarted.
			- window_size: int, >=1, default: 50000(bp), calculate averaged thetaPi/site in such size of unoverlaping windows. 
			- start_site: int, >=0, default: 0, which means start calculating from the first site in the sequence. 1 from position 1.
		eg: \$ perl ThetaPi_v3.pl sample -m0.90 -w100 -s1 > sample.pi
		Input data format:
			#Chr Position Ind1 Ind2 Ind3 ... Indn
			chr01 303 A A T ... A
			chr01 499 T N T ... A
			... ...
			...
		 Output data format: chromosome start_loci end_loci Average_Pi_of_the_window

		Author: Shujun Ou (oushujun\@msu.edu)
		Versions: 
				v3.0 Shujun Ou 05/16/2017 Adaptive to other data format; add support to diploid genotypes instead of haploid.
				v2.4 Shujun Ou 11/26/2014 Improve the end position obtaining scheme
				v2.3 Shujun Ou 02/19/2014 Add the -s para, improve the start and end positions obtaining process
				v2.2 Shujun Ou 02-15-2013
				v1.0 Shujun Ou 2012/10/30 The first version
		";
if ($#ARGV<0){die "$usage"}

my $window=1000;
my $missing=0.50; #default window size is 50kb and missing cutoff is 0.50
my $firstloc=0; 
foreach (@ARGV){  #customer windos size and missing rate reset
	if (m/^-m([0-9.]+)$/i){$missing=$1}
	elsif (m/^-w([0-9]+)$/i){$window=$1}
	elsif (m/^-s([0-9]+)$/i){$firstloc=$1}
	}

##Store the infile into hash, with location as keys
open FILE, "<$ARGV[0]" or die "No such file: $ARGV[0]\n";
my $chr="unknown";
my %file;
foreach (<FILE>) {
	next if /^#/; #skip annotations
	next if /Position/i; #skip annotations
	next if /loci/i; #skip annotations
	$chr=$1 if s/^(\D+\S+)\s+//i; #remove chromosome name from a line
	my $key=0;
	$key=$1 if s/(^[0-9]+)\s+//; #This is the position of the SNP
	$file{$key}=$_;
	}
close FILE;

##get the first location, the last location, and the ploidy level of the infile
my $start=0;
my $end=0;
my $ploidy=1;
($start, $end)=(sort{$a<=>$b}(keys %file))[0,-1];
my $sample=(split /\s+/, $file{$start})[0];
$ploidy=($sample=~tr/\///)+1;

die "The input file may be malformatted\n" unless $end>0 and $start>0;

print "Chr\tFrom\tTo\tThetaPi/site\n";
$firstloc=$start if $firstloc==0; #start Pi calculation from the first SNP if not specified the start site
my $lastloc=$end+$window; #this is the position of the SNP

my $bed=$firstloc;
my $high=$bed+$window;
while ($high<=$lastloc) {
	my ($window_pi, $c, $exist);
	$window_pi=$c=$exist=0;
	my $win_start=$bed;
	my $win_end=$high;
	while ($bed<$high) {
		my $site_pi=0;
		if (exists $file{$bed}) {
			$exist++; #loci existance control
			my $site;
			$site=$file{$bed};
			$site=~s/chr[0-9]+//gi;
			$site=~s/([0-9]+)//g;
			$site=~s/\s+//g;
			$site=~s/N/-/gi;
			$site=~s/\?/-/g;
			$site=~s/\./-/g;
			my ($count_A, $count_T, $count_C, $count_G, $count_N);
			$count_A=$count_T=$count_C=$count_G=$count_N=0;
			$count_A=$site=~tr/Aa//;
			$count_T=$site=~tr/Tt//;
			$count_C=$site=~tr/Cc//;
			$count_G=$site=~tr/Gg//;
			$count_N=$site=~tr/-//;

			#reduce polyploid to haploid; heterozygosity is also taken care of.
			$count_A=$count_A/$ploidy;
			$count_T=$count_T/$ploidy;
			$count_C=$count_C/$ploidy;
			$count_G=$count_G/$ploidy;
			$count_N=$count_N/$ploidy;

			#If there is any missing value, reduce the population size for that site
			my $n=$count_A+$count_T+$count_C+$count_G;
			my $miss=$count_N/($count_N+$n);
			if ($miss>$missing){next} #missing rate control
			if ($n>1){ #n is the number of bases in a line, missing not included
				$c++;
				my ($kA, $kT, $kC, $kG);
				$kA=$kT=$kC=$kG=0;
				$kA=($count_A*($count_A-1))/($n*($n-1));
				$kT=($count_T*($count_T-1))/($n*($n-1));
				$kC=($count_C*($count_C-1))/($n*($n-1));
				$kG=($count_G*($count_G-1))/($n*($n-1));
				$site_pi=1-($kA+$kT+$kC+$kG);
				$window_pi+=$site_pi;
					}
			}
	} continue {
	$bed++;
		}

	if ($c>0){ #c is the number of data lines in a window. The line with all missing is not included.
		$window_pi=sprintf("%.4f", $window_pi/$c); #get the mean pi value in the window
		print "$win_start\t$win_end\t$window_pi\n";
		} elsif ($exist > 0) { #$exist is all lines in a window, lines with all missing are included
		print "$win_start\t$win_end\t0\n";
		} else {
		print "$win_start\t$win_end\t0\n"; #NA lines calculated as 0, for data processing convenience
		}
	$bed=$high;
	$high+=$window;
	}
