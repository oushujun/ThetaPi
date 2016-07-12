##This program is worte for ThetaPi calculation, the principle equation could refer to libsequence, under the tag of PolySNP.c/ThetaPi.
##NOTICE: this program will only for homozygous SNP calculation, heterozygous situations have not yet been taken care of.

#!/usr/bin/perl -w
use strict;

my $usage="This program is written for thetaPi calculation for homozygous SNP data
		Usage: \$ perl ThetaPi.ori.pl infile missing_rate(o*) window_size(o) start_site(o) > outfile
			* o: optional
			- missing_rate: [0, 1], default: 0.50, sequence with missing rate higher than this value will be discarted.
			- window_size: int, >=1, default: 50000(bp), calculate averaged thetaPi/site in such size of unoverlaping windows. 
			- start_site: int, >=0, default: 0, which means start calculating from the first site in the sequence. 1 from position 1.
		eg: \$ perl ThetaPi.ori.pl sample -m0.90 -w100 -s1 > sample.pi
		Input data format:
			#Position Ind1 Ind2 Ind3 ... Indn
			303 A A T ... A
			499 T N T ... T
			... ...
			...
		 Output data format: start_loci end_loci Average_Pi_of_the_window

		Author: Shujun Ou (oushujun\@msu.edu)
		Versions: 
				v2.4 Shujun Ou 11/26/2014 Improve the end position obtaining scheme
				v2.3 Shujun Ou 02/19/2014 Add the -s para, improve the start and end positions obtaining process
				v2.2 Shujun Ou 02-15-2013
				v1.0 Shujun Ou 2012/10/30 The first version
		";
if ($#ARGV<0){die "$usage"}

my $window=50000;
my $missing=0.50; #default window size is 50kb and missing cutoff is 0.50
my $firstloc=0;
foreach (@ARGV){  #customer windos size and missing rate reset
	if (m/-m([0-9.]+)/i){$missing=$1}
	elsif (m/-w([0-9]+)/i){$window=$1}
	elsif (m/-s([0-9]+)/i){$firstloc=$1}
	}

print "From\tTo\tThetaPi/site\n";

##get the first location and the last location of the infile
my $end;
my $i=1;
die "No such file: $ARGV[0]\n" unless open TEST, "< $ARGV[0]";
close TEST;
while (){
	$end=`tail -$i $ARGV[0]`;
	$end=~s/chr[0-9]+\s//gi;
	$end=(split /\s+/,$end)[0]; #this is the position of the SNP
	last if $end>0;
	$i++;
	}	
die "no last loc found\n" unless $end>0;
my $lastloc=$end+$window; #this is the position of the SNP

my $start=`head -2 $ARGV[0]`;
$start=~s/^#.*\n?//;
$start=~s/.*Position.*\n?//i;
$start=~s/chr[0-9]+\s//gi;#remove 'chr01' from a line

if ($firstloc==0){
	$firstloc=(split /\s+/, $start)[0]; #start Pi calculation from the first SNP(second row of the file)
}

##Transfer the infile into hash, with location as keys
open FILE2, "<$ARGV[0]";
my %file;
foreach (<FILE2>) {
	if (/^#/){next}
	if (/Position/i){next}
	s/^chr[0-9]+\s//gi;
	my $key;
if (s/(^[0-9]+)//){$key=$1}; #This is the position of the SNP
	$file{$key}=$_;
	}
close FILE2;

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
			my ($count_A, $count_T, $count_C, $count_G, $count_N);
			$count_A=$count_T=$count_C=$count_G=$count_N=0;
			$count_A=$site=~tr/Aa//;
			$count_T=$site=~tr/Tt//;
			$count_C=$site=~tr/Cc//;
			$count_G=$site=~tr/Gg//;
			$count_N=$site=~tr/-//;
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
