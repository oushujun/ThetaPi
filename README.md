# ThetaPi
Nucleotide diversity (Tajima's Pi) estimation using population SNP data. Refer to Nei and Li (1979) and libsequence::PolySNP.c/ThetaPi.
Works for homozygous SNPs, heterozygous SNPs have not yet been taken take care of.
Applies missing rate screening for input data. Population size of a SNP is adjusted by the presence of individuals without missing.


Usage: $ perl ThetaPi.ori.pl infile missing_rate(o*) window_size(o) start_site(o) > outfile
\* o: optional
\- missing_rate: [0, 1], default: 0.50, sequence with missing rate higher than this value will be discarted.
\- window_size: int, >=1, default: 50000(bp), calculate averaged thetaPi/site in such size of unoverlaping windows.
\- start_site: int, >=0, default: 0, which means start calculating from the first site in the sequence. 1 from position 1.
eg: \$ perl ThetaPi.ori.pl sample -m0.90 -w100 -s1 > sample.pi
Input data format:
#Position Ind1 Ind2 Ind3 ... Indn
303 A A T ... A
499 T N T ... T
... ...
...
Output data format: start_loci end_loci Average_Pi_of_the_window

