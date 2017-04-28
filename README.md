# What's this
This is a PERL script for nucleotide diversity (Tajima's Pi) estimation using population SNP data. Concepts and equations refer to Nei and Li (1979) and libsequence::PolySNP.c/ThetaPi.

Works for homozygous SNPs, heterozygous SNPs have not yet been taken care of.

Applies missing rate screening for input data. Population size of a SNP is adjusted by the presence of individuals without missing.

# How to use
$ perl ThetaPi.ori.pl infile missing_rate(o*) window_size(o) start_site(o) > outfile

	* o: optional

	- missing_rate: [0, 1], default: 0.50, sequence with missing rate higher than this value will be discarted.

	- window_size: int, >=1, default: 50000(bp), calculate averaged thetaPi/site in such size of unoverlaping windows.

	- start_site: int, >=0, default: 0, which means start calculating from the first site in the sequence. 1 from position 1.

eg: $ perl ThetaPi_v2.pl example.sd1.gt.txt -m0.8 -w100 -s0 > example.sd1.gt.txt.pi

# Input and Output
For the genotype file input, lines started with '#' is the annotation line, which will be ignored by the script. Each row is for one locus, either polymorphic or invariable, with the locus coordinate stated at the first position of the row. Between rows, the coordinates do not need to be continuous but must be in the ascending numerical order. Genotypes are denoted using the IUPAC rule with one letter denotes the genotype. Each column represents the sequence of one individual. You can put as many individuals in this table and all of them together will be considered a single population.

Input file looks like (please see the example.sd1.gt.txt file as an example):

	#Position Ind1 Ind2 Ind3 ... Indn

	303 A A T ... A

	499 T N T ... T

	... ...

	...

Output data format: 

	start_loci end_loci Average_Pi_of_the_window

# Reference
Hu, Bin, Wei Wang, Shujun Ou, Jiuyou Tang, Hua Li, Ronghui Che, Zhihua Zhang, et al. “Variation in NRT1.1b Contributes to Nitrate-Use Divergence between Rice Subspecies.” Nat Genet 47, no. 7 (2015): 834-8.

Xu, Fan., Jun. Fang, Shujun Ou, Shaopei Gao, Fengxia Zhang, Lin. Du, Yunhua Xiao, et al. “Variations in CYP78A13 Coding Region Influence Grain Size and Yield in Rice.” Plant Cell Environ 38, no. 4 (2014): 800-11.

libsequence::PolySNP.c/ThetaPi. https://molpopgen.github.io/libsequence/doc/html/classSequence_1_1PolySNP.html#a9c26d27c7b0aaf1faf859272871b3cf5

Wikipedia. https://en.wikipedia.org/wiki/Nucleotide_diversity

