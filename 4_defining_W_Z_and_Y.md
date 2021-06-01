# Defining W, Z, and Y

The three tads we have are siblings with the following genotypes: WW, WY, and WZ. The mother of the tads was the daughter of the father of the tads. The father was WY. The mother was WZ. So the mother should have the same W chromosome as the father.  If there is recombination on the distal tips of the W and Z and/or the distal tips of the W and the Y, then the WW tad might have heterozygosity in this region.

# Workflow

First step is to evaluate the genotypes that are mapped to v10.  Probably the easiest way to do this is to first export them to a tab-delimited file and then write a perl script to generate a synthetic W, Z, and Y chromosome.  This would be interesting for estimation of genetic distances (prediction is that Z and Y are more similar to each other than either is to the W) and also to explore where differences are distributed (prediction is that differences between Z and Y should be highest immediately surrounding the sex-linked region; differences between the Z/Y and the W should be more broadly distributed).

The coordinates of these chromosomes will not match the reference due to insertion deletion events. For this reason, it probably would be wise to generate paml files directly from the tab delimited file for each gene. This would allow me to combine exons for each gene, check for and quantify frameshifts, and then use paml to quantify NS and S changes.

I think I will begin by exporting a tab-delimited file and then using a Perl script to identify sites in the WW individual and subtracting this variation to generate a Z or Y in the WZ and WY individuals. I'd have to develop some rules for resolving heterozygosity in the W. There almost certainly will be heterozygous positions due to repetitive regions.  Probably the most conservative way to deal with this is to insert Ns for anything with W heterozygosity.  

I could also use the Z- and Y-specific kmer assemblies to verify SNPs. I could generate two types of data: one with each chr based on comparisons only between the called genotypes from the WW, WZ, and WY individuals and the other based on this plus verification with the Z- and Y-specific kmer assemblies.

To accomplish the second goal of analyzing NS and S variation in alleles on each chr, I could export the coding regions of each gene (with exons concatenated) as a vcf file, convert to tab, and then use the same script to output paml formatted files for analysis.

Things to check on
* does paml accept Ns
* can I export multiple regions into a concatenated vcf file

# subset a vcf file
```
module load bcftools/1.11
bcftools view -r Chr7:1-30000000 ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_filtered.vcf.gz -o ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_filtered_pos1_30Mb.vcf
```

# Make a tab file from a subsetted vcf file
```
module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.16
vcf-to-tab < ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_filtered_pos1_30Mb.vcf > ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_filtered_pos1_30Mb.tab
```

# Make a phy file out of a tab file
I wrote a perl script for this:
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );


#  This program reads in a tab delimited genotype file generated
#  from vcftools that has a WW, WZ and WY individual
#  and prints out a fasta file with an interpretation for the 
#  W, Z, and Y chromosome sequence

# Data should be from one chromosome.

# the first individual is XT10_WZ
# the second individual is XT11_WW
# the third individual is XT7_WY

# to execute type Makes_W_Z_Y_from_tab.pl inputfile.tab output.fasta  


my $inputfile = $ARGV[0];
my $outputfile2 = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}


unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2\n";
	exit;
}
print "Creating output file: $outputfile2\n";

my $W;
my $Z;
my $Y;
my $allele_1;
my $allele_2;
my $allele_3;
my $allele_4;
my $allele_5;
my $allele_6;
my @temp;
my @lengths;
my $max;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		print $temp[8],"\n";
		# load the alleles
			# load each of the six alleles
				$allele_1 = $temp[3];
				$allele_2 = $temp[4];
				$allele_3 = $temp[5];
				$allele_4 = $temp[6];
				$allele_5 = $temp[7];
				$allele_6 = $temp[8];
			# now see which is the longest
				@lengths=();
				push(@lengths,length($allele_1));
				push(@lengths,length($allele_2));
				push(@lengths,length($allele_3));
				push(@lengths,length($allele_4));
				push(@lengths,length($allele_5));
				push(@lengths,length($allele_6));
				$max = max @lengths;
			# add the data to the chrs
				if($max == 1){ # this is just a SNP
					####################
					# first consider W
					####################
					if($allele_3 eq $allele_4){ # the W genotype is homozygous
						if(($allele_3 eq 'A')||($allele_3 eq 'C')||($allele_3 eq 'G')||($allele_3 eq 'T')){
							$W = $W.$allele_3;
						}
						elsif($allele_3 eq '*'){
							$W = $W.'-';
						}
						else{
							$W = $W.'N';
						}
					}
					elsif($allele_3 ne $allele_4){ # the W genotype is heterozygous
						if((($allele_3 eq 'A')||($allele_3 eq 'C')||($allele_3 eq 'G')||($allele_3 eq 'T'))
							&&
							(($allele_4 eq 'A')||($allele_4 eq 'C')||($allele_4 eq 'G')||($allele_4 eq 'T'))
							){ # both alleles are nucleotides
							# insert an IUPAC symbol
							if((($allele_3 eq 'A')&&($allele_4 eq 'C'))||(($allele_3 eq 'C')&&($allele_4 eq 'A'))){
								$W = $W.'M';
							}	
							elsif((($allele_3 eq 'A')&&($allele_4 eq 'T'))||(($allele_3 eq 'T')&&($allele_4 eq 'A'))){
								$W = $W.'W';
							}	
							elsif((($allele_3 eq 'A')&&($allele_4 eq 'G'))||(($allele_3 eq 'G')&&($allele_4 eq 'A'))){
								$W = $W.'R';
							}	
							elsif((($allele_3 eq 'C')&&($allele_4 eq 'G'))||(($allele_3 eq 'G')&&($allele_4 eq 'C'))){
								$W = $W.'S';
							}	
							elsif((($allele_3 eq 'C')&&($allele_4 eq 'T'))||(($allele_3 eq 'T')&&($allele_4 eq 'C'))){
								$W = $W.'Y';
							}	
							elsif((($allele_3 eq 'G')&&($allele_4 eq 'T'))||(($allele_3 eq 'T')&&($allele_4 eq 'G'))){
								$W = $W.'K';
							}	
						}
						else{ # a het with a dot means unknown, or any other weirdness will be set to be unknown
							$W = $W.'N';
						}

					}	
					####################
					# Now consider the Z
					####################
					if($allele_1 eq $allele_2){ # the Z genotype is homozygous
						# Confirm that the alleles are the same as the W and only use them if they are
						if($allele_1 eq $allele_3){
							if(($allele_1 eq 'A')||($allele_1 eq 'C')||($allele_1 eq 'G')||($allele_1 eq 'T')){
								$Z = $Z.$allele_1;
							}
							elsif($allele_1 eq '*'){
								$Z = $Z.'-';
							}
							else{
								$Z = $Z.'N';
							}
						}
						else{ # we don't know what the Z is because both alleles are different from the W
							$Z = $Z.'N';
						}
					}
					elsif($allele_1 ne $allele_2){ # the Z genotype is heterozygous
						if((($allele_1 eq 'A')||($allele_1 eq 'C')||($allele_1 eq 'G')||($allele_1 eq 'T'))
							&&
							(($allele_2 eq 'A')||($allele_2 eq 'C')||($allele_2 eq 'G')||($allele_2 eq 'T'))
							){ # both Z alleles are nucleotides
							# check if only one is the same as the W
							if((($allele_1 ne $allele_3)||($allele_1 ne $allele_4))&&
								(($allele_2 eq $allele_3)||($allele_2 eq $allele_4))){
								$Z = $Z.$allele_1;
							}
							elsif((($allele_2 ne $allele_3)||($allele_2 ne $allele_4))&&
								(($allele_1 eq $allele_3)||($allele_1 eq $allele_4))){
								$Z = $Z.$allele_2;
							}
							else{ # the WZ genotype is completely different from the WW genotype
								# so we are unable to make a call
								$Z = $Z.'N';
							}
						}
						elsif(($allele_1 eq '*')||($allele_2 eq '*')){ # one allele on the Z is a gap
							if((($allele_1 ne $allele_3)||($allele_1 ne $allele_4))&&
								(($allele_2 eq $allele_3)||($allele_2 eq $allele_4))){
								if($allele_1 eq '*'){
									$Z = $Z.'-'; # the Z has a gap and the W does not
								}
								else{
									$Z = $Z.$allele_1; # the W has a gap and the Z does not
								}	
							}
							elsif((($allele_2 ne $allele_3)||($allele_2 ne $allele_4))&&
								(($allele_1 eq $allele_3)||($allele_1 eq $allele_4))){
								if($allele_2 eq '*'){
									$Z = $Z.'-'; # the Z has a gap and the W does not
								}
								else{
									$Z = $Z.$allele_2; # the W has a gap and the Z does not
								}	
							}
							else{ # the WZ genotype is completely different from the WW genotype
								$Z = $Z.'N';
							}
						}
						else{ # a homozygous dot means unknown, or any other weirdness will be set to be unknown
							$Z = $Z.'N';
						}
					}	

					####################
					# Now consider the Y
					####################
					if($allele_5 eq $allele_6){ # the Y genotype is homozygous
						# Confirm that the alleles are the same as the W and only use them if they are
						if($allele_5 eq $allele_3){
							if(($allele_5 eq 'A')||($allele_5 eq 'C')||($allele_5 eq 'G')||($allele_5 eq 'T')){
								$Y = $Y.$allele_5;
							}
							elsif($allele_5 eq '*'){
								$Y = $Y.'-';
							}
							else{
								$Y = $Y.'N';
							}
						}
						else{ # we don't know what the Y is because both alleles are different from the W
							$Y = $Y.'N';
						}
					}
					elsif($allele_5 ne $allele_6){ # the Y genotype is heterozygous
						if((($allele_5 eq 'A')||($allele_5 eq 'C')||($allele_5 eq 'G')||($allele_5 eq 'T'))
							&&
							(($allele_6 eq 'A')||($allele_6 eq 'C')||($allele_6 eq 'G')||($allele_6 eq 'T'))
							){ # both Y alleles are nucleotides
							# check if only one is the same as the W
							if((($allele_5 ne $allele_3)||($allele_5 ne $allele_4))&&
								(($allele_6 eq $allele_3)||($allele_6 eq $allele_4))){
								$Y = $Y.$allele_5;
							}
							elsif((($allele_6 ne $allele_3)||($allele_6 ne $allele_4))&&
								(($allele_5 eq $allele_3)||($allele_5 eq $allele_4))){
								$Y = $Y.$allele_6;
							}
							else{ # the WY genotype is completely different from the WW genotype
								# so we are unable to make a call
								$Y = $Y.'N';
							}
						}
						elsif(($allele_5 eq '*')||($allele_6 eq '*')){ # one allele on the Y is a gap
							if((($allele_5 ne $allele_3)||($allele_5 ne $allele_4))&&
								(($allele_6 eq $allele_3)||($allele_6 eq $allele_4))){
								if($allele_5 eq '*'){
									$Y = $Y.'-'; # the Y has a gap and the W does not
								}
								else{
									$Y = $Y.$allele_5; # the W has a gap and the Y does not
								}	
							}
							elsif((($allele_6 ne $allele_3)||($allele_6 ne $allele_4))&&
								(($allele_5 eq $allele_3)||($allele_5 eq $allele_4))){
								if($allele_6 eq '*'){
									$Y = $Y.'-'; # the Y has a gap and the W does not
								}
								else{
									$Y = $Y.$allele_6; # the W has a gap and the Y does not
								}	
							}
							else{ # the WY genotype is completely different from the WW genotype
								$Y = $Y.'N';
							}
						}
						else{ # a homozygous dot means unknown, or any other weirdness will be set to be unknown
							$Y = $Y.'N';
						}
					}	
				}	# end of check for one nucleotide in all genotypes
	} # end if to check for header of input file
} # end while	
close DATAINPUT;				

# print out the lengths of each chr to see if we missed anything
print "Length of W :",length($W),"\n";
print "Length of Z :",length($Z),"\n";
print "Length of Y :",length($Y),"\n";

# OK print out the fasta file

print OUTFILE2 "3 ",length($W),"\n";
print OUTFILE2 "W_chr     ";
print OUTFILE2 $W,"\n";
print OUTFILE2 "Z_chr     ";
print OUTFILE2 $Z,"\n";
print OUTFILE2 "Y_chr     ";
print OUTFILE2 $Y,"\n";
close OUTFILE2;
```

here is an example command:
```
./Makes_W_Z_Y_from_tab.pl ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_pos8Mb_12Mb.tab ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_pos8Mb_12Mb_W_Z_Y.phy
```
#Calculate distances
```
module load StdEnv/2020  intel/2020.1.217 phylip/3.698
dnadist ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_pos8Mb_12Mb_W_Z_Y.phy ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_pos8Mb_12Mb_W_Z_Y.dnadist
```
