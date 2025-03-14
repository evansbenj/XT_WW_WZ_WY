# Nigeria shared heterozygous positions

Based on heterozygosity, it appears that the samples Emmanuela collected were probably allotetraploids. The mtDNA confirms this for two of the samples because they have X. calcaratus mtDNA. But two others are X. tropicalis. One plausible explanation is that these samples are still allotetraploid but that they are hybrids with X. tropicalis mtDNA.

To test this I am quantifying the number of shared heterozygous positions in the samples. Because they are mapped to the trop genome, we expect the allotetraploids to have divergent sites between subgenomes that are called as heterozygous. 

I want to make a Venn diagram (or Euler diagram) that (minimally) shows the number of shared heterozygous positions among these four samples across the genome. Probably even better would be to do this including the other X. calcaratus sample from Cameroon.

Here' is a Perl script that quantifies this:
```perl
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
#  Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load StdEnv/2020 perl/5.30.2
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and counts sites that are heterozygous 
#  subsets of samples (e.g. all are het, or only some)
#  the easiest way to do this is to count homoz sites and assume the hets are
#  all identical (this can be checked too)

# to execute type Parse_tab_check_for_shared_heterozygotes.pl temp.tab 10110000011000000011 output.out 
# where 10110000011000000011 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) included, or (0) skipped

# I am doing this to check for evidence of polyploidy in the nigeria trop samples
# which have high nucleotide diversity

my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile1 = $ARGV[2];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE1, ">$outputfile1"))  {
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n";

my @sexes = split("",$ARGV[1]);

my @temp;
my $y;
my $x;
my $counter=0;
my $diverged=0;
my $number_of_included_individuals=0;
my $patterns=0;
my @patterns;
my @include;
my @unique_included_nucleotides;
my @individual_homoz_counter; # this counts how many positions are homoz in each individual and heteroz in all other individuals (not very useful)
my @individual_het_counter; # this counts how many positions are heterozygous in each individual, irrespective of
							#  the genotypes of other individuals 
my $homoz_indiv;
my $hets=0;
my $keyz;
my $other_counter=0;

my %Venn; # this will keep track of the overlap of heterozygous positions
# for 4 individuals these are the keys:
# 1,2,3,4 the first 4 values should be individual het positions for each sampl
# 12,13,14,23,24,34: the next 6 should be each possible pair
# 123,124,234: the next 3 should be each set of 3
# 1234: and the last 1 is shared hets for all individuals

for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 1){
		$number_of_included_individuals +=1;
	}	
}
print "This includes ",$number_of_included_individuals," individuals\n";
for ($x = 0 ; $x < $number_of_included_individuals ; $x++ ) {
	$patterns[$x]=0;
	$individual_homoz_counter[$x]=0;
}



while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if($temp[0] ne '#CHROM'){
		if($#temp ne (($#sexes+1)*2)+2){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n";
		}

		# parse the bases in all genotypes 
		@include=();
		$counter=0;
		for ($y = 3 ; $y <= $#temp; $y=$y+2 ) {
			if(($temp[$y] ne ".")&&($temp[$y+1] ne ".")){
				if($sexes[$counter] == 1){
					push(@include, $temp[$y]);
					push(@include, $temp[$y+1]);
				}	
			}
			$counter+=1;
		}	
		# OK I should have all the bases loaded for included genotypes
		
		@unique_included_nucleotides = uniq @include;
		#print "yo ",$#unique_included_nucleotides," ",$#include,"\n";
		if(($#unique_included_nucleotides == 1)&&($#include == ($number_of_included_individuals*2)-1)){
			# we can count homoz and het genotypes because all individuals have data
			# and there is at least one heterozygous position and the variation is biallelic
				$hets=0;
				$diverged=0;
				$keyz="";
				$other_counter=0;
				# go through each position
				for ($x = 0 ; $x <= $#include ; $x=$x+2 ) {
					$other_counter+=1;
					#print "hey ",$include[$x]," ",$include[$x+1],"\n";
					if($include[$x] ne $include[$x+1]){
						$hets+=1;
						$keyz=$keyz.$other_counter;
					}
					else{ # this individual is homoz
						$homoz_indiv = $x/2;
					}
				}
				# now update the hash
				$Venn{$keyz}+=1;
				if($hets == 3){ # tabulate which individual is homozygous while the others are all het
					$individual_homoz_counter[$homoz_indiv]+=1;
				}
				$patterns[$hets-1]+=1;
		}
	}
} # end while
#print "@patterns\n";
#print "@individual_homoz_counter\n";

my $number_of_keys = scalar keys %Venn;
my $countr = 0;

foreach my $key (sort keys %Venn){
	if($key ne ""){
		@temp=split "",$key;
		$countr+=1;
		#print "key $key hello @temp hi $temp[$#temp]\n";
		print "`";
		for ($y =0 ; $y < $#temp; $y++) {
			print $temp[$y]."&";
		}
		print $temp[$#temp],'`=';
		if($countr < ($number_of_keys-1)){
			print $Venn{$key},",";
		}
		else{
			print $Venn{$key},"\n";
		}	
		#print OUTFILE1 $key," ",$Venn{$key},"\n";
	}
}	
print OUTFILE1 "@patterns\n";
close OUTFILE1;
```


Then this output could be copied and pasted into this R code to print an UpSet plot:
```R
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2021_XT_sex_linked_markers/2024_Parsetab_sharedhets")
library("UpSetR")
options(scipen=999)
# example of expression input
#expressionInput <- c(one = 2, two = 1, three = 2, `one&two` = 1, `one&three` = 4, 
#                     `two&three` = 1, `one&two&three` = 2)

# 1:EUA0331	(NG2)
# 2:EUA0333	(NG3)
# 3:EUA0334 (NG4)
# 4:EUA0335 (NG5)
# 5:Xcal  (Cal)

expressionInput <- c(`NG2`=1214684,
                     `NG2&NG3`=1847532,
                     `NG2&NG3&NG4`=6969040,
                     `NG2&NG3&NG4&NG5`=49389,
                     `NG2&NG3&NG4&NG5&Cal`=126092,
                     `NG2&NG3&NG4&Cal`=17871976,
                     `NG2&NG3&NG5`=12877,
                     `NG2&NG3&NG5&Cal`=19806,
                     `NG2&NG3&Cal`=2898475,
                     `NG2&NG4`=946774,
                     `NG2&NG4&NG5`=15152,
                     `NG2&NG4&NG5&Cal`=11110,
                     `NG2&NG4&Cal`=895576,
                     `NG2&NG5`=13173,
                     `NG2&NG5&Cal`=5940,
                     `NG2&Cal`=595158,
                     `NG3`=1424512,
                     `NG3&NG4`=1194663,
                     `NG3&NG4&NG5`=18085,
                     `NG3&NG4&NG5&Cal`=15841,
                     `NG3&NG4&Cal`=1367525,
                     `NG3&NG5`=14828,
                     `NG3&NG5&Cal`=7599,
                     `NG3&Cal`=768914,
                     `NG4`=10795399,
                     `NG4&NG5`=833968,
                     `NG4&NG5&Cal`=35534,
                     `NG4&Cal`=1355129,
                     `NG5`=1514606,
                     `NG5&Cal`=45505,
                     `Cal`=9704927)


upset(fromExpression(expressionInput), 
      #keep.order = T,
      order.by = "freq",
      show.numbers = "no",
      set_size.show = F,
      sets = c("NG2", "NG3", "NG4","NG5","Cal"))


# you need to save the plot from the RStudio plot window
# a size of 6.5 x 9 works well

```
