# This file looks for W-linked, Z-linked, and Y-linked SNPs
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
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and searches for sites that are consistent with W-, Y-, and Z-linked inheritance

# A W-specific variant is present in all females and none or some males; 
	# it can be homozygous or heterozygous in females
	# it can have any genotype in males because the W and Y variation may not overlap (so a W-linked variant could be linked to a Y-linked variant); however, this makes it have no pattern
	# instead we will define it as never homozygous in males
# A Y-specific variant is present in some males and never in females
	 # it can be heterozygous in males but never homozygous 
	 # that SNP is never in females
# A Z-specific variant is present some males and some females; 
	# it can be homozygous or heterozygous in males 
	# it can be heterozygous but never homozygous in females

# to execute type Parse_tab.pl inputfile.tab 1111100110000111100011100110010100002200 interesting_sites.out proportion
# where 1111100110000111100011100110010100002200 refers to whether or not each individual in the ingroup 
# in the vcf file is (0) male, (1) female, and or (2) skipped

# to skip the tadpoles use this: 111111000000122221
# to include the tadpoles use this: 111111000000111001

# Example for XT_WGS (ignoring tads)
# perl Parse_tab_XT.pl XT_Chr7_1_22416950.tab 111111000000122221 interesting_sites.out

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
print OUTFILE1 "CHR\tPOS\tTYPE\thet_females\thet_males\tn_FEMs\tn_MALS\n";



my @sexes = split("",$ARGV[1]);

my @males=();
my @females=();
my @MFcombined=();
my @temp;
my @unique_male_nucleotides;
my @unique_female_nucleotides;
my @unique_MFcombined_nucleotides;
my @W_linked_nucleotides;
my @Y_linked_nucleotides;
my @Z_linked_nucleotides;
my $not_W=0;
my $not_Y=0;
my $not_Z=0;
my $y;
my $x;
my $counter=0;
my $diverged=0;
my $diverged_2=0;
my $number_of_male_individuals_genotyped=0;
my $number_of_female_individuals_genotyped=0;

for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 0){
		$number_of_male_individuals_genotyped +=1;
	}	
}	
for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s) and  ", $number_of_male_individuals_genotyped," males\n";

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if($temp[0] ne '#CHROM'){
		if($#temp ne (($#sexes+1)*2)+2){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n";
		}

		# parse the bases in all genotypes in each sex
		@males=();
		@females=();
		@MFcombined=();
		$counter=0;
		for ($y = 3 ; $y <= $#temp; $y=$y+2 ) {
			if(($temp[$y] ne ".")&&($temp[$y+1] ne ".")){ # allowing asterisks (*), which are deletions
				if($sexes[$counter] == 0){
						push(@males, $temp[$y]);
						push(@males, $temp[$y+1]);
				}
				elsif($sexes[$counter] == 1){
					push(@females, $temp[$y]);
					push(@females, $temp[$y+1]);
				}
				push(@MFcombined, $temp[$y]);	
				push(@MFcombined, $temp[$y+1]);	
			}
			$counter+=1;
		}	
		# OK I should have all the bases loaded for non-missing genotypes for each male and each female
		# and combined as well
		
		# find out what and how many unique nucleotides are in each sex
		@unique_male_nucleotides = uniq @males;
		@unique_female_nucleotides = uniq @females;
		# find out how many unique nucleotides there are in the combined genotypes
		@unique_MFcombined_nucleotides = uniq @MFcombined;
		
		#print @females," ",@males,"\n";
		#print $#unique_male_nucleotides," ",$#unique_female_nucleotides,"\n";
		# looks fine
		if(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides != -1)){
			# this means that there is at least one genotype in each sex
			# we can compare homoz and het genotypes because both sexes have data
			
			# First check for W-linked SNPs; 
				# all males homozygous 
				# all females are heterozygous or homozygous for the non-male variant
			if(($#unique_male_nucleotides == 0)&&($#unique_female_nucleotides > 0)){
				# the males are all homozygous for $unique_male_nucleotides[0]
				# now check if any females are homozygous for this variant 
				$not_W=0;
				$diverged=0;
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) {
					if(($females[$x] eq $females[$x+1])&&($females[$x] eq $unique_male_nucleotides[0])){
						# this female is homoz for the same variant as in the males
						# which means this position is not W-linked
						$not_W=1;
					}
					elsif($females[$x] ne $females[$x+1]){
						$diverged+=1;  # this is the number of heterozygous females
					}
					
				}
				if($not_W==0){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tW_linked\t",$diverged,"\t0\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# W-linked variation
				}								
			} # end of check for W-linked SNPs
			
			# Now check for Y-linked SNPs; 
				# only found in males in heterozygous genotypes
				# never in females 
				# never homozygous in males
			if(($#unique_male_nucleotides > 0)&&($#unique_female_nucleotides == 0)){
				# the males have variation and the females do not
				# now check if any females are homozygous for this variant 
				$not_Y=0;
				$diverged=0;
				for ($x = 0 ; $x <= $#males ; $x=$x+2 ) {
					if($males[$x] eq $males[$x+1]){
						# this male is homoz, and potentially ZZ
						# check if this is not the variant that is found in 
						# females, which would mean that
						# this position is not Y-linked (because
						# we can't be homoz for YY)
						if($males[$x] ne $unique_female_nucleotides[0]){
							$not_Y=1;
						}	
					}
					else{
						$diverged+=1; # this is the number of heterozygous males
					}
				}
				if($not_Y==0){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tY_linked\t0\t$diverged\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# W-linked variation
				}								
			} # end of check for Y-linked SNPs

			# Now check for Z-linked SNPs; 
				# never homozy in females 
			if($#unique_female_nucleotides > 0){
				# the females have variation
				# now check if any females are homozygous for this variant 
				$not_Z=0;
				$diverged=0;
				$diverged_2=0;
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) {
					if($females[$x] eq $females[$x+1]){
						# this female is homoz, so this
						# is not a Z-linked SNP
						$not_Z=1;
					}
					else{
						$diverged+=1; # this is the number of heteroz females
					}
				}
				# check  how many males are heteroz
				for ($x = 0 ; $x <= $#males ; $x=$x+2 ) {
					if($males[$x] ne $males[$x+1]){
						# this male is het, possibly because he is XY
						$diverged_2+=1; # this is the number of heteroz females
					}
				}
				
				if($not_Z==0){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tZ_linked\t$diverged\t$diverged_2\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# W-linked variation
				}								
			} # end of check for Z-linked SNPs 
		} # end of check that there is at least one genotype in each sex
	} # end of check to see if we are at the first line	
} # end while	
close OUTFILE1;



```
