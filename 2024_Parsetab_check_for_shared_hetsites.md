# Path
Filtered tab files (and gvcfs) are here:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_after_filtering
```

# Check for shared het sites
I want to quantify the number of sites that have het genotypes shared in all samples from Nigeria, or a portion of them

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
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and counts sites that are heterozygous 
#  subsets of samples (e.g. all are het, or only some)
#  the easiest way to do this is to count homoz sites and assume the hets are
#  all identical (this can be checked too)

# to execute type Parse_tab_check_for_shared_heterozygotes.pl temp.tab 00110000011000000000 output.out 
# where 00110000011000000000 refers to whether or not each individual in the ingroup 
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
my $hets=0;

for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 1){
		$number_of_included_individuals +=1;
	}	
}
print "This includes ",$number_of_included_individuals," individuals\n";
for ($x = 0 ; $x < $number_of_included_individuals ; $x++ ) {
	$patterns[$x]=0;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if($temp[0] ne '#CHROM'){
		if($#temp ne (($#sexes+1)*2)+2){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n";
		}

		# parse the bases in all genotypes in each sex
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
				$hets=0;
				# there is at least one heterozygous position and the variation is biallelic
				# now check for the extreme case where all females are heterozygous and all males are homoz
				# this is rare because we expect some genotypes in females to be undercalled, even in sex-linked regions
				$diverged=0;
				for ($x = 0 ; $x <= $#include ; $x=$x+2 ) {
					#print "hey ",$include[$x]," ",$include[$x+1],"\n";
					if($include[$x] ne $include[$x+1]){
						$hets+=1;
					}
				}
				$patterns[$hets-1]+=1;
		}
	}
} # end while
print "@patterns\n";
print OUTFILE1 "@patterns\n";
close OUTFILE1;
```

