# Artifical allotetraploid

One possibility is that X. calcaratus is derived from allotetraploidization between the Liberia strain and another trop. To test this possibility, I wrote the script below to make an artifical allotetraploid out of positions that were homozygous for different variants in Liberia and trop.

```perl
#!/usr/bin/env perl
use strict;
use warnings;
#use lib qw(~/perl_modules);
#use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
# Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.36.1

# actually don't need to do these:
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and generates a fake allotetraploid using genotypes from two of the diploid sequences

# It also will count and compare the number of pseudoheterozygous positions that are generated in the 
# artificial allotetraploid and compare this to heterozygous counts in X. cal and X. mel

# It additionally will count the number of homozygous apomorphies in Liberia and the other diploid
# to (possibly?) permit an expectation to be calculated for the number of overlapping pseudoheterozygous
# positions

# if Liberia and trop were the two ancestors of X. calcaratus, then there should be a large number of shared pseudohet 
# sites in the artificial allotet and X. cal (more than X. mell)

# Or if Lib and trop were each ancestors of X. cal and X. mel but one was maternal for one and paternal for the other
# then the shared pseudohet sites with the artificial allotet should be high for both X. cal and X. mel


# to execute type 2025_artificial_allotetraploid.pl inputfile.tab a_x_y_z output.out 
# where a, x, y, and z are the columns that have the trop, liberia, calcaratus, and mellotrop genotypes respectively.
# these numbers shoudl start at 1 for the first sample in the file, etc.

# perl 2025_artificial_allotetraploid.pl inputfile.tab 11_18_19_20 output.out
# this uses NG5 (from Nigeria) and Liberia plus cal and mel


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

my @samples = split("_",$ARGV[1]);
my @temp;
my @cal;
my @mel;
my @trop;
my @liberia;
my $number_of_pseudohets_in_articificial_allo=0;
my $number_of_shared_pseudohets_between_artificial_and_Xcal = 0;
my $number_of_shared_pseudohets_between_artificial_and_Xmel = 0;
my $number_of_shared_pseudohets_between_artificial_and_Xcal_and_Xmel = 0;
my $number_of_shared_pseudohets_between_Xcal_and_Xmel_but_not_artificial = 0;
my $number_of_hets_in_trop=0;
my $number_of_hets_in_lib=0;
my $number_of_hets_in_cal=0;
my $number_of_hets_in_mel=0;

my @check;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/\t/,$line);
	if($temp[0] ne '#CHROM'){
		@trop = split("/",$temp[$samples[0]+2]);
		@liberia = split("/",$temp[$samples[1]+2]);
		@cal = split("/",$temp[$samples[2]+2]);
		@mel = split("/",$temp[$samples[3]+2]);
		if(($trop[0] ne '.')&&($trop[1] ne '.')&&($trop[0] ne '*')&&($trop[1] ne '*')&&
		($liberia[0] ne '.')&&($liberia[1] ne '.')&&($liberia[0] ne '*')&&($liberia[1] ne '*')&&
		($cal[0] ne '.')&&($cal[1] ne '.')&&($cal[0] ne '*')&&($cal[1] ne '*')&&
		($mel[0] ne '.')&&($mel[1] ne '.')&&($mel[0] ne '*')&&($mel[1] ne '*')
		)
		{
			# focus only on positions that are homozygous
			# for different SNPs in trop and liberia 
			if($trop[0] ne $trop[1]){
				$number_of_hets_in_trop+=1;
			}		
			if($liberia[0] ne $liberia[1]){
				$number_of_hets_in_lib+=1;
			}		
			if($cal[0] ne $cal[1]){
				$number_of_hets_in_cal+=1;
			}		
			if($mel[0] ne $mel[1]){
				$number_of_hets_in_mel+=1;
			}		
			if(
			   (($trop[0] eq $trop[1])&&($liberia[0] eq $liberia[1]))&&
		 	  (($trop[0] ne $liberia[0]))
		 	  ){   
				  # trop and liberia are homoz for different SNPs
				# keep track of these positions 
				$number_of_pseudohets_in_articificial_allo+=1;
					# check if this position is also heterozygous for the same SNPs in cal
					if(
					(($cal[0] ne $cal[1])&&($cal[0] eq $trop[0])&&($cal[1] eq $liberia[0]))
			  		||
			  		(($cal[0] ne $cal[1])&&($cal[0] eq $liberia[0])&&($cal[1] eq $trop[0]))
			  		)
			  		{
						$number_of_shared_pseudohets_between_artificial_and_Xcal+=1;
					}
					# check if this position is also heterozygous for the same SNPs in cal
					if(
			  		(($mel[0] ne $mel[1])&&($mel[0] eq $trop[0])&&($mel[1] eq $liberia[0]))
			  		||
			  		(($mel[0] ne $mel[1])&&($mel[0] eq $liberia[0])&&($mel[1] eq $trop[0]))
			  		)
			  		{
						$number_of_shared_pseudohets_between_artificial_and_Xmel+=1;
						# check if this position is also heterozygous for the the same SNPs in cal
						if(
				  		(($cal[0] ne $cal[1])&&($cal[0] eq $trop[0])&&($cal[1] eq $liberia[0]))
				  		||
				  		(($cal[0] ne $cal[1])&&($cal[0] eq $liberia[0])&&($cal[1] eq $trop[0]))
				  		)
				  		{
							$number_of_shared_pseudohets_between_artificial_and_Xcal_and_Xmel+=1;
						}
			  		}
				}
			else{
				# for positions that are not homozygous for different SNPs in Lib and trop
				# check if there is a shared pseudohet between X. cal and X. mel
				# for positions where trop and liberia are not homoz for different variants
				if(
			 	 (($mel[0] ne $mel[1])&&($cal[0] ne $cal[1])&&(($mel[0] eq $cal[0])&&($mel[1] eq $cal[1])))
			  	||
			  	(($mel[0] ne $mel[1])&&($cal[0] ne $cal[1])&&(($mel[0] eq $cal[1])&&($mel[1] eq $cal[0])))
			  	 )
			  	{
					$number_of_shared_pseudohets_between_Xcal_and_Xmel_but_not_artificial+=1;
			  	}
			}
		} # end of if to check for gaps and missing data		
	} # end of check to see if we are at the first line	
} # end while

print OUTFILE1 "#_hets_in_trop ",$number_of_hets_in_trop,"\n";
print OUTFILE1 "#_hets_in_lib ",$number_of_hets_in_lib,"\n";
print OUTFILE1 "#_hets_in_cal ",$number_of_hets_in_cal,"\n";
print OUTFILE1 "#_hets_in_mel ",$number_of_hets_in_mel,"\n";
print OUTFILE1 "#_of_pseudohets_in_articificial_allo: ",$number_of_pseudohets_in_articificial_allo,"\n";
print OUTFILE1 "#_of_shared_pseudohets_between_artificial_and_Xcal: ",$number_of_shared_pseudohets_between_artificial_and_Xcal,"\n";
print OUTFILE1 "#_of_shared_pseudohets_between_artificial_and_Xmel: ",$number_of_shared_pseudohets_between_artificial_and_Xmel,"\n";
print OUTFILE1 "#_of_shared_pseudohets_between_artificial_and_Xcal_and_Xmel: ",$number_of_shared_pseudohets_between_artificial_and_Xcal_and_Xmel,"\n";
print OUTFILE1 "#_of_shared_pseudohets_between_Xcal_and_Xmel_but_not_artificial: ",$number_of_shared_pseudohets_between_Xcal_and_Xmel_but_not_artificial,"\n";

close OUTFILE1;


print "#_hets_in_trop ",$number_of_hets_in_trop,"\n";
print "#_hets_in_lib ",$number_of_hets_in_lib,"\n";
print "#_hets_in_cal ",$number_of_hets_in_cal,"\n";
print "#_hets_in_mel ",$number_of_hets_in_mel,"\n";
print "#_of_pseudohets_in_articificial_allo: ",$number_of_pseudohets_in_articificial_allo,"\n";
print "#_of_shared_pseudohets_between_artificial_and_Xcal: ",$number_of_shared_pseudohets_between_artificial_and_Xcal,"\n";
print "#_of_shared_pseudohets_between_artificial_and_Xmel: ",$number_of_shared_pseudohets_between_artificial_and_Xmel,"\n";
print "#_of_shared_pseudohets_between_artificial_and_Xcal_and_Xmel: ",$number_of_shared_pseudohets_between_artificial_and_Xcal_and_Xmel,"\n";
print "#_of_shared_pseudohets_between_Xcal_and_Xmel_but_not_artificial: ",$number_of_shared_pseudohets_between_Xcal_and_Xmel_but_not_artificial,"\n";

```
