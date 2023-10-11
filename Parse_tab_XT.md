# Paths
mello
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/Germany_genome/raw_data
```
Liberia
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/Liberia/trim_fq
```
tads
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT7_WY_trim_noadapters
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT11_WW_trim_noadapters
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT10_WZ_trim_noadapters
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XTR_1_ZY/trimmed
```
newseqs
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/raw_data/combined/gvcfs
```


# This script looks for W-linked, Z-linked, and Y-linked SNPs in XT
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
my @male_specific_nucleotides=();
my $male_het_nuc1=0;
my $male_het_nuc2=0;
my $female_hom_nuc1=0;

my $f_diverged=0;

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

		if(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides != -1)){
			# this means that there is at least one genotype in each sex
			# we can compare homoz and het genotypes because both sexes have data
			# First check for W-linked SNPs; 
				# all males homozygous 
				# all females are heterozygous or homozygous for a variant
				# that is never homoz in a male
			if(($#unique_male_nucleotides == 0)&&($#unique_female_nucleotides > 0)){
				# this is one scenario where all males are ZZ
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
			} # end of check for scenario 1 for W-linked SNPs
			elsif(($#unique_male_nucleotides > 0)&&($#unique_female_nucleotides > 0)){
				# this is a second scenario where at least some males could be WY or ZY
				# need to make sure no males are WW
				# so if a fem is homoz, need to check if any of the males are homoz for this position as well
				$not_W=0;
				$diverged=0;
				$male_het_nuc1=0;
				$male_het_nuc2=0;
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) {
					if($females[$x] eq $females[$x+1]){
						# this female is homoz
						# need to go through each male to check if this nucleotide is also homoz in any male
						for ($y = 0 ; $y <= $#males ; $y=$y+2 ) {
							if(($males[$y] eq $males[$y+1])&&($males[$y] eq $females[$x])){
								# this is a male that was homoz for a nucleotide that was also homoz in a female
								# so it is not a W-linked position
								$not_W=1;
							}
						}
					}
					elsif($females[$x] ne $females[$x+1]){
						# this female is heteroz
						# go through each male to check if both of these nucleotides are homoz in different males
						# this could be something like this: female: AT, male1: AA, male2: TT (which is not W-linked)
						for ($y = 0 ; $y <= $#males ; $y=$y+2 ) {
							if(($males[$y] eq $males[$y+1])&&($males[$y] eq $females[$x])){
								$male_het_nuc1=1;
							}
							elsif(($males[$y] eq $males[$y+1])&&($males[$y] eq $females[$x])){
								$male_het_nuc2=1;
							}	
						}
						$diverged+=1;  # this is the number of heterozygous females
					}
				}
				# ok I went through each female and checked if there are any pairs of males that are 
				# homoz for both nucleotides
				if(($male_het_nuc1==1)&&($male_het_nuc2==1)){
					$not_W=1;
				}
				if($not_W==0){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tW_linked\t",$diverged,"\t0\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# W-linked variation
				}								
			}

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
						# homoz females, which would mean that
						# this position is not Y-linked (because
						# we can't be homoz for YY)
						if($males[$x] ne $unique_female_nucleotides[0]){
							$not_Y=1;
						}	
					}
					else{
						# this male is heteroz and one of the SNPs is never
						# homoz in females
						$diverged+=1; # this is the number of heterozygous males
					}
				}
				if($not_Y==0){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tY_linked\t0\t$diverged\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Y-linked variation
				}								
			} # end of check for Y-linked SNPs with females all homoz
			if(($#unique_male_nucleotides > 0)&&($#unique_female_nucleotides > 0)){
				# the males and females both have variation 
				# need to check if there are unique SNPs in males (e.g. females are WW and WZ and males are ZZ, WY, and ZY)
				# figure out which SNPs are uniquely male and then check if any males are homoz for these SNPs
				# but also allow for all females to be WW and all males to be ZZ and ZY
				$not_Y=0;
				$diverged=0;
				for ($x = 0 ; $x <= $#males; $x=$x+2 ) {
					if($males[$x] ne $males[$x+1]){
						# this male is heteroz, and potentially ZY
						# check if one of the variants is never found in 
						# females, which would mean that
						# this position is not Y-linked (because
						# we can't be homoz for YY)
						for ($y = 0 ; $y <= $#females ; $y=$y+2 ) {
							if(
							(($males[$x] eq $females[$y])||($males[$x] eq $females[$y+1]))
							&&
							(($males[$x+1] eq $females[$y])||($males[$x+1] eq $females[$y+1]))
							){
								# neither of these nucleotides is male specific
								$not_Y=1;
							}							
						}
						$diverged+=1; # this is the number of heterozygous males	
					}
				}
				# check for two different homoz genotypes in females and do not count this
				# position as a Y-linked one if you find this (e.g. fems: C/C, T/T and males C/C and T/C)
				$female_hom_nuc1=0;
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) {
					if($females[$x] eq $females[$x+1]){
						if($female_hom_nuc1 ne 0){
							$female_hom_nuc1=$females[$x];
						}
						elsif($female_hom_nuc1 ne $females[$x]){
							$not_Y=1;
						}
					}
				}
				$f_diverged=0;
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) {
					if($females[$x] ne $females[$x+1]){
						$f_diverged+=1;
					}
				}	
				if($not_Y==0){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tY_linked\t$f_diverged\t$diverged\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Y-linked variation
				}								
			}
			
			
							
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
# Plotting

```R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Run like this:
# Rscript --vanilla W_Y_Z_specific_overlay_computecanada.R Nigeria_only_221122220022222222.out
library (ggplot2)
library(reshape2) # this facilitates the overlay plot
setwd("./")


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please specify the (input file).n", call.=FALSE)
}

# read in the data and name the df based on the file name
locations <- read.table(args[1], header = T, sep="\t")

subset_WY_specific <- locations[locations$TYPE!="Z_linked",]


png(filename = paste(args[1],"_W_Y_Z_specific.png",sep=""),w=1200, h=800,units = "px", bg="transparent")
  ggplot(locations, aes(x=POS/1000000, fill=TYPE)) +
  geom_histogram(alpha=0.2, position="identity", binwidth = 0.5) + 
  xlab("Position (Mb)") + ylab("Density") +
  theme_classic(base_size = 22)
dev.off()

png(filename = paste(args[1],"_W_Y_specific.png",sep=""),w=1200, h=800, units = "px", bg="transparent")
  ggplot(subset_WY_specific, aes(x=POS/1000000, fill=TYPE)) +
  geom_histogram(alpha=0.2, position="identity", binwidth = 0.5) + 
  xlab("Position (Mb)") + ylab("Density") +
  theme_classic(base_size = 22)
dev.off()


# make a density plot of #het_fems/#fems and #het_males/#males

subset_W_homoz <- locations[(locations$het_females/locations$n_FEMs==0) & (locations$TYPE == "W_linked"),]
subset_Z_homoz <- locations[(locations$het_males/locations$n_MALS==0) & (locations$TYPE == "Z_linked"),]

combined <- rbind(subset_W_homoz,subset_Z_homoz)
png(filename = paste(args[1],"_hetfems_hetmale.png",sep=""),w=1200, h=800, units = "px", bg="transparent")
  ggplot(combined, aes(x=POS/1000000, fill=TYPE)) +
  geom_histogram(alpha=0.2, position="identity", binwidth = 0.5) + 
  xlab("Position (Mb)") + ylab("Density") +
  scale_fill_manual('Legend Name', labels=c('W-linked fem hets', 'Z-linked male hets'), values=c('pink', 'green')) +
  theme_classic(base_size = 22) 
dev.off()

```

# Looking at interesting regions

In this directory (which doesn't have the new calcaratus seq):
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_before_filtering/
```

first subset a vcf
```
module load bcftools
bcftools view -R Chr7_9577608_9593589.bed -o Chr7_9577608_9593589.vcf combined_Chr7.g.vcf.gz_Chr7_GenotypedSNPs.vcf.gz_filtered.vcf.gz
```
or (somewhere else):
```
bcftools view -R Chr7_9550000_9600000.bed -o Chr7_9550000_9600000.vcf combined_Chr7.g.vcf.gz_Chr7_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz
```
now convert to tab
```
module load vcftools
zcat file.vcf.gz | vcf-to-tab > out.tab
```

Parsetab can be run using this sbatch script:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/ben_scripts/2022_Parse_tab.sh
```
