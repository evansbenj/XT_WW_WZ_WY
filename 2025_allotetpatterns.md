# Path
Filtered tab files (and gvcfs) are here:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_after_filtering
```

# AB'BB, AB'AB, ABB'B, and AAB'B sites
Here is a script to count these different types of sites

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
#  by vcftools (vcf2tab) and searches for sites that have an ABBA and BABA site pattern
# where the taxa are in this order: trop, liberia, mello1, mello2
# the genotypes of mello1 and mello2 is based on a heterozygous call when mello data are mapped to trop

# if trop is more closely related to mello1 than to liberia, we expect more BABA sites than ABBA
# if liberia is more closely related to mello1 than to trop, we expect more ABBA sites than BABA
# if trop and liberia are more closely related to each other than to mello1, then we expect and equal frequency of ABBA and BABA


# to execute type Liberia_ABBA_BABA.pl inputfile.tab a_x_y_z output.out 
# where a, x, y, and z are the columns that have the trop, liberia, calcaratus, and mellotrop genotypes respectively.
# these numbers shoudl start at 1 for the first sample in the file, etc.

# perl Liberia_ABBA_BABA.pl inputfile.tab 11_18_19_20 output.out

# There are 4 site patterns: 
#         AB'BB and AB'AB (this is with the first allotet as the focal)
# 	      and ABB'B and AAB'B (this is with the second allotet as the focal)
#		  the first A is the first diploid genotype; the last B is the second diploid genotype
# 		  in these patterns, B' is a heterozygous genotype; A and B are homozygous genotypes


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
my @mello;
my @trop;
my @liberia;
my @all_of_them;
# AB'BB and AB'AB
my $ABpBB=0;
my $ABpAB=0;
# ABB'B and AAB'B
my $ABBpB=0;
my $AABpB=0;
my @check;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/\t/,$line);
	if($temp[0] ne '#CHROM'){
		@cal = split("/",$temp[$samples[2]+2]);
		@mello = split("/",$temp[$samples[3]+2]);
		if(
		   (($cal[0] ne $cal[1])||($mello[0] ne $mello[1]))&&
		   (($cal[0] eq $cal[1])||($mello[0] eq $mello[1]))
		   ){ 
			# either cal or mel genotype is heteroz, but not both
			@trop = split("/",$temp[$samples[0]+2]);
			@liberia = split("/",$temp[$samples[1]+2]);
			# first make sure both diploids are homozygous
			if(($trop[0] eq $trop[1])&&($liberia[0] eq $liberia[1])&&($trop[0] ne $liberia[0])){
				# now check if there are more than one variant
				# both trop samples are homozygous for different variants
				@all_of_them = (@trop, @liberia, @cal, @mello);
				@check = grep /\*/, @all_of_them;
				# print "hello $#check\n";
				if((uniq(@all_of_them) eq 2)&&($#check == -1)){
				# there are only two variants	
					# count sites
					# first count AB'BB and AB'AB  - these have the first allotet heteroz
					if($cal[0] ne $cal[1]){
						if($mello[0] eq $liberia[0]){
							# the second allotet homozygous genotype is the same as the second diploid
							$ABpBB+=1;
							print "ABpBB ",@trop," ",@cal," ",@mello," ",@liberia,"\n";
						}
						elsif($mello[0] eq $trop[0]){
							# the second allotet homozygous genotype is the same as the first diploid
							$ABpAB+=1;
							print "ABpAB ",@trop," ",@cal," ",@mello," ",@liberia,"\n";
						}	
					}
					# now count ABB'B and AAB'B  - these have the second allotet heteroz
					elsif($mello[0] ne $mello[1]){
						if($cal[0] eq $liberia[0]){
							# the first allotet homozygous genotype is the same as the second diploid
							$ABBpB+=1;
							print "ABBpB ",@trop," ",@cal," ",@mello," ",@liberia,"\n";
						}
						elsif($cal[0] eq $trop[0]){
							# the first allotet homozygous genotype is the same as the first diploid
							$AABpB+=1;
							print "AABpB ",@trop," ",@cal," ",@mello," ",@liberia,"\n";
						}	
					}
				}				
			}
		}
	} # end of check to see if we are at the first line	
	else{
		print OUTFILE $temp[$samples[0]+2],"\t",$temp[$samples[1]+2],"\t",$temp[$samples[2]+2],"\t",$temp[$samples[3]+2],"\n";
		print $temp[$samples[0]+2],"\t",$temp[$samples[1]+2],"\t",$temp[$samples[2]+2],"\t",$temp[$samples[3]+2],"\n";
	}
} # end while

print "Number of AB'BB sites: ",$ABpBB,"\n";
print "Number of AB'AB sites: ",$ABpAB,"\n";
print "Number of ABB'B sites: ",$ABBpB,"\n";
print "Number of AAB'B sites: ",$AABpB,"\n";

print OUTFILE1 "Number of AB'BB sites: ",$ABpBB,"\n";
print OUTFILE1 "Number of AB'AB sites: ",$ABpAB,"\n";
print OUTFILE1 "Number of ABB'B sites: ",$ABBpB,"\n";
print OUTFILE1 "Number of AAB'B sites: ",$AABpB,"\n";

close OUTFILE1;


sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
```

# Apomorphy script

I also developed an apomorphy script

``` perl
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
# Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and searches for sites that have an ABB and BAB site pattern
# where the taxa are in this order: trop, liberia, allotet
# all genotypes must be homozygous
# the genotypes of allotet is based on a homozygous call when mello data are mapped to trop

# if trop is more closely related to allotet than to liberia, we expect more BAB sites than ABB
# if liberia is more closely related to allotet than to trop, we expect more ABB sites than BAB
# if trop and liberia are more closely related to each other than to allotet, then we expect and equal frequency of ABB and BAB

# to execute type Liberia_ABBA_BABA.pl inputfile.tab x_y_z output.out 
# where x, y, and z are the columns that have the trop, liberia, and allotet genotypes respectively.
# these numbers should start at 1 for the first sample in the file, etc.

# perl Liberia_ABBA_BABA.pl inputfile.tab 11_18_19 output.out


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
#print OUTFILE1 "CHR\tPOS\tTYPE\n";

my @samples = split("_",$ARGV[1]);
my @temp;
my @allotet;
my @trop;
my @liberia;
my $ABB=0;
my $BAB=0;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/\t/,$line);
	if($temp[0] ne '#CHROM'){
		my @allotet = split("/",$temp[$samples[2]+2]);
		if(($allotet[0] eq $allotet[1])&&($allotet[0] ne '.')&&($allotet[0] ne '*')){ # mello genotype is homoz
			my @trop = split("/",$temp[$samples[0]+2]);
			my @liberia = split("/",$temp[$samples[1]+2]);
			# for now focus on homoz genotypes in trop and liberia
			if(($trop[0] eq $trop[1])&&($liberia[0] eq $liberia[1])&&($trop[0] ne '.')&&
			($liberia[0] ne '.')&&($trop[0] ne '*')&&($liberia[0] ne '*')){
				if(($trop[0] ne $liberia[0])&&($trop[0] eq $allotet[0])){
					# this is a BAB site
					$BAB+=1;
					print "BAB $temp[0] $temp[1] @trop @liberia @allotet\n";
				}
				elsif(($trop[0] ne $liberia[0])&&($liberia[0] eq $allotet[0])){
					# this is a ABB site
					$ABB+=1;
					print "ABB $temp[0] $temp[1]  @trop @liberia @allotet\n";
				}
				
			}
		}
	} # end of check to see if we are at the first line	
	else{
		print $temp[$samples[0]+2],"\t",print $temp[$samples[1]+2],"\t",print $temp[$samples[2]+2],"\t",
	}
} # end while
print OUTFILE1 "Number of trop1,trop2,Allotet BAB sites: ",$BAB,"\n";
print OUTFILE1 "Number of trop1,trop2,Allotet ABB sites: ",$ABB,"\n";
print "Number of trop1,trop2,Allotet BAB sites: ",$BAB,"\n";
print "Number of trop1,trop2,Allotet ABB sites: ",$ABB,"\n";
close OUTFILE1;

```

# ABB_BAB_with_windows
```
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
# Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and searches for sites that have an ABB and BAB site pattern
# where the taxa are in this order: trop, liberia, allotet
# all genotypes must be homozygous
# the genotypes of allotet is based on a homozygous call when mello data are mapped to trop

# if trop is more closely related to allotet than to liberia, we expect more BAB sites than ABB
# if liberia is more closely related to allotet than to trop, we expect more ABB sites than BAB
# if trop and liberia are more closely related to each other than to allotet, then we expect and equal frequency of ABB and BAB

# to execute type Liberia_ABBA_BABA.pl inputfile.tab x_y_z output.out 
# where x, y, and z are the columns that have the trop, liberia, and allotet genotypes respectively.
# these numbers should start at 1 for the first sample in the file, etc.

# perl Liberia_ABBA_BABA.pl inputfile.tab 11_18_19 output.out

# EUA335 11
# AMNH17271 12
# lib: 18
# cal: 19
# mel: 20

# ./2025_Liberia_autapomorphy_windowz.pl all.tab 11_18_20 EUA335_lib_mel
# ./2025_Liberia_autapomorphy_windowz.pl all.tab 11_18_19 EUA335_lib_cal

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
#print OUTFILE1 "CHR\tPOS\tTYPE\n";

unless (open(OUTFILE2, ">".$outputfile1."_windowz"))  {
	print "I can\'t write to ".$outputfile1."_windowz\n";
	exit;
}
print "Creating output file: ".$outputfile1."_windowz\n";
#print OUTFILE1 "CHR\tPOS\tTYPE\n";

my @samples = split("_",$ARGV[1]);
my @temp;
my @allotet;
my @trop;
my @liberia;
my $ABB=0;
my $BAB=0;
my $ABB_window=0;
my $BAB_window=0;
my $current_chr="Chr10"; # all.tab starts with Chr10
my $previous_chr="Chr10";
my $window=5000000; # use 5Mb windows
my $switch=0;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/\t/,$line);
	if($temp[0] ne '#CHROM'){
		$current_chr=$temp[0];
		my @allotet = split("/",$temp[$samples[2]+2]);
		if(($allotet[0] eq $allotet[1])&&($allotet[0] ne '.')&&($allotet[0] ne '*')){ # mello genotype is homoz
			my @trop = split("/",$temp[$samples[0]+2]);
			my @liberia = split("/",$temp[$samples[1]+2]);
			# for now focus on homoz genotypes in trop and liberia
			if(($trop[0] eq $trop[1])&&($liberia[0] eq $liberia[1])&&($trop[0] ne '.')&&
			($liberia[0] ne '.')&&($trop[0] ne '*')&&($liberia[0] ne '*')){
				if(($trop[0] ne $liberia[0])&&($trop[0] eq $allotet[0])){
					# this is a BAB site
					$BAB+=1;
					$BAB_window+=1;
					# print "BAB $temp[0] $temp[1] @trop @liberia @allotet\n";
				}
				elsif(($trop[0] ne $liberia[0])&&($liberia[0] eq $allotet[0])){
					# this is a ABB site
					$ABB+=1;
					$ABB_window+=1;
					# print "ABB $temp[0] $temp[1]  @trop @liberia @allotet\n";
				}
			}
			if(($temp[1] > $window)||($current_chr ne $previous_chr)){
				print OUTFILE2 $previous_chr,"\t",$window,"\t",$ABB_window,"\t",$BAB_window,"\n";
				print $previous_chr,"\t",$window,"\t",$ABB_window,"\t",$BAB_window,"\n";
				if($current_chr ne $previous_chr){
					$window = 5000000;
				}
				else{
					$window += 5000000;
				}	
				$ABB_window=0;
				$BAB_window=0;
				$previous_chr = $current_chr;
			}
		}
	} # end of check to see if we are at the first line	
	else{
		if($switch == 0){
			print $temp[$samples[0]+2],"\t",$temp[$samples[1]+2],"\t",$temp[$samples[2]+2],"\n";
			print OUTFILE2 "Chr\twindow\tABB\tBAB\n";
		}
		$switch=1;
	}
} # end while

# print the last line
print OUTFILE2 $previous_chr,"\t",$window,"\t",$ABB_window,"\t",$BAB_window,"\n";
print $previous_chr,"\t",$window,"\t",$ABB_window,"\t",$BAB_window,"\n";



print OUTFILE1 "Number of trop1,trop2,Allotet BAB sites: ",$BAB,"\n";
print OUTFILE1 "Number of trop1,trop2,Allotet ABB sites: ",$ABB,"\n";
print "Number of trop1,trop2,Allotet BAB sites: ",$BAB,"\n";
print "Number of trop1,trop2,Allotet ABB sites: ",$ABB,"\n";
close OUTFILE1;
close OUTFILE2;
```
# Block Bootstrap for ABB_BAB
```
#!/usr/bin/env perl
use strict;
use warnings;

use Statistics::Descriptive;

# This program reads in the output of a script for Libera ABB and BAB sites
# and calculates the standard error of the weighted mean value of with weightings based 
# on the sum of the number of ABB and BAB sites in each window.

# From SI of Green et al. 2010 Neanderthal genome paper:
# "By computing the variance of the statistic over the entire genome M times leaving each
# block of the genome in turn, and then multiplying by M and taking the square root, we can obtain an
# approximately normally distributed standard error using the theory of the jackknife. "

# I can do this for ABB and BAB sites

my $inputfile = $ARGV[0];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my @temp; 
my @ABB=();
my @BAB=();
my $counter=0;
my $counter2=0;
my $ABB=0;
my $BAB=0;
my $jack_ABB;
my $jack_BAB;
my $jack_ABB_boot;
my $jack_BAB_boot;
my @jackABBarray;
my @jackBABarray;
my @jack_ABB;
my @jack_BAB;
my $y;
my $x;


while ( my $line = <DATAINPUT>) {
	@temp=split('\t',$line);
	if($temp[0] ne 'Chr'){
		if(($temp[2] !~ /nan/)&&($temp[3] !~ /nan/)){
			push(@ABB,$temp[2]);
			push(@BAB,$temp[3]);
			$counter+=1;
		}
	}	
}		

close DATAINPUT;

# check that array lengths are ok
if(($#ABB ne $#BAB)||($#ABB ne ($counter-1))){
	print "Problem!\n";
}

for ($y = 0 ; $y < $counter; $y++ ) {
	$jack_ABB+=$ABB[$y]; # this is the global stat
	$jack_BAB+=$BAB[$y]; # this is the global stat
	for ($x = 0 ; $x < $counter; $x++ ) {
		# leave out one row for each jackknfe replicates
		if($y != $x){
			# load arrays with values for each window for later variance calculations
			$jack_ABB_boot += $ABB[$x];
			$jack_BAB_boot += $BAB[$x];
			# keep track of the number of windows
			$counter2+=1;			
		}
	}
	# store the stats from each bootstrap replicate
	push(@jackABBarray,$jack_ABB_boot);
	push(@jackBABarray,$jack_BAB_boot);
	# reset variables
	$jack_ABB_boot=0;
	$jack_BAB_boot=0;
	$counter2=0;
}


# check that array lengths are ok
if($#jackABBarray ne ($counter-1)){
	print "Problem2!\n";
}

# now calculate the variance of the jackknife replicates
my $jack_ABBmean=0;
my $jack_BABmean=0;

# first we need the mean
for ($x = 0 ; $x < $counter; $x++ ) {
	$jack_ABBmean+=$jackABBarray[$x];
	$jack_BABmean+=$jackBABarray[$x];
}
$jack_ABBmean=$jack_ABBmean/($counter);
$jack_BABmean=$jack_BABmean/($counter);

# now get the variance
my $jack_ABBvar=0;
my $jack_BABvar=0;
for ($x = 0 ; $x < $counter; $x++ ) {
	$jack_ABBvar+=($jack_ABBmean-$jackABBarray[$x])**2;
	$jack_BABvar+=($jack_BABmean-$jackBABarray[$x])**2;
}
# for the sample variance, divide by (n-1)
$jack_ABBvar=$jack_ABBvar/($counter-1);
$jack_BABvar=$jack_BABvar/($counter-1);

# get the standard error
my $sterr_ABB = sqrt($counter*$jack_ABBvar);
my $sterr_BAB = sqrt($counter*$jack_BABvar);

# print the CIs
print "The number of ABB sites is ",sprintf("%.0f",$jack_ABB)," (",sprintf("%.0f",($jack_ABB-1.96*$sterr_ABB))," - ",sprintf("%.0f",($jack_ABB+1.96*$sterr_ABB)),")\n";
print "The number of BAB sites is ",sprintf("%.0f",$jack_BAB)," (",sprintf("%.0f",($jack_BAB-1.96*$sterr_BAB))," - ",sprintf("%.0f",($jack_BAB+1.96*$sterr_BAB)),")\n";

print "The standard error of ABB is ",sprintf("%.0f",$sterr_ABB),"\n";
print "The standard error of BAB is ",sprintf("%.0f",$sterr_BAB),"\n";
print "The 95\%CI of the weighted ABB is ",
sprintf("%.0f",($jack_ABB-1.96*$sterr_ABB))," - ",sprintf("%.0f",($jack_ABB+1.96*$sterr_ABB)),"\n";
print "The 95\%CI of the weighted BAB is ",
sprintf("%.0f",($jack_BAB-1.96*$sterr_BAB))," - ",sprintf("%.0f",($jack_BAB+1.96*$sterr_BAB)),"\n";

```
