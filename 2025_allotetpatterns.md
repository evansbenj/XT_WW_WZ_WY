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
