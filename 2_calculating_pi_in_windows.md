# Calculating pi in windows for each genome

I'm going to try three approaches - vcftools, tabfile + myscript, and general_genomics. The first works from a vcf file, which I now have.  The last works from a geno file. 

First phase with Beagle:
```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=def-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_
depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/beagle.18May20.d20.jar gt=${1} out=${1}_phased.vcf
.gz impute=true
```


Then make the geno file like this (2020_make_geno_from_vcf.sh):
```
#!/bin/sh
#SBATCH --job-name=makegeno
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=2gb
#SBATCH --output=makegeno.%J.out
#SBATCH --error=makegeno.%J.err
#SBATCH --account=def-ben

# sbatch 2020_make_geno_from_vcf.sh path_and_name_of_vcf.gz_file

python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dep
th_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i ${1} 
| gzip > ${1}.geno.gz
```

then it is necessary to swap asterisks with Ns:
```
gunzip ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz
sed -i 's/\*/N/g' ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno
gzip -c ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno > ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz
```

Here is an sbatch script that runs the popgenWindows calculation of Fst:
```
#!/bin/sh
#SBATCH --job-name=popgen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=2gb
#SBATCH --output=popgen.%J.out
#SBATCH --error=popgen.%J.err
#SBATCH --account=def-ben


# sbatch ./popgenWindows_allautosomes.sh pop1 pop2

# populations
# bru papio hec mau nem sum nig nge tog ton

module --force purge
module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2

## declare an array variable
declare -a chrs=("chr01" "chr02a" "chr02b" "chr03" "chr04" "chr05" "chr06" "chr07" "chr08" "chr09" "chr10" "chr11
" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19")

for file in "${chrs[@]}"; do  

    echo python3 ./popgenWindows.py -g ./VCF_processing/${file}.geno.gz -o ./VCF_processing/popgenWindows_${file}
_${1}_${2}_windowstats.csv -w 100000 -m 100 -s 100000 -p ${1} -p ${2} -f phased -T 10 --popsFile pops.txt --write
FailedWindows --windType coordinate

python3 ./popgenWindows.py -g ./VCF_processing/${file}.geno.gz -o ./VCF_processing/popgenWindows_${file}_${1}_${2
}_windowstats.csv -w 100000 -m 100 -s 100000 -p ${1} -p ${2} -f phased -T 10 --popsFile pops.txt --writeFailedWin
dows --windType coordinate

done
```


The sbatch scripts are not working for some reason probably related to module incompatibility.  It runs really quickly directly though:

```
./2020_popgenWindows.sh ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz XT7_WY XT10_WZ

python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/popgenWindows.py -g ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz -o ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz_XT7_WY_XT10_WZ.csv -m 1 -p XT7_WY -p XT10_WZ -f phased -T 10 --popsFile pops.txt --writeFailedWindows -w 10000 -s 10000 -m 10 --windType coordinate

./2020_popgenWindows.sh ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz XT10_WZ XT11_WW

python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/popgenWindows.py -g ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz -o ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz_XT10_WZ_XT11_WW.csv -m 1 -p XT10_WZ -p XT11_WW -f phased -T 10 --popsFile pops.txt --writeFailedWindows -w 10000 -s 10000 -m 10 --windType coordinate

./2020_popgenWindows.sh ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz XT11_WW XT7_WY

python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/popgenWindows.py -g ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz -o ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz_XT11_WW_XT7_WY.csv -m 1 -p XT11_WW -p XT7_WY -f phased -T 10 --popsFile pops.txt --writeFailedWindows -w 10000 -s 10000 -m 10 --windType coordinate
```

Update Nov 2021 - I now have 4 individuals sequenced and genotyped.
example commandline:
```
./2021_general_genomics_popgen_4pops.sh ../combined_and_genotyped_vcf_filez_SNPsonly/MandF_Chr7.g.vcf.gz_Chr7_SNPs_pos1_30Mb.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz XT10_WZ XT11_WW XT1_ZY XT7_WY
```
with 2021_general_genomics_popgen_4pops.sh being this:

```
#!/bin/sh
#SBATCH --job-name=popgen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=2gb
#SBATCH --output=popgen.%J.out
#SBATCH --error=popgen.%J.err
#SBATCH --account=def-ben


# sbatch ./2021_general_genomics_popgen_2pops.sh genofile pop1 pop2

module --force purge
module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2


python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_fi
ltered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/popgenWind
ows.py -g ${1} -o ${2}_${3}_${4}_${5}_windowstats.csv -w 10000 -m 100 -s 10000 -p ${2} -p ${3} -p ${4} -
p ${5} -f phased -T 10 --popsFile pops.txt --writeFailedWindows --windType coordinate
```




The tab file approach is made with vcftools and my script is this (Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts_.pl):
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use Sort::Naturally;


#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and calculates polymorphism stats in windows

#	This analysis will allow for missing data but in doing so assumes that the missing data are randomly distributed.

# 	If this were not the case, for example if a diverged sample had lots of missing data, 
#	then the estimates would be biased (downward biased in that example)

# 	This program will be compatible with files with one or multiple outgroup columns

# 	It will also accomodate data from multiple species and calculate the stats only from selected columns

# to execute type Boot_from_tab_diverge_poly_2018_allowmissingdata.pl inputfile.tab 1111100110000111100011100110010100000000 
# 3_4_22_23_25_26 nigra_poly_and_diverge.txt  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_4_22_23_25_26 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case), (ii) the column number of the first individual in the ingroup 
# (4 in this case), and (iii) the sample number that contain the data from the individuals you want to 
# include (22, 23, 25, and 26 in this case), which are the four nigra samples itemized below.

# IMPORTANT: (i) and (ii) are columns beginning with 1 but (iii) is based on the individual samples such
# as enumerated below

# for example, with a tab file with only the baboon sequence in the 4th column, here is the input command:

# tonk
# Boot_from_tab_diverge_poly_2018_allowmissingdata.pl temp.tab 1100000100011010101001111110110101001011 3_4_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39_40 temp.poly_by_windows

my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile2 = $ARGV[3];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}


unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2\n";
	exit;
}
print "Creating output file: $outputfile2\n";


my @sexes = split("",$ARGV[1]);
my @whotoinclude = split("_",$ARGV[2]);


my @temp;
my @temp1;
my $previous= 0;
my $aDNA_total_number_of_sites=0; # this is the total number of sites with at least two alleles genotyped.

# variables for sliding window
my $sliding_window=5000000;
my $current_window=0;
my $window_counter=0;
my $current_chromosome="blah";

my $pi_aDNA=0;
my %pi_aDNA_by_window=(); # this is for the pi in sliding windows
my %aDNA_sites_by_window=(); # this is the number of sites that are genotyped in each window

my $x_uniq=0;
my $diff=0;

# variables for bootstrapping
my @pi_aDNA=();
my $aDNA_divergence=0;
my $string;
my $m;
my $n;
my $w;
my $y;
my $x;
my @unique;

my $number_of_individuals_genotyped=($#whotoinclude - 1);

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped=0;
for ($y = 2 ; $y <= $#whotoinclude ; $y++ ) {
	if($sexes[$whotoinclude[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	

print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";

my $pi_counter=0;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		# base the genomic location on the outgroup	
		if($temp[0] ne $current_chromosome){
			$current_chromosome = $temp[0];
			$current_window = 0;
			$window_counter+=1;
		}
		until($temp[1] < ($current_window+$sliding_window)){ # this should never be invoked for transcriptomes
			$current_window = $current_window+$sliding_window;
			$window_counter+=1;
			print "Something wrong with sliding window\n";
		}
		# load the alleles
		$string=();
		for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
			# load the first allele
			# This setup will not count alleles that are heterozygous for a weird base
			# This aims to minimize biases where weird genotypes are called as ref
			if(
				((uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'A')||
				(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'T')||
				(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'C')||
				(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'G'))
				&&
				((uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'A')||
				(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'T')||
				(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'C')||
				(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'G'))						
				){
				# load the first allele
				$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
				$string=$string.$w;
				# now load the second allele
				$w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
				$string=$string.$w;
			}
			# if this fails, there could be genotypes at individuals that are not included in the analysis
		}
		if(defined($string)){ 
			if(length($string)>1){ # we need at least two alleles to estimate diversity
				@temp1=split('',$string);
				if(
					(uc $temp[$whotoinclude[0]-1] eq "A")||
					(uc $temp[$whotoinclude[0]-1] eq "C")||
					(uc $temp[$whotoinclude[0]-1] eq "T")||
					(uc $temp[$whotoinclude[0]-1] eq "G")
					){ # the outgroup must be defined
						$aDNA_total_number_of_sites+=1; # this is the total number of sites with at least two alleles genotyped
														# and the outgroup defined
						$aDNA_sites_by_window{$current_chromosome}+=1; # this is the count of sites genotyped per sliding windown
						$x_uniq = uniq @temp1; 			# this is the number of unique variants at this position
						if($x_uniq > 1){
							# calculate pi for this site
							$diff=0; # this is for pairwise comparisons within a site
							for ($y = 0 ; $y < $#temp1 ; $y++ ) {
								for ($x = ($y+1) ; $x <= $#temp1 ; $x++ ) {
									if($temp1[$y] ne $temp1[$x]){
										$diff+=1;
									}
									$pi_counter+=1;
								}
							}
							$pi_aDNA+=$diff/$pi_counter; 	# this will be standardized later by dividing by the total number of sites
															# $aDNA_total_number_of_sites
							$pi_aDNA_by_window{$current_chromosome}+=$diff/$pi_counter; # this is the pi per window; 
																											# will standardize by number of sites later
							$pi_counter=0;	
						}
						else{
							# there is no diversity at this position
							$pi_aDNA+=0;
							$pi_aDNA_by_window{$current_chromosome}+=0;
							$pi_counter=0;	
						}
				}
			}
			else{
				print "Something weird happened - only one allele:\n",$line,"\n";
			}
		}
		# if this fails, there could be genotypes in individuals that are not included
	}		
	elsif($temp[0] eq '#CHROM'){
		for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
			print "Individual ",$temp[$whotoinclude[$y+2]+$whotoinclude[1]-2]," is a ";
			if($sexes[$whotoinclude[$y+2]-1] == 1){
				print "female\n";
			} 
			elsif($sexes[$whotoinclude[$y+2]-1] == 0){
				print "male\n";
			} 
		}
	}	
} # end while


print OUTFILE2 "Chr\tNum_sites\tpi_window\n";



# Now standardize and print the pi value within windows
foreach my $key (nsort ((keys %pi_aDNA_by_window))){
	print OUTFILE2 $key,"\t",$aDNA_sites_by_window{$key},"\t";
	if($aDNA_sites_by_window{$key}>0){
		print OUTFILE2 sprintf("%.5f",$pi_aDNA_by_window{$key}/$aDNA_sites_by_window{$key}),"\n"; # this is pi per site
	}
	else{
		print OUTFILE2 "NAN\n";
	}	
}	


close DATAINPUT;
close OUTFILE2;
```
