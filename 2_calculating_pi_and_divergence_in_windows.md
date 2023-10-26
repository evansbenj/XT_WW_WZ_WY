# Calculating pi in windows for each genome
I'm using angsd for this.

bam files are here:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/raw_data/combined
```
First extract a new bam file that has only chr7 for each sample:
```
sbatch /home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/ben_scripts/2023_samtools_subset_bam.sh  /home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT7_WY_trim_noadapters/XT7_WY_no_adapt__sorted.bam_rg.bam XT7_WY
```

# Try using exons only to reduce noise from mismapped reads and repettive regions
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/ben_scripts/2023_samtools_subset_bam_using_bed.sh
```

# make a file (mellobam_path.txt) that has the path to the bam file (do this for each one).
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/ben_scripts/2023_angsd.sh
```
```
#!/bin/sh
#SBATCH --job-name=angsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=2gb
#SBATCH --output=angsd.%J.out
#SBATCH --error=angsd.%J.err
#SBATCH --account=def-ben


module load StdEnv/2020 angsd/0.939

#angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta -G
L 1 -out ${1}_angsd_out

angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta -GL
 1 -out ${1}_out
realSFS ${1}_out.saf.idx -P 24 -fold 1 > ${1}_out.sfs
realSFS saf2theta ${1}_out.saf.idx -sfs ${1}_out.sfs -outname ${1}_out
thetaStat do_stat ${1}_out.thetas.idx -win 50000 -step 10000  -outnames ${1}_theta.thetasWindow.gz
```
This generated an output file with the suffix ".pestPG" that has the windowed stats.

I need to divide these by the number of sites to get the pi per site...

# Divergence using general_genomics
Directory (on graham):
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_before_filtering
```

First phase with Beagle using this file `2023_Beagle_phasing.sh`:
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/ben_scripts/2023_Beagle_phasing.sh
```
```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=256gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=rrg-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by
_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/beagle.22Jul22.46e.jar gt=${1} out=${1}_phased.v
cf.gz impute=true 
```


Then make the geno file like this (2020_make_geno_from_vcf.sh):
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/ben_scripts/2020_make_geno_from_vcf.sh
```
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

Here is an sbatch script that runs the popgenWindows calculation:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_after_filtering/2023_general_genomics_popgen_14pops.sh
```
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


# sbatch ./2023_general_genomics_popgen_14pops.sh genofile pop1 pop2 pop3 ... pop14
# sbatch ./2023_general_genomics_popgen_14pops.sh genofile GE_F1 IC_F1 NI_F1 NI_F2 SL_F1 SL_F2 REF_F1 GE_M1 GW_M1 NI_M1 NI_
M2 SL_M1 SL_M2 LI_F1

module --force purge
module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2
#module load StdEnv/2023 python/3.11.5

#python3 /home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_after_filtering/genomics_general/popgenWindows.py -g ${1} -o $
{1}_windowstats.csv -w 10000 -m 100 -s 10000 -p ${2} -p ${3} -p ${4} -p ${5} -p ${6} -p ${7} -p ${8} -p ${9} -p ${10} -p ${
11} -p ${12} -p ${13} -p ${14} -p ${15} -f phased -T 10 --popsFile XT_pops.txt --writeFailedWindows --windType coordinate

python3 /home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_after_filtering/genomics_general/popgenWindows.py -g ${1} -o ${
1}_windowstats.csv -w 10000 -m 100 -s 10000 -p GE_F1 F_Ghana_WZ_BJE4687_combined__sorted.bam -p IC_F1 F_IvoryCoast_xen228_c
ombined__sorted.bam -p NI_F1 F_Nigeria_EUA0331_combined__sorted.bam -p NI_F2 F_Nigeria_EUA0333_combined__sorted.bam -p SL_F
1 F_SierraLeone_AMNH17272_combined__sorted.bam -p SL_F2 F_SierraLeone_AMNH17274_combined__sorted.bam -p REF_F1 JBL052_conca
tscafs_sorted.bam -p GE_M1 M_Ghana_WY_BJE4362_combined__sorted.bam -p GW_M1 M_Ghana_ZY_BJE4360_combined__sorted.bam -p NI_M
1 M_Nigeria_EUA0334_combined__sorted.bam -p NI_M2 M_Nigeria_EUA0335_combined__sorted.bam -p SL_M1 M_SierraLeone_AMNH17271_c
ombined__sorted.bam -p SL_M2 M_SierraLeone_AMNH17273_combined__sorted.bam -p LI_F1 all_ROM19161_sorted.bam -f phased -T 10 
--writeFailedWindows --windType coordinate
```
**** Not used below ****

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
