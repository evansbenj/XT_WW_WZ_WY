# Liberia and others
working directory:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/raw_data/combined
```
and
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY
```
including these four tad genomes:
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT7_WY_trim_noadapters
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT11_WW_trim_noadapters
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT10_WZ_trim_noadapters
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XTR_1_ZY/trimmed
```
and mellotrop from Germany
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/Germany_genome/raw_data
```
and XT from ref genome
```
/home/ben/projects/rrg-ben/ben/2020_XT_v10_raw_data/XT-v10_rawdata/
```

In total we have 20 WGS samples:
```
Fieldnumber	Museum_or_barcode	species	Country	Locality	Sex	Sex chromosome genotype if known
AMNH17271	2644.055	Xenopus tropicalis	Sierra Leone	near Freetown	Male	
AMNH17272	unknown	Xenopus tropicalis	Sierra Leone	near Freetown	Female	
AMNH17273	unknown	Xenopus tropicalis	Sierra Leone	near Freetown	Male	
AMNH17274	unknown	Xenopus tropicalis	Sierra Leone	near Freetown	Female	
xen228	no specimen	Xenopus tropicalis	Ivory Coast	Adiopodoume	Female	
BJE4362	3D6.1D5980178C	Xenopus tropicalis	Ghana	Ametozofe (Ghana East)	Male	WY
BJE4687	3D6.1D5980179E	Xenopus tropicalis	Ghana	Ametozofe (Ghana East)	Female	WZ
BJE4360	alive	Xenopus tropicalis	Ghana	Ankasa West	male	ZY
EUA0331	Z23682	Xenopus tropicalis	Nigeria	Benin	Female	
EUA0333	Z23684	Xenopus tropicalis	Nigeria	Benin	Female	
EUA0334	Z23685	Xenopus tropicalis	Nigeria	Benin	Male	
EUA0335	Z23686	Xenopus tropicalis	Nigeria	Benin	Male	
ROM19161	ROM19161	Xenopus tropicalis	Liberia		Female	
XTR_1_ZY	tadpole	Xenopus tropicalis	Ghana		probM	ZY
XT7_WY	tadpole	Xenopus tropicalis	Ghana		probM	WY
XT11_WW	tadpole	Xenopus tropicalis	Ghana		probF	WW
XT10_WZ	tadpole	Xenopus tropicalis	Ghana		probF	WZ
JBL052 adult Xenopus tropicalis lab F probWZ
mello_GermSeq_sorted.bam_rg_rh	BJE3652	X. mellotropicalis
calcaratus_3D6.1D5980.1953	3D6.1D5980.1953	X. calcaratus
```

# XT genome with concatenated scaffolds:
```
/home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat_goodnamez.fasta
```

# XT genome seq
Extract reads from bam file; first sorting by query using samtools
```
samtools sort -n -o aln.qsort.bam aln.bam
```
then extract paired reads using bedtools
```
$ bedtools bamtofastq -i aln.qsort.bam \
                      -fq aln.end1.fq \
                      -fq2 aln.end2.fq
```
as detailed here:
https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html

# Align and genotyping
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts/2020_align_paired_fq_to_ref.sh
```
```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=64:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2018_Austin_XB_genome/Austin_genome/Xbo.v1.fa.gz patht
ofqfilez
# or for XL genome use:
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa.gz pathtofq
filez

module load bwa/0.7.17
module load samtools/1.10


for file in ${2}/*trim_R1.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	#${file::-9}
	echo bwa mem ${1} ${file::-13}trim_R1.fq.gz ${file::-13}trim_R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - 
-o ${file::-9}_sorted.bam
	bwa mem ${1} ${file::-13}trim_R1.fq.gz ${file::-13}trim_R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${
file::-13}_sorted.bam
	samtools index ${file::-13}_sorted.bam
  fi
done

```
# readgroups
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts/2021_picard_add_read_groups.sh
```
```
#!/bin/sh
#SBATCH --job-name=readgroups
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=readgroups.%J.out
#SBATCH --error=readgroups.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2021_picard_add_read_groups.sh /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapt
ed_by_species_across_three_plates/clivii/ 

module load picard/2.23.3

for file in ${1}*_sorted.bam
do
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${file} O=${file}_rg.bam RGID=4 RGLB=$(basename $file) RGPL=I
LLUMINA RGPU=$(basename $file) RGSM=$(basename $file)
done

module load StdEnv/2020 samtools/1.12


for file in ${1}*_sorted.bam_rg.bam
do
    samtools index ${file}
done


```

# haplotype caller
```
/home/ben/projects/rrg-ben/ben/2021_M_f_aurea/ben_scripts/2021_HaplotypeCaller.sh
```
```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --output=HaplotypeCaller.%J.out
#SBATCH --error=HaplotypeCaller.%J.err
#SBATCH --account=def-ben


# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 

# execute like this:
# sbatch 2021_HaplotypeCaller.sh /home/ben/projects/rrg-ben/ben/2021_rhemac_v10/rheMac10.fa.sa path chr

module load nixpkgs/16.09 gatk/4.1.0.0

for file in ${2}*_sorted.bam_rg.bam
do
    gatk --java-options -Xmx8G HaplotypeCaller  -I ${file} -R ${1} -L ${3} -O ${file}_${3}.g.vcf -ERC GVCF
done
```
# CombineGVCF
```
sbatch /home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts/2021_CombineGVCFs.sh /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta ./ Chr7:1-133565930
```
# GenotypeGVCFs
```
sbatch /home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts/2021_GenotypeGVCFs.sh /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta allsites_Chr7:1-133565930.g.vcf.gz
```

# Liberia apomorphies (Liberia_apomorphies.pl)
This script tabulates apomorphies in pairwise comparisons between three taxa:
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
#  Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and searches for apomorphies in each of three taxa

# to execute type Liberia_apomorphies.pl inputfile.tab 12000000000000000300 apomorphies.out 
# where 12000000000000000300 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) individ1, (2) individ2, (3) individ3, or (0) not included

# ./Liberia_apomorphies.pl combined_Chr10.g.vcf.gz_Chr10_GenotypedSNPs.vcf.gz_filtered.vcf.gz.tab 00100000000000000230 apomorphies_00100000000000000230.out

# ./Liberia_apomorphies.pl combined_Chr10.g.vcf.gz_Chr10_GenotypedSNPs.vcf.gz_filtered.vcf.gz.tab 00001000000000000230 apomorphies_00001000000000000230.out

# ./Liberia_apomorphies.pl combined_Chr10.g.vcf.gz_Chr10_GenotypedSNPs.vcf.gz_filtered.vcf.gz.tab 00101000000000000030 apomorphies_00101000000000000030.out


my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile = $ARGV[2];


unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";
print OUTFILE "indiv1\tindiv2\tindiv3\n";

my @sexes = split("",$ARGV[1]);
my @temp;
my $number_of_individuals_included=0;
my $divergence_12_3=0;
my $divergence_13_2=0;
my $divergence_23_1=0;
my @indiv1=();
my @indiv2=();
my @indiv3=();
my $y;
my @unique_indiv1_nucleotides;
my @unique_indiv2_nucleotides;
my @unique_indiv3_nucleotides;
my $counter = 0;


for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if(($sexes[$y] == 1)||($sexes[$y] == 2)||($sexes[$y] == 3)){
		$number_of_individuals_included +=1;
	}	
}	
print "This number should be three: ",$number_of_individuals_included,".\n";

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if($temp[0] ne '#CHROM'){
		if($#temp ne (($#sexes+1)*2)+2){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n";
		}

		# parse the bases in all genotypes in each sex
		@indiv1=();
		@indiv2=();
		@indiv3=();

		$counter=0;
		for ($y = 3 ; $y <= $#temp; $y=$y+2 ) {
			if(($temp[$y] ne ".")&&($temp[$y+1] ne ".")){ # allowing asterisks (*), which are deletions
				if($sexes[$counter] == 1){
						push(@indiv1, $temp[$y]);
						push(@indiv1, $temp[$y+1]);
				}
				elsif($sexes[$counter] == 2){
					push(@indiv2, $temp[$y]);
					push(@indiv2, $temp[$y+1]);
				}
				elsif($sexes[$counter] == 3){
					push(@indiv3, $temp[$y]);
					push(@indiv3, $temp[$y+1]);
				}
			}
			$counter+=1;
		}	
		# OK I should have all the bases loaded for non-missing genotypes for each male and each female
		# and combined as well
		
		# find out what and how many unique nucleotides are in each sex
		@unique_indiv1_nucleotides = uniq @indiv1;
		@unique_indiv2_nucleotides = uniq @indiv2;
		@unique_indiv3_nucleotides = uniq @indiv3;
		if(($#unique_indiv1_nucleotides == 0)&&($#unique_indiv2_nucleotides == 0)&&($#unique_indiv3_nucleotides == 0)){
			# this means that all individuals are homozygous
			if(($unique_indiv1_nucleotides[0] eq $unique_indiv2_nucleotides[0])&&
			($unique_indiv1_nucleotides[0] ne $unique_indiv3_nucleotides[0])){
				$divergence_12_3+=1
			}	
			if(($unique_indiv1_nucleotides[0] eq $unique_indiv3_nucleotides[0])&&
			($unique_indiv1_nucleotides[0] ne $unique_indiv2_nucleotides[0])){
				$divergence_13_2+=1
			}	
			if(($unique_indiv2_nucleotides[0] eq $unique_indiv3_nucleotides[0])&&
			($unique_indiv2_nucleotides[0] ne $unique_indiv1_nucleotides[0])){
				$divergence_23_1+=1
			}	
		} # end of check that there is at least one genotype in each sex
	} # end of check to see if we are at the first line	
	else{ # print the names of the included samples to the outfile
		for ($y = 0 ; $y <= $#sexes ; $y++ ) {
			if(($sexes[$y] == 1)||($sexes[$y] == 2)){
				print OUTFILE $temp[$y+3],"\t";
			}	
			if(($sexes[$y] == 3)){
				print OUTFILE $temp[$y+3],"\n";
			}	
		}	
	}
} # end while
print OUTFILE "divergence_12_3 ",$divergence_12_3,"\n";
print OUTFILE "divergence_13_2 ",$divergence_13_2,"\n";
print OUTFILE "divergence_23_1 ",$divergence_23_1,"\n";
close OUTFILE;
```

# Working directory:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_vcfs_before_filtering/
```
