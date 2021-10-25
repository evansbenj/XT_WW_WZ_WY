# Kmers with meryl

Working in this directory on graham:
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/ben_scripts
```

I installed meryl from github (https://github.com/marbl/meryl), but first had to load StdEnv/2020

```
module load StdEnv/2020
```

On graham, meryl is here:
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin
```
I made meryl kmer databases like this:

```
#!/bin/sh
#SBATCH --job-name=meryl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl.%J.out
#SBATCH --error=meryl.%J.err
#SBATCH --account=def-ben


# sbatch 2020_meryl_make_kmerdb.sh fastqfile
# ../raw_data/NS.1462.004.IDT_i7_102---IDT_i5_102.XT7_WY_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_102---IDT_i5_102.XT7_WY_R2.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_trim.R1_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_trim.R2_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_WZ_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_WZ_R2.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_trim.R1_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_trim.R2_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_WW_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_WW_R2.fastq.gz


/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl count ${1} threads=4 memory=128 k=29 ou
tput ${1}_meryldb.out
```

## Make a union-sum of the kmer-dbs of the forward and reverrse reads:
```
#SBATCH --job-name=meryl                                                                                                        
#SBATCH --nodes=1                                                                                                               
#SBATCH --ntasks-per-node=1                                                                                                     
#SBATCH --time=48:00:00                                                                                                         
#SBATCH --mem=128gb                                                                                                             
#SBATCH --output=meryl.%J.out                                                                                                   
#SBATCH --error=meryl.%J.err                                                                                                    
#SBATCH --account=def-ben                                                                                                       


# sbatch 2020_meryl_union_kmer_dbs.sh db1 db2 out 

# sbatch 2020_meryl_union_kmer_dbs.sh ../raw_data/XT10_WZ_trim.R1.fq.gz_meryldb.out ../raw_data/XT10_WZ_trim.R2.fq.gz_meryldb.out\
# ../raw_data/XT10_WZ_R1R2_meryldb.out

# sbatch 2020_meryl_union_kmer_dbs.sh ../raw_data/XT11_WW_trim.R1.fq.gz_meryldb.out ../raw_data/XT11_WW_trim.R2.fq.gz_meryldb.out\
# ../raw_data/XT11_WW_R1R2_meryldb.out

# sbatch 2020_meryl_union_kmer_dbs.sh ../raw_data/XT7_WY_trim.R1.fq.gz_meryldb.out ../raw_data/XT7_WY_trim.R2.fq.gz_meryldb.out .\
# ./raw_data/XT7_WY_R1R2_meryldb.out

/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl union-sum ${1} ${2} threads=4 memory=128 k=29 output \
${3}
```
## Subtract kmers to get kmer db with sex chr specific kmers:

```
sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT7_WY_R1R2_meryldb.out ../raw_data/XT11_WW_R1R2_meryldb.out ../raw_data/XT7_WY_minus_XT11_WW_putative_Y_specific.out
sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT10_WZ_R1R2_meryldb.out ../raw_data/XT11_WW_R1R2_meryldb.out ../raw_data/XT10_WZ_minus_XT11_WW_putative_Z_specific.out
```
where `2020_meryl_difference_kmer_dbs.sh` is:
```
#!/bin/sh
#SBATCH --job-name=meryl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=132gb
#SBATCH --output=meryl.%J.out
#SBATCH --error=meryl.%J.err
#SBATCH --account=def-ben


# sbatch 2020_meryl_difference_kmer_dbs.sh db1 db2 out

# sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT7_WY_R1R2_meryldb.out ../raw_data/XT11_WW_R1R2_meryldb.out ../raw_data/XT7_WY_minus_XT11_WW_putative_Y_specific.out
# sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT10_WZ_R1R2_meryldb.out ../raw_data/XT11_WW_R1R2_meryldb.out ../raw_data/XT10_WZ_minus_XT11_WW_putative_Z_specific.out

/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl difference ${1} ${2} threads=4 memory=128 k=29 output ${3}
```
After this try 
```
sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT7_WY_minus_XT11_WW_putative_Y_specific.out ../raw_data/XT10_WZ_R1R2_meryldb.out ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_really_Y_specific.out  

sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT10_WZ_minus_XT11_WW_putative_Z_specific.out ../raw_data/XT7_WY_R1R2_meryldb.out ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out

sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT11_WW_R1R2_meryldb.out ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out ../raw_data/XT11_WW_R1R2_minus_Z_specific_putative_W_specific_meryldb.out

sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT11_WW_R1R2_minus_Z_specific_putative_W_specific_meryldb.out ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_really_Y_specific.out ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out

```

## getting fastq from bam for v10 ZW

```
module load nixpkgs/16.09  intel/2018.3 bamutil/1.0.14
```
```
bam bam2FastQ --in JBL052.bam --outBase XXX
```
## make kmer db from forward and rev reads:
```
sbatch 2020_meryl_make_kmerdb.sh ../../2020_XT_v10_raw_data/XT-v10_rawdata/JBL052__1.fastq
```
```
sbatch 2020_meryl_make_kmerdb.sh ../../2020_XT_v10_raw_data/XT-v10_rawdata/JBL052__2.fastq
```

## combine kmer dbs:
```
sbatch 2020_meryl_union_kmer_dbs.sh ../../2020_XT_v10_raw_data/XT-v10_rawdata/JBL052__1.fastq_meryldb.out ../../2020_XT_v10_raw_data/XT-v10_rawdata/JBL052__2.fastq_meryldb.out ../../2020_XT_v10_raw_data/XT-v10_rawdata/JBL052_R1R2.fastq_meryldb.out
```

## now substract v10_ZW
```
sbatch 2020_meryl_difference_kmer_dbs.sh ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_really_Y_specific.out ../raw_data/XT7_WY_R1R2_meryldb.out ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out
```

## Putative kmer specific dbs
```
../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out
../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out
../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out
```

## print kmer-specific dbs for Z and Y
```
sbatch 2021_meryl_print_kmer_dbs.sh ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out
sbatch 2021_meryl_print_kmer_dbs.sh ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out
sbatch 2021_meryl_print_kmer_dbs.sh ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out
sbatch 2021_meryl_print_kmer_dbs.sh ../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out
```

## filter to include only kmers with more than 2 observations and less than 100 observations (which is way more than the coverage):
```
awk '{ if (($2 > 2)&&($2 < 100)) { print } }' ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed.out > ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out
```
```
awk '{ if (($2 > 2)&&($2 < 100)) { print } }' ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out > ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out
```
```
awk '{ if (($2 > 2)&&($2 < 100)) { print } }' ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out_printed.out > ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out_printed.out_printed_filtered_gt_2_lt_100.out
```
```
awk '{ if (($2 > 2)&&($2 < 100)) { print } }' ../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out_printed.out > ../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out_printed.out_printed_filtered_gt_2_lt_100.out
```


## Make a multifasta file out of the kmer list for use with cookiecutter
```
awk '{print ">\n"$1}' ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out > ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa
```
```
awk '{print ">\n"$1}' ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out > ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out_seqs.fa
```
```
awk '{print ">\n"$1}' ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out_printed.out_printed_filtered_gt_2_lt_100.out > ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa
```
```
awk '{print ">\n"$1}' ../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out_printed.out_printed_filtered_gt_2_lt_100.out > ../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa
```

## extract reads with sex-chr specific kmers
```
sbatch 2021_cookiecutter_extract.sh ../raw_data/XT7_WY_trim_noadapters/XT7_WY_trim_no_adapt.R1.fq ../raw_data/XT7_WY_trim_noadapters/XT7_WY_trim_no_adapt.R2.fq ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_fq_filez
```
```
sbatch 2021_cookiecutter_extract.sh ../raw_data/XT10_WZ_trim_noadapters/XT10_WZ_trim_no_adapt.R1.fq ../raw_data/XT10_WZ_trim_noadapters/XT10_WZ_trim_no_adapt.R2.fq ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out_seqs.fa ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out_seqs.fa_fq_filez
```
```
sbatch 2021_cookiecutter_extract.sh ../raw_data/XT11_WW_trim_noadapters/XT11_WW_trim_no_adapt.R1.fq ../raw_data/XT11_WW_trim_noadapters/XT11_WW_trim_no_adapt.R2.fq ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa ../raw_data/XT11_WW_R1R2_minus_Z_specific_minus_Y_specific_putative_really_W_specific_meryldb.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa_fq_filez
```
```
sbatch 2021_cookiecutter_extract.sh ../raw_data/XT11_WW_trim_noadapters/XT11_WW_trim_no_adapt.R1.fq ../raw_data/XT11_WW_trim_noadapters/XT11_WW_trim_no_adapt.R2.fq ../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa ../raw_data/XT11_WW_minus_XT1_ZY_putative_W_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_seqs.fa_fq_filez
```

# convert fq files to trinity format (with /1 and /2 after the for and rev reads)
```
awk '{ if (NR%4==1) { print $1"_"$2"/1" } else { print } }' Read1.fastq > rename_Read1.fastq
awk '{ if (NR%4==1) { print $1"_"$2"/2" } else { print } }' Read2.fastq > rename_Read2.fastq
```

# Assemble
```
sbatch 2021_trinity.sh ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_fq_filez/PE_SE_combined_for_trinity.fastq ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_fq_filez/XT7_WY_trim_no_adapt.R2.filtered_trinity.fastq ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_fq_filez
```

# Align assembly to ref
```
sbatch 2021_align_multifasta_to_ref.sh /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_fq_filez/trinity.Trinity.fasta ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_fq_filez/XT7_WY_Y_specific_XTv10
```
```
sbatch 2021_align_multifasta_to_ref.sh /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out_seqs.fa_fq_filez/trinity_XT10_WZ.Trinity.fasta ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out_seqs.fa_fq_filez/XT10_WZ_Z_specific_XTv10
```

# Here are the assemblies based on kmer reads
ChrY
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2_lt_100.out_fq_filez/trinity.Trinity.fasta
```
ChrZ
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2_lt_100.out_seqs.fa_fq_filez/trinity_XT10_WZ.Trinity.fasta
```

# Make blast dbs and blast some candidates
```
module load nixpkgs/16.09 gcc/7.3.0 blast+/2.9.0
makeblastdb -in trinity_XT10_WZ.Trinity.fasta -dbtype nucl -out trinity_XT10_WZ.Trinity.fasta_blastable
blastn -query PZP.fasta -db trinity_XT10_WZ.Trinity.fasta_blastable -out PZP.fasta_to_Z_kmer_assembly.out
```

# BELOW NOT USED



# get ave, max, min
```
awk '{print $2}' ../raw_data/XT10_WZ_minus_XT11_WW_minus_XT7_WY_putative_really_Z_specific.out_printed_filtered_gt_2.out | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print "avg " total/count," | max "max," | min " min}'
```
```
awk '{print $2}' ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2.out | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print "avg " total/count," | max "max," | min " min}'
```

# make this output into a multifasta file for cookie cutter
```
awk '{print ">\n"$1}' ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2.out > ../raw_data/XT7_WY_minus_XT11_WW_minus_XT10_WZ_putative_reallyreally_Y_specific.out_printed.out_printed_filtered_gt_2.out_seqs.fa
```



# Extract chr7 from v10 ref
```
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(Chr7)}print if $c' XENTR_10.0_genome.fasta > XT_Chr7.fasta
```

# Focus on first 20million bp

```
head -n 333333 XT_Chr7.fasta > XT_Chr7_first20mil.fasta
```

# Split into 200,000 bp bits
```
./Split_fasta.pl /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XT_Chr7_first20mil.fasta 3333 /home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/Chr7_bits
```
Because there are 60 bp/line, the 3333 means that each bit will have 3333*60=~200,000 bp.

Where `Split_fasta.pl` is this:
```
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a fasta file and splits it into several
# smaller files of length equal to $ARGV[1].  
# This is useful for doing kmer quantification of different parts of a chromosome
# Output files will be printed to $ARGV[2], which is the "output_directory"
# First calculate how many lines the input fasta file is. The divide this number
# by the number of bits you want to divide it into to get the $ARGV[1] length

# Run this script like this
# Split_fasta.pl input.fasta length output_directory 


my $inputfile = $ARGV[0]; # This is the input file
my $length = $ARGV[1]; # this how many lines each subset file will have
my $outputdirectory = $ARGV[2]; # This is where the bits will be printed

unless (open DATAINPUT, $inputfile) {
	print 'Can not find the input file.\n';
	exit;
}

my $filecounter=0;
my $linecounter=0;
my $header;
my @temp;
my $outputfile;
my $line;

# Read in datainput file

while ( my $line = <DATAINPUT>) {
	if($line =~ '>'){
		# save the fasta HEADER for later to print in each output file
		$filecounter+=1;
		# parse header
		chomp($line);
		@temp=split('>',$line);
		$header=$temp[1];
		unless (open(OUTFILE, ">".$outputdirectory."\/".$header."_".$filecounter))  {
			print "I can\'t write to $outputfile\n";
			exit;
		}
		print OUTFILE ">".$header."_".$filecounter."\n";
		$linecounter=0;
	}
	elsif($linecounter < $length){
		print OUTFILE $line;
		$linecounter+=1;
	}	
	else{	
		$filecounter+=1;
		unless (open(OUTFILE, ">".$outputdirectory."\/".$header."_".$filecounter))  {
			print "I can\'t write to $outputfile\n";
			exit;
		}
		print OUTFILE ">".$header."_".$filecounter."\n";
		print OUTFILE $line;
		$linecounter=0;
	}
}		
```

# Make little kmer db out of each bit
```
for file in ../Chr7_bits/*
do
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl count "$file" threads=4 memory=128 k=29 output "$file"_meryldb.out
done
```
