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

# Make a union-sum of the kmer-dbs of the forward and reverrse reads:
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
