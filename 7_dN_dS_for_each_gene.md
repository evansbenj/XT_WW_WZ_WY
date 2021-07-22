# First get the coordinates of each exon within each gene

I wrote a perl script that does this; it will output coordinates separately for each transcript, including separate transcript files for each gene if multiple transcripts exist:
```
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in gff file and outputs
# the coordinates of all CDS for each gene

# to execute type ./Get_coordinates_of_CDS_in_each_gene.pl XENTR_10.0_Xenbase_Chr7_lt_30Mb.gff3 exon_coordinates.txt 


my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my @temp;
my @temp1;
my @Gene_name;
my @Gene_ID;
my $Gene_name;
my $Gene_ID;
my $Transcript_ID;
my %gene_hash;
my $gene_counter=0;
my $exon_counter;
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	#print $temp[1],"\n";
	if($temp[2] eq 'gene'){ # we need to parse the gene name and geneID
		# We are assuming that all CDs are preceded by a "gene" annotation
		@temp1=split(';',$temp[8]);
		@Gene_name=split('=',$temp1[1]);
		$Gene_name=$Gene_name[1];
		@Gene_ID=split('ID:',$temp1[2]);
		$Gene_ID=$Gene_ID[1];
		$gene_counter+=1;
		$exon_counter=1;
	}
	elsif($temp[2] eq 'CDS'){ # we need to save the coordinates and phase of this CDS
		# first get the transcript id
		@temp1=split('protein_id=',$temp[8]);
		if(exists $temp1[1]){
			$Transcript_ID=$temp1[1];
		}
		else{
			$Transcript_ID=$Gene_ID;
		}
		$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[0]=$temp[3]; # start coordinate
		$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[1]=$temp[4];	# stop coordinate
		$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[2]=$temp[7];	# phase sort of
								#'0' indicates that the first base of the feature is the first base of a codon, 
								#'1' that the second base is the first base of a codon,
								#'2' that the third base is the first base of a codon,
		$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[4]=$temp[0]; # chromosome						
		if($temp[6] eq '+'){
			$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[3]=1;	# forward orientation

		}
		elsif($temp[6] eq '-'){
			$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[3]=-1;	# reverse orientation
		}
		else{
			print "WTF ",$temp[6],"\n";
		}
		$exon_counter+=1;
	}
} # end while	
close DATAINPUT;	
# OK, now all the CDS are in a hash
my $switch=0;
# print a bed file for each transcript (sometimes there will be multiple transcripts for the same gene)
foreach my $key (sort keys %gene_hash){
	foreach my $transcript_id (sort keys %{$gene_hash{$key}}){
		foreach my $exon (sort {$a <=> $b} keys %{$gene_hash{$key}{$transcript_id}}){
			if($switch == 0){ # open a new file for this transcript
				$switch=1;
				if($gene_hash{$key}{$transcript_id}{$exon}[3] eq 1){ # assume all exons are in this orientation
					unless (open(OUTFILE, ">".$key."_".$transcript_id.".coord"))  {
						print "I can\'t write to $outputfile\n";
						exit;
					}
					print "Creating output file: ".$key."_".$transcript_id.".coord\n";
				}
				else{ # assume all exons are in this orientation
					unless (open(OUTFILE, ">".$key."_".$transcript_id."_rc.coord"))  {
						print "I can\'t write to $outputfile\n";
						exit;
					}
					print "Creating output file: ".$key."_".$transcript_id."_rc.coord\n";
				}	
			}	
			print OUTFILE $gene_hash{$key}{$transcript_id}{$exon}[4]," ";
			print OUTFILE $gene_hash{$key}{$transcript_id}{$exon}[0]," ";
			print OUTFILE $gene_hash{$key}{$transcript_id}{$exon}[1],"\n";
		}
		close OUTFILE;
		$switch=0; # reset for next transcript
	}		
}	

```

# Now subset the exons from the chr7 vcf file, and concatenate them into individual vcfs for each gene

This is done like this:
```
./Run_bcftools_with_lots_of_inputs.pl ../dNdS/gene_beds_chr7_1_30Mb
```

Where `Run_bcftools_with_lots_of_inputs.pl` is:
```
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in coordinate files from a directory
# and feeds them into a bash script that extracts sections using bcftools

# to execute type ./Run_bcftools_with_lots_of_inputs.pl path_to_bed_files


my $inputfile = $ARGV[0];
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}

my @files = glob($inputfile.'/*coord');

foreach ( @files ) {
	print $_,"\n";
	system( "./2021_bcftools_extract_sections_from_vcf.sh $_")
}
```
and `2021_bcftools_extract_sections_from_vcf.sh` is:
```
#!/bin/sh
#SBATCH --job-name=bcftools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:10:00
#SBATCH --mem=2gb
#SBATCH --output=bcftools.%J.out
#SBATCH --error=bcftools.%J.err
#SBATCH --account=def-ben

# execute like this: ./2021_bcftools_extract_sections_from_vcf.sh path_and_filename_of_coordinate_file

module load StdEnv/2020 gcc/9.3.0 bcftools/1.11
bcftools view -R ${1} ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR.vcf.gz -o ${1}.vcf
```

# Generate tab file for each gene

This is done like this:
```
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in coordinate files from a directory
# and feeds them into a bash script that extracts sections using bcftools

# before executing load modules
# module load nixpkgs/16.09 intel/2018.3 vcftools/0.1.16
# to execute type ./Make_lots_of_tab_files.pl path_to_vcf_files


    
my $inputfile = $ARGV[0];
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}

my @files = glob($inputfile.'/*vcf');

foreach ( @files ) {
 #   print $_,"\n";
	system( "vcf-to-tab < $_ > $_.tab")
}
```



These tab file can be converted to W, Z, and Y sequences and piped to paml using a modification of this script:
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );


#  This program reads in a tab delimited genotype file generated
#  from vcftools that has a WW, WZ and WY individual
#  and prints out a fasta file with an interpretation for the 
#  W, Z, and Y chromosome sequence

# Data should be from one chromosome.

# the first individual is XT10_WZ
# the second individual is XT11_WW
# the third individual is XT7_WY

# to execute type Makes_W_Z_Y_from_tab.pl inputfile.tab output_paml_in 


my $inputfile = $ARGV[0];
my $outputfile2 = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}


unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2\n";
	exit;
}
print "Creating output file: $outputfile2\n";

my $W="";
my $Z="";
my $Y="";
my $allele_1;
my $allele_2;
my $allele_3;
my $allele_4;
my $allele_5;
my $allele_6;
my @allele1;
my @allele2;
my @allele3;
my @allele4;
my @allele5;
my @allele6;
my @temp;
my @lengths;
my $max;
my $a;
my $b;
my $start;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		# print $temp[8],"\n";
		# load the alleles
			# load each of the six alleles
			$allele_1 = $temp[3];
			$allele_2 = $temp[4];
			$allele_3 = $temp[5];
			$allele_4 = $temp[6];
			$allele_5 = $temp[7];
			$allele_6 = $temp[8];
			# remove trailing positions if the string is more than one base long
					
			#	$allele_1 = substr($allele_1,0,1);
			#	$allele_2 = substr($allele_2,0,1);
			#	$allele_3 = substr($allele_3,0,1);
			#	$allele_4 = substr($allele_4,0,1);
			#	$allele_5 = substr($allele_5,0,1);
			#	$allele_6 = substr($allele_6,0,1);

			# now see which is the longest
			@lengths=();
			push(@lengths,length($allele_1));
			push(@lengths,length($allele_2));
			push(@lengths,length($allele_3));
			push(@lengths,length($allele_4));
			push(@lengths,length($allele_5));
			push(@lengths,length($allele_6));
			$max = max @lengths;
			# print $max;
			# split each allele into an array
			@allele1=split(//,$allele_1);
			@allele2=split(//,$allele_2);
			@allele3=split(//,$allele_3);
			@allele4=split(//,$allele_4);
			@allele5=split(//,$allele_5);
			@allele6=split(//,$allele_6);	
			if(($#allele1+1) < $max){ # we need padding for this allele; cycle through indexes
									  # to one less than the $max length, which will give $max indexes.
				$start=$#allele1+1;
				for($a=$start; $a<$max; $a++){
					$allele1[$a] = "-";
				}
			}	
			if(($#allele2+1) < $max){ # we need padding for this allele; cycle through indexes
									  # to one less than the $max length, which will give $max indexes.
				$start=$#allele2+1;
				for($a=$start; $a<$max; $a++){
					$allele2[$a] = "-";
				}
			}	
			if(($#allele3+1) < $max){ # we need padding for this allele; cycle through indexes
									  # to one less than the $max length, which will give $max indexes.
				$start=$#allele3+1;
				for($a=$start; $a<$max; $a++){
					$allele3[$a] = "-";
				}
			}	
			if(($#allele4+1) < $max){ # we need padding for this allele; cycle through indexes
									  # to one less than the $max length, which will give $max indexes.
				$start=$#allele4+1;
				for($a=$start; $a<$max; $a++){
					$allele4[$a] = "-";
				}
			}	
			if(($#allele5+1) < $max){ # we need padding for this allele; cycle through indexes
									  # to one less than the $max length, which will give $max indexes.
				$start=$#allele5+1;
				for($a=$start; $a<$max; $a++){
					$allele5[$a] = "-";
				}
			}	
			if(($#allele6+1) < $max){ # we need padding for this allele; cycle through indexes
									  # to one less than the $max length, which will give $max indexes.
				$start=$#allele6+1;
				for($a=$start; $a<$max; $a++){
					$allele6[$a] = "-";
				}
			}	

			# now the lengths should all be the same
			# print $#allele1," ",$#allele2," ",$#allele3," ",$#allele4," ",$#allele5," ",$#allele6,"\n";

			# add the data to the chrs
				for($a=0; $a<=$#allele1; $a++){
					# Reassign each position to the allele
					$allele_1=$allele1[$a];
					$allele_2=$allele2[$a];
					$allele_3=$allele3[$a];
					$allele_4=$allele4[$a];
					$allele_5=$allele5[$a];
					$allele_6=$allele6[$a];
					####################
					# first consider W
					####################
					if($allele_3 eq $allele_4){ # the W genotype is homozygous
						if(($allele_3 eq 'A')||($allele_3 eq 'C')||($allele_3 eq 'G')||($allele_3 eq 'T')){
							$W = $W.$allele_3;
						}
						elsif($allele_3 eq '*'){
							$W = $W.'-';
						}
						else{
							$W = $W.'N';
						}
					}
					elsif($allele_3 ne $allele_4){ # the W genotype is heterozygous
						if((($allele_3 eq 'A')||($allele_3 eq 'C')||($allele_3 eq 'G')||($allele_3 eq 'T'))
							&&
							(($allele_4 eq 'A')||($allele_4 eq 'C')||($allele_4 eq 'G')||($allele_4 eq 'T'))
							){ # both alleles are nucleotides
							# insert an IUPAC symbol
							if((($allele_3 eq 'A')&&($allele_4 eq 'C'))||(($allele_3 eq 'C')&&($allele_4 eq 'A'))){
								$W = $W.'M';
							}	
							elsif((($allele_3 eq 'A')&&($allele_4 eq 'T'))||(($allele_3 eq 'T')&&($allele_4 eq 'A'))){
								$W = $W.'W';
							}	
							elsif((($allele_3 eq 'A')&&($allele_4 eq 'G'))||(($allele_3 eq 'G')&&($allele_4 eq 'A'))){
								$W = $W.'R';
							}	
							elsif((($allele_3 eq 'C')&&($allele_4 eq 'G'))||(($allele_3 eq 'G')&&($allele_4 eq 'C'))){
								$W = $W.'S';
							}	
							elsif((($allele_3 eq 'C')&&($allele_4 eq 'T'))||(($allele_3 eq 'T')&&($allele_4 eq 'C'))){
								$W = $W.'Y';
							}	
							elsif((($allele_3 eq 'G')&&($allele_4 eq 'T'))||(($allele_3 eq 'T')&&($allele_4 eq 'G'))){
								$W = $W.'K';
							}	
						}
						else{ # a het with a dot means unknown, or any other weirdness will be set to be unknown
							$W = $W.'N';
						}

					}	
					####################
					# Now consider the Z
					####################
					if($allele_1 eq $allele_2){ # the Z genotype is homozygous
						# Confirm that the alleles are the same as the W and only use them if they are
						if($allele_1 eq $allele_3){
							if(($allele_1 eq 'A')||($allele_1 eq 'C')||($allele_1 eq 'G')||($allele_1 eq 'T')){
								$Z = $Z.$allele_1;
							}
							elsif($allele_1 eq '*'){
								$Z = $Z.'-';
							}
							else{
								$Z = $Z.'N';
							}
						}
						else{ # we don't know what the Z is because both alleles are different from the W
							$Z = $Z.'N';
						}
					}
					elsif($allele_1 ne $allele_2){ # the Z genotype is heterozygous
						if((($allele_1 eq 'A')||($allele_1 eq 'C')||($allele_1 eq 'G')||($allele_1 eq 'T'))
							&&
							(($allele_2 eq 'A')||($allele_2 eq 'C')||($allele_2 eq 'G')||($allele_2 eq 'T'))
							){ # both Z alleles are nucleotides
							# check if only one is the same as the W
							if((($allele_1 ne $allele_3)||($allele_1 ne $allele_4))&&
								(($allele_2 eq $allele_3)||($allele_2 eq $allele_4))){
								$Z = $Z.$allele_1;
							}
							elsif((($allele_2 ne $allele_3)||($allele_2 ne $allele_4))&&
								(($allele_1 eq $allele_3)||($allele_1 eq $allele_4))){
								$Z = $Z.$allele_2;
							}
							else{ # the WZ genotype is completely different from the WW genotype
								# so we are unable to make a call
								$Z = $Z.'N';
							}
						}
						elsif(($allele_1 eq '*')||($allele_2 eq '*')){ # one allele on the Z is a gap
							if((($allele_1 ne $allele_3)||($allele_1 ne $allele_4))&&
								(($allele_2 eq $allele_3)||($allele_2 eq $allele_4))){
								if($allele_1 eq '*'){
									$Z = $Z.'-'; # the Z has a gap and the W does not
								}
								else{
									$Z = $Z.$allele_1; # the W has a gap and the Z does not
								}	
							}
							elsif((($allele_2 ne $allele_3)||($allele_2 ne $allele_4))&&
								(($allele_1 eq $allele_3)||($allele_1 eq $allele_4))){
								if($allele_2 eq '*'){
									$Z = $Z.'-'; # the Z has a gap and the W does not
								}
								else{
									$Z = $Z.$allele_2; # the W has a gap and the Z does not
								}	
							}
							else{ # the WZ genotype is completely different from the WW genotype
								$Z = $Z.'N';
							}
						}
						else{ # a homozygous dot means unknown, or any other weirdness will be set to be unknown
							$Z = $Z.'N';
						}
					}	

					####################
					# Now consider the Y
					####################
					if($allele_5 eq $allele_6){ # the Y genotype is homozygous
						# Confirm that the alleles are the same as the W and only use them if they are
						if($allele_5 eq $allele_3){
							if(($allele_5 eq 'A')||($allele_5 eq 'C')||($allele_5 eq 'G')||($allele_5 eq 'T')){
								$Y = $Y.$allele_5;
							}
							elsif($allele_5 eq '*'){
								$Y = $Y.'-';
							}
							else{
								$Y = $Y.'N';
							}
						}
						else{ # we don't know what the Y is because both alleles are different from the W
							$Y = $Y.'N';
						}
					}
					elsif($allele_5 ne $allele_6){ # the Y genotype is heterozygous
						if((($allele_5 eq 'A')||($allele_5 eq 'C')||($allele_5 eq 'G')||($allele_5 eq 'T'))
							&&
							(($allele_6 eq 'A')||($allele_6 eq 'C')||($allele_6 eq 'G')||($allele_6 eq 'T'))
							){ # both Y alleles are nucleotides
							# check if only one is the same as the W
							if((($allele_5 ne $allele_3)||($allele_5 ne $allele_4))&&
								(($allele_6 eq $allele_3)||($allele_6 eq $allele_4))){
								$Y = $Y.$allele_5;
							}
							elsif((($allele_6 ne $allele_3)||($allele_6 ne $allele_4))&&
								(($allele_5 eq $allele_3)||($allele_5 eq $allele_4))){
								$Y = $Y.$allele_6;
							}
							else{ # the WY genotype is completely different from the WW genotype
								# so we are unable to make a call
								$Y = $Y.'N';
							}
						}
						elsif(($allele_5 eq '*')||($allele_6 eq '*')){ # one allele on the Y is a gap
							if((($allele_5 ne $allele_3)||($allele_5 ne $allele_4))&&
								(($allele_6 eq $allele_3)||($allele_6 eq $allele_4))){
								if($allele_5 eq '*'){
									$Y = $Y.'-'; # the Y has a gap and the W does not
								}
								else{
									$Y = $Y.$allele_5; # the W has a gap and the Y does not
								}	
							}
							elsif((($allele_6 ne $allele_3)||($allele_6 ne $allele_4))&&
								(($allele_5 eq $allele_3)||($allele_5 eq $allele_4))){
								if($allele_6 eq '*'){
									$Y = $Y.'-'; # the Y has a gap and the W does not
								}
								else{
									$Y = $Y.$allele_6; # the W has a gap and the Y does not
								}	
							}
							else{ # the WY genotype is completely different from the WW genotype
								$Y = $Y.'N';
							}
						}
						else{ # a homozygous dot means unknown, or any other weirdness will be set to be unknown
							$Y = $Y.'N';
						}
					}	
				}	# end of loop across all bases in this position
	} # end if to check for header of input file
} # end while	
close DATAINPUT;				

# print out the lengths of each chr to see if we missed anything
# print "Length of W :",length($W),"\n";
# print "Length of Z :",length($Z),"\n";
# print "Length of Y :",length($Y),"\n";

# reversecomplement if needed
if(substr($inputfile,-17) eq "_rc.coord.vcf.tab"){
	my $W_revcom = reverse $W;
	my $Z_revcom = reverse $Z;
	my $Y_revcom = reverse $Y;
	$W_revcom =~ tr/ACGTacgt/TGCAtgca/;
	$Z_revcom =~ tr/ACGTacgt/TGCAtgca/;
	$Y_revcom =~ tr/ACGTacgt/TGCAtgca/;
	$W = $W_revcom;
	$Z = $Z_revcom;
	$Y = $Y_revcom;
}


# OK print out the fasta file

print OUTFILE2 "3 ",length($W),"\n";
print OUTFILE2 "W_chr     ";
print OUTFILE2 $W,"\n";
print OUTFILE2 "Z_chr     ";
print OUTFILE2 $Z,"\n";
print OUTFILE2 "Y_chr     ";
print OUTFILE2 $Y,"\n";



close OUTFILE2;
my $n = length($W)/3;

if( $n != int(length($W)/3) ){
	print "This infile not multiples of 3: ",$inputfile," ",length($W),"\n";
}


```
I piped all the tab files to this perl script like this:
` ./Make_lots_of_paml_files.pl ../dNdS/gene_beds_chr7_1_30Mb/vcfs_indiv_genes_0_30Mb`
where `Make_lots_of_paml_files.pl ` looks like this:
```
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in coordinate files from a directory
# and feeds them into a bash script that extracts sections using bcftools

# before executing load modules
# module load nixpkgs/16.09 intel/2018.3 vcftools/0.1.16
# to execute type ./Make_lots_of_tab_files.pl path_to_vcf_files


    
my $inputfile = $ARGV[0];
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}

my @files = glob($inputfile.'/*tab');

foreach ( @files ) {
 #   print $_,"\n";
	system( "./Makes_W_Z_Y_from_tab.pl $_ $_\.paml_in")
}	

```

# Analysis with paml

I ran lots of paml runs using the perl script below.
```
#!/usr/bin/env perl
use strict;
use warnings;


# this program will read in lots of file names that are
# paml input files and run them by coping them to a temp
# input file, running the program, and renaming the output file

# before executing load modules
# module load StdEnv/2020  gcc/9.3.0 paml/4.9j
# to execute type ./Execute_lots_of_paml_runs.pl path_to_paml_input_files


    
my $inputfile = $ARGV[0];
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}

my @files = glob($inputfile.'/*paml_in');

foreach ( @files ) {
    system("echo; echo $_; echo");
    system("scp $_ temp.in");
    system("echo Y | codeml ./codeml.ctl");
    system("mv temp.mlc $_.mlc");	
}	
```
Check for empty files (which are due to input files that are not multiples of 3 because of indel weirdness in the vcf and tab files:
```
find . -name '*.mlc' -size 0 | wc -l
```
