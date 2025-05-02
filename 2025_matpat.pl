# Making a matpat output from angsd

Two scripts. The first one uses sex-specific homozygosity and the second one uses genotypes from parents
```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a vcf file with genotypic information from
# a family and identifies positions that
# are homozygous in all sons with genotypes and het in at least five daughters (mat sites because daughters 
# are heteroz) and positions that are
# are homozygous in all daughters with genotypes and het in at least five sons (pat sites because sons 
# are heteroz) 
# This is for dataset 
# where there is no information from one or both parent (e.g. the dad for calcaratus)
# module load StdEnv/2023 perl/5.36.1
# execute like this:
# ./Gets_matonly_positions_from_vcf_file_homoz_daughters_and_sons.pl vcf 21111111111111111111000000000000000 matout patout

# epitrop RADseq: 111111111120000000002 (10 daughters; 9 sons; exclude parents)
# ./Gets_matpat_positions_from_vcf_file_homoz_daughters_and_sons.pl all_combined.g.vcf.gz_genotyped.vcf.gz 111111111120000000002 epi_mat_allchrs_allsonshomoz_atleastfivedaughterssheteroz.out epi_pat_allchrs_alldaughtersshomoz_atleastfivesonssheteroz.out

# calcaratus RADseq: 21111111111111111111000000000000000 (19 daughters, 15 sons, exclude mom)
# ./Gets_matpat_positions_from_vcf_file_homoz_daughters_and_sons.pl all_combined.g.vcf.gz_genotyped.vcf.gz 21111111111111111111000000000000000  cal_mat_allchrs_allsonshomoz_atleastfivedaughterssheteroz.out cal_pat_allchrs_alldaughtersshomoz_atleastfivesonssheteroz.out

# mel RADseq 21111111112000000000
# ./Gets_matpat_positions_from_vcf_file_homoz_daughters_and_sons.pl combined_Chr7:1-133565930.g.vcf.gz_Chr7:1-133565930_SNPs.vcf.gz 21111111112000000000  mel_mat_allchrs_allsonshomoz_atleastfivedaughterssheteroz_Chr7.out mel_pat_allchrs_alldaughtersshomoz_atleastfivesonssheteroz_Chr7.out

# trop C659 RADseq 110000000022110000000011111111001111111
# ./Gets_matpat_positions_from_vcf_file_homoz_daughters_and_sons.pl combined_Chr7.g.vcf.gz_Chr7_GenotypedSNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz 110000000022110000000011111111001111111  C659_mat_allchrs_allsonshomoz_atleastfivedaughterssheteroz_Chr7.out C659_pat_allchrs_alldaughtersshomoz_atleastfivesonssheteroz_Chr7.out

# where 21111111111111111111000000000000000 is whether the individuals are the excluded(2), a daughter(1), or a son(0)
# and matout and patout are the output files (with chromsoome name in them) 

my $vcf = $ARGV[0];
my $sexes = $ARGV[1];
my $outputfile = $ARGV[2];
my $outputfile2 = $ARGV[3];
my @columns;
my @pat;
my @pat1;
my $x;

my @sexes = split("",$sexes);

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2  $!\n\n";
	exit;
}
print "Creating output file: $outputfile2\n";



if ($vcf =~ /.gz$/) {
	#open DATAINPUT, '<:gzip', $vcf or die "Could not read from $vcf: $!";
	open(DATAINPUT, "gunzip -c $vcf |") || die "can’t open pipe to $vcf";
}
else {
	open DATAINPUT, $vcf or die "Could not read from $vcf: $!";
}

my $number_of_samples=0;
my $switch1=0;
my $switch2=0;
my $num_daughters=0;
my $num_het_daughters=0;
my $num_sons=0;
my $num_het_sons=0;


while ( my $line = <DATAINPUT>) {
	@columns=split("	",$line);
		#print $line,"\n";
		if($columns[0] =~ m/^#/){ # this is a commented line
			if($columns[0] eq '#CHROM'){ # this is the first line
				$number_of_samples = scalar(@columns)-9;
				print "Number of samples ",$number_of_samples,"\n";
			}
		}
		else{ # this is the genotype data
			$switch1=0;
			$switch2=0;
			$num_daughters=0;
			$num_het_daughters=0;
			$num_sons=0;
			$num_het_sons=0;
			for ($x = 1 ; $x <= $number_of_samples; $x++ ){ # cycle through each sample
				if($sexes[$x-1] != 2){ # only consider samples that are included
					@pat = split(":",$columns[$x+8]); # this is the whole genotype and other info from the vcf file
					@pat1 = split(/[\|\/]/,$pat[0]); # this is only the genotype
					# check if the daughters are homozygous or missing
					if(($pat[0] ne './.')&&($pat[0] ne '.|.')&&
					($pat[0] ne '.')&&($sexes[$x-1] == 1)){ # this is a daughter with a genotype
					# this is a daughter
					$num_daughters+=1;
						if(($pat1[0] eq $pat1[1])&&($switch1 != 9)){
							# this daughter is homoz 
							$switch1 = 1;
							#print "goodpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
						}
						elsif(($pat1[0] ne $pat1[1])&&($switch1 != 9)){
							# this daughter is the first heteroz
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$switch1 = 9;
							$num_het_daughters+=1;
						}	
						elsif(($pat1[0] ne $pat1[1])&&($switch1 == 9)){
							# this daughter is heteroz, but not the first one
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$num_het_daughters+=1;
						}	
					}
					# check if the sons are heterozygous or missing
					elsif(($pat[0] ne './.')&&($pat[0] ne '.|.')&&
					($pat[0] ne '.')&&($sexes[$x-1] == 0)){
					# this is a son
					$num_sons+=1;
						if(($pat1[0] eq $pat1[1])&&($switch2 != 9)){
							# this son is homoz 
							$switch2 = 1;
							#print "goodpats ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
						}
						elsif(($pat1[0] ne $pat1[1])&&($switch2 != 9)){
							# this son is the first heteroz
							#print "badpats ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$switch2 = 9;
							$num_het_sons+=1;
						}
						elsif(($pat1[0] ne $pat1[1])&&($switch2 == 9)){
							# this son is heteroz, but not the first one
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$num_het_sons+=1;
						}	
					}	
				}
				#print "\n";
			}
			#if(($switch1 == 1)&&($switch2 == 9)&&($num_sons == $num_het_sons)){
				# all the daughters with genotypes are homoz
				# all the sons with genotypes are heterozygous
			#if(($switch1 == 1)&&($switch2 == 9)){	
				# all the daughters with genotypes are homoz
				# at least one son with a genotypes is heterozygous
			if(($switch1 == 1)&&($switch2 == 9)&&($num_het_sons >=5)){		
				# all the daughters with genotypes are homoz
				# at least five sons with a genotypes are heterozygous
				print OUTFILE2 $columns[0],"\t",$columns[1],"\n";
			}
			#elsif(($switch2 == 1)&&($switch1 == 9)&&($num_daughters == $num_het_daughters)){
				# all the sons with genotypes are homoz
				# all the daughters with genotypes are heterozygous
			#elsif(($switch2 == 1)&&($switch1 == 9)){
				# all the sonz with genotypes are homoz
				# at least one daughter with a genotype is heterozygous
			elsif(($switch2 == 1)&&($switch1 == 9)&&($num_het_daughters >=5)){				
				# all the sonz with genotypes are homoz
				# at least five daughters with a genotype are heterozygous
				print OUTFILE $columns[0],"\t",$columns[1],"\n";
			}
			
		} # end else
} # end while
close DATAINPUT;
close OUTFILE;
close OUTFILE2;

```


```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a angsd output file and two other files
# that have lists of mat and pat sites

# execute like this:
# ./Makes_a_mat_and_pat_from_angsdoutput.pl allo_family1_out_additive_F1.lrt0.gz allo_family1_matout.txt allo_family1_patout.txt allo_family1_mat_angsdout allo_family1_pat_angsdout

# where matsites are the sample positions of mother-only het sites 
# and patsites are the sample positions of father-only het sites
# and matangsdout and patangsdout are the filtered mat and pat angsd output files 

my $angsd = $ARGV[0];
my $matsites = $ARGV[1];
my $patsites = $ARGV[2];
my $outputfile = $ARGV[3];
my $outputfile2 = $ARGV[4];
my @columns;
my @mat;
my @pat;
my @mat1;
my @pat1;


unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2  $!\n\n";
	exit;
}
print "Creating output file: $outputfile2\n";



#### first open the matsites and patsites and load this into a hash
my %matsiteshash;
my %patsiteshash;

unless (open DATAINPUT, $matsites) {
	print "Can not find the matsites file, jackass.\n";
	exit;
}
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@columns=split(/\s+/,$line);
	$matsiteshash{$columns[0]."_".$columns[1]} = 1;
	#print $columns[0]."_".$columns[1],"\n";
}	
close DATAINPUT;

unless (open DATAINPUT2, $patsites) {
	print "Can not find the patsites file, jackass.\n";
	exit;
}
while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@columns=split(/\s+/,$line);
	$patsiteshash{$columns[0]."_".$columns[1]} = 1;
}	
close DATAINPUT2;



# now load the angsd output 
if ($angsd =~ /.gz$/) {
	#open DATAINPUT, '<:gzip', $angsd or die "Could not read from $angsd: $!";
	open(DATAINPUT3, "gunzip -c $angsd |") || die "can’t open pipe to $angsd";
}
else {
	open DATAINPUT3, $angsd or die "Could not read from $angsd: $!";
}


while ( my $line = <DATAINPUT3>) {
	@columns=split(/\s+/,$line);
		#print $columns[0]."_".$columns[1],"\n";
		if(defined $matsiteshash{$columns[0]."_".$columns[1]}){
			#print "hello\n";
			$matsiteshash{$columns[0]."_".$columns[1]} = $line;
		}
		elsif(defined $patsiteshash{$columns[0]."_".$columns[1]}){
			$patsiteshash{$columns[0]."_".$columns[1]} = $line;
		} # end else
} # end while
close DATAINPUT3;


# now print the filtered mat outputs
my @keys = sort { $a cmp $b } keys %matsiteshash;

foreach my $key ( @keys ){
		if($matsiteshash{$key} ne "1"){
    		print OUTFILE $matsiteshash{$key};
		}
    }
close OUTFILE;

# now print the filtered pat outputs
@keys = sort { $a cmp $b } keys %patsiteshash;

foreach my $key ( @keys ){
		if($patsiteshash{$key} ne "1"){
    		print OUTFILE2 $patsiteshash{$key};
		}	
    }
close OUTFILE2;
```
