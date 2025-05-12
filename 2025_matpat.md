# Making a matpat output from angsd

MatPat using genotypes from parents. The first script is when we have data for both parents; the second one is for Xcal where we only have data for the mom.


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


```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a vcf file with genotypic information from
# a family and identifies positions that
# are heterozygous in only the mother (mat sites).
# We are not identifying pat positions but will print out a blank pat file
# This is for dataset 
# where there is no information from one parent (the dad for calcaratus)
# module load StdEnv/2023 perl/5.36.1
# execute like this:
# ./Gets_matonly_positions_from_vcf_file_nopat.pl vcf mat 20000000000000000000111111111111111 matout patout

# where mat is the sample position of the mother and
# 20000000000000000000111111111111111 is whether the individuals are the mother(2), a daughter(0), or a son(1)
# and matout and patout are the output files (with chromsoome name in them) 

my $vcf = $ARGV[0];
my $mat = $ARGV[1];
my $sexes = $ARGV[2];
my $outputfile = $ARGV[3];
my $outputfile2 = $ARGV[4];
my @columns;
my @mat;
my @pat;
my @mat1;
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



#### Prepare the input file  with values for svl

#unless (open DATAINPUT, $vcf) {
#	print "Can not find the vcf file, jackass.\n";
#	exit;
#}


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


while ( my $line = <DATAINPUT>) {
	@columns=split("	",$line);
		#print $line,"\n";
		if($columns[0] =~ m/^#/){
			if($columns[0] eq '#CHROM'){
				print $columns[$mat+8],"\n";
				$number_of_samples = scalar(@columns)-9;
				print "Number of samples ",$number_of_samples,"\n";
			}
		}
		else{
			@mat = split(":",$columns[$mat+8]);
			#print $mat[0],"\n";
			$switch1=0;
			$switch2=0;
			# select positions that are not missing in the mother
			if(($mat[0] ne './.')&&($mat[0] ne '.|.')&&($mat[0] ne '.')){
				@mat1 = split(/[\|\/]/,$mat[0]);
				if($mat1[0] ne $mat1[1]){ # this is a mat site
					print OUTFILE $columns[0],"\t",$columns[1],"\n";
					print "mat\t",$columns[0],"\t",$columns[1],"\t",$mat1[0],"\t",$mat1[1],"\n";
				}
				elsif($mat1[0] eq $mat1[1]){ # this could be a pat site if there are heterozygous offspring
					for ($x = 1 ; $x <= $number_of_samples; $x++ ){ # cycle through each sample
						if($x ne $mat){ # only consider offspring, not the mother
							@pat = split(":",$columns[$x+8]);
							@pat1 = split(/[\|\/]/,$pat[0]);
							# check if the daughters are homozygous or missing
							if(($pat[0] ne './.')&&($pat[0] ne '.|.')&&
							($pat[0] ne '.')&&($sexes[$x-1] == 0)){ # this is a daughter with a genotype
								if(($mat1[0] eq $pat1[0])&&($pat1[0] eq $pat1[1])&&($switch1 != 9)){
									# the daughter is homoz for the maternal variant and other 
									# daughters inspected so far are too
									$switch1 = 1;
									print "goodpatd ",$columns[0],"\t",$columns[1]," ",$mat[0]," ",$pat[0]," ";
								}
								else{
									# the daughter is homoz for a non-maternal variant or heteroz
									print "badpatd ",$columns[0],"\t",$columns[1]," ",$mat[0]," ",$pat[0]," ";
									$switch1 = 9;
								}	
							}
							# check if the sons are heterozygous or missing
							elsif(($pat[0] ne './.')&&($pat[0] ne '.|.')&&
							($pat[0] ne '.')&&($sexes[$x-1] == 1)){
								# this is a son
								if(($pat1[0] ne $pat1[1])&&($switch2 != 9)){
									# this son is heterozygous and other sons inspected so far are too.
									$switch2 = 1;
								}
								else{
									print "badpats ",$columns[0],"\t",$columns[1]," ",$pat[0]," ";
									# $switch2 = 9;
								}
							}	
						}
						print "\n";
					}
					if(($switch1 == 1)&&($switch2 == 1)){
						# all the daughters are homoz for the maternal variant
						# and all the sons are heterozygous
						#print OUTFILE2 $columns[0],"\t",$columns[1],"\n";
						#print "pat\t",$columns[0],"\t",$columns[1],"\t",$mat1[0],"\t",$mat1[1],"\t",$pat1[0],"\t",$pat1[1],"\n";
					}
				}
			}
		} # end else
} # end while
close DATAINPUT;
close OUTFILE;
close OUTFILE2;
```
