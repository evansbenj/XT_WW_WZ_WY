# Checking for Mendelian inheritance in RADseq data

Some matpat analyses raised the possibility that the epitropicalis father was mixed up or mislabeled. To test this I wrote a script that quanitifies the number of positions that are heterozygous in each parent that are also hets in at least one offspring:
```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a vcf file with genotypic information from
# a family and counts positions that are heterozygous in each parent
# and also heterozygous in at least one offspring

# module load StdEnv/2023 perl/5.36.1
# execute like this:
# ./Quantifies_Mendelian_inheritance_in_both_parents vcf 111111111130000000004 

# epitrop RADseq: 111111111120000000002 (10 daughters; 9 sons; exclude parents)
# ./Quantifies_Mendelian_inheritance_in_both_parents.pl all_combined.g.vcf.gz_genotyped.vcf.gz 111111111130000000004 


# where 111111111120000000003 is whether the individuals are the excluded(2), mom(3), dad(4), a daughter(1), or a son(0)


my $vcf = $ARGV[0];
my $sexes = $ARGV[1];
my @columns;
my @pat;
my @pat1;
my $x;
my $mom_sample;
my $dad_sample;
my $mat_het=0;
my $pat_het=0;
my $number_of_samples=0;

my @sexes = split("",$sexes);



if ($vcf =~ /.gz$/) {
	#open DATAINPUT, '<:gzip', $vcf or die "Could not read from $vcf: $!";
	open(DATAINPUT, "gunzip -c $vcf |") || die "can’t open pipe to $vcf";
}
else {
	open DATAINPUT, $vcf or die "Could not read from $vcf: $!";
}


my $num_mat_hets=0;
my $num_pat_hets=0;
my $num_mat_offspring_hets=0;
my $num_pat_offspring_hets=0;
my $mat_switch=0;
my $pat_switch=0;


while ( my $line = <DATAINPUT>) {
	@columns=split("	",$line);
		#print $line,"\n";
		if($columns[0] =~ m/^#/){ # this is a commented line
			if($columns[0] eq '#CHROM'){ # this is the first line
				$number_of_samples = scalar(@columns)-9;
				print "Number of samples ",$number_of_samples,"\n";
				# find out which sample is the mom and dad
				for ($x = 0 ; $x < $number_of_samples; $x++ ){
					print $sexes[$x]," ";
					if($sexes[$x] == 3){
						print "hello1 $x\n";
						$mom_sample=$x;
					}
					elsif($sexes[$x] == 4){
						print "hello2 $x\n";
						$dad_sample=$x;
					}
				}	
				print $mom_sample,"\t",$dad_sample,"\n";
				print "mat_hets\tpat_hets\n";
			}
		}
		else{ # this is the genotype data
			# determine whether mom and dad are hets
			$mat_het=0;
			$pat_het=0;
			$mat_switch=0;
			$pat_switch=0;
			@pat = split(":",$columns[$mom_sample+9]); # this is the whole genotype and other info from the vcf file
			@pat1 = split(/[\|\/]/,$pat[0]); # this is only the genotype
				if($pat1[0] ne $pat1[1]){
					$mat_het=1;
					$num_mat_hets+=1;
				}	
				@pat = split(":",$columns[$dad_sample+9]); # this is the whole genotype and other info from the vcf file
				@pat1 = split(/[\|\/]/,$pat[0]); # this is only the genotype
				if($pat1[0] ne $pat1[1]){
					$pat_het=1;
					$num_pat_hets+=1;
				}	
				if($mat_het != $pat_het){ # mom or dad is het but not both
					for ($x = 0 ; $x < $number_of_samples; $x++ ){ # cycle through each sample
						if(($sexes[$x] != 2)&&($sexes[$x] != 3)&&($sexes[$x] != 4)){ # only consider offspring
						@pat = split(":",$columns[$x+9]); # this is the whole genotype and other info from the vcf file
						@pat1 = split(/[\|\/]/,$pat[0]); # this is only the genotype
						# check if the daughters are homozygous or missing
							if($pat1[0] ne $pat1[1]){
								# this offspring is heteroz
								if($mat_het==1){
									# flip the mat switch because at least one offspring is het
									$mat_switch=9;									
								}
								elsif($pat_het==1){
									# flip the pat switch because at least one offspring is het
									$pat_switch=9;
								}
							}	
						}
					}
					# now add to hets only if the switch is flipped
					if($mat_switch == 9){
						$num_mat_offspring_hets+=1;
					}
					if($pat_switch == 9){
						$num_pat_offspring_hets+=1;
					}
					
				}	
			} # end else
} # end while
close DATAINPUT;

print $num_mat_offspring_hets," offspring hets out of ",$num_mat_hets," mat hets ",$num_mat_offspring_hets/$num_mat_hets,"\n";
print $num_pat_offspring_hets," offspring hets out of ",$num_pat_hets," pat hets ",$num_pat_offspring_hets/$num_pat_hets,"\n";

```
