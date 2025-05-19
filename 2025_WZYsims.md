# Simulating WZY 

```perl
#!/usr/bin/env perl
use strict;
use warnings;
# use lib qw(~/perl_modules);
# use List::MoreUtils qw/ uniq /;

# This script will run a simulation that explores how sex chromosome
# genotypes change over generations of random mating

# It will start with 50:50 ratio of WW and WZ females and
# a 33% abundance of males with ZZ, ZY, and WY sex chromosomes

# it will print out the following in a file:
# the frequecies of each sex chromosome and
# the frequencies of each sex chromosome genotype

# to execute:
# ./WZY_sims.pl number_of_generations outputfilename

my $number_of_generations = $ARGV[0];
my $outputfile1 = $ARGV[1];

# unless (open DATAINPUT, $inputfile) {
#	print "Can not find the input file.\n";
#	exit;
# }

unless (open(OUTFILE1, ">$outputfile1"))  {
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n";

print OUTFILE1 "Generation WW WZ ZZ ZY WY W Z Y SexRatio\n";


my $population_size=6000; # use a number that is divisible by 2 and also by 3
my $a;
my $x;
my $y;
my @mothers = ();
my @fathers = ();
my $mother;
my $father;
my @daughters = ();
my @sons = ();
my @genotype_countz = ();
my @sex_chr_countz = ();
my $num;


# female genotypes: 0 = WW, 1 = WZ
# male genotypes: 2 = ZZ, 3 = ZY, 4 = YY

# first populate the parent array with an equal frequency of each genotype
# mothers
for ($y = 0 ; $y < $population_size/2 ; $y++ ) {
	push(@mothers, 0); # 0 = WW
	push(@mothers, 1); # 1 = WZ
}	
# fathers
for ($y = 0 ; $y < $population_size/3 ; $y++ ) {
	push(@fathers, 2); # 2 = ZZ
	push(@fathers, 3); # 3 = ZY
	push(@fathers, 4); # 4 = WY
}




# now run the simulation
for ($a = 0 ; $a <= $number_of_generations ; $a++ ) {
	# first generate offspring genotypes
	# allow a sex ratio skew
	@daughters = ();
	@sons = ();
	for ($y = 0 ; $y < $population_size ; $y++ ) {
		# choose parents randomly with replacement - approximates a Poisson distribution of reproductive variation
		$mother = $mothers[rand @mothers];
		#print "hellow ",$mother,"\n";
		$father = $fathers[rand @fathers];
		# make offspring
		if(($mother == 0)&&($father == 2)){ # WW x ZZ the child is a WZ daughter
			push(@daughters,1); # WZ daughter
		}	
		elsif(($mother == 0)&&($father == 3)){ # WW x ZY the child is either a WZ daughter or a WY son
			# flip a coin
			$x = int(rand(2)); # this should generate a random 0 or 1
			if($x == 0){
				push(@daughters,1); # WZ daughter
			}
			elsif($x == 1){
				push(@sons,4); # WY son
			}
		}
		elsif(($mother == 0)&&($father == 4)){ # WW x WY the child is either a WW daughter or a WY son
			# flip a coin
			$x = int(rand(2)); # this should generate a random 0 or 1
			if($x == 0){ 
				push(@daughters,0); # WW daughter
			}
			elsif($x == 1){
				push(@sons,4); # WY son
			}
		}
		elsif(($mother == 1)&&($father == 2)){ # WZ x ZZ the child is a WZ daughter or a ZZ son
			$x = int(rand(2)); # this should generate a random 0 or 1
			if($x == 0){ 
				push(@daughters,1); # WZ daughter
			}
			elsif($x == 1){
				push(@sons,2); # WY son
			}
		}
		elsif(($mother == 1)&&($father == 3)){ # WZ x ZY the child is either a WZ daughter, or a ZZ, WY, or ZY son
			# flip a coin
			$x = int(rand(4)); # this should generate a random 0 or 1 or 2 or 3
			if($x == 0){ 
				push(@daughters,1); # WZ daughter
			}
			elsif($x == 1){
				push(@sons,2); # ZZ son
			}	
			elsif($x == 2){
				push(@sons,4); # WY son
			}
			elsif($x == 3){
				push(@sons,3); # ZY son
			}
		}	
		elsif(($mother == 1)&&($father == 4)){ # WZ x WY the child is a WW or WZ daughter or a WY or a ZY son
			$x = int(rand(4)); # this should generate a random 0 or 1 or 2 or 3
			if($x == 0){ 
				push(@daughters,0); # WW daughter
			}
			elsif($x == 1){
				push(@daughters,1); # WZ daughter
			}	
			elsif($x == 2){
				push(@sons,4); # WY son
			}
			elsif($x == 3){
				push(@sons,3); # ZY son
			}
		}	
	}
	# offspring are generated, now assign them to be parents of the next generation
	#print "hello @daughters goodby @sons";
	
	@mothers = @daughters;
	@fathers = @sons;
	# print "Wow @daughters\n";
	# print information to outfile
	@genotype_countz = (); # 0 = WW, 1 = WZ, 2 = ZZ, 3 = ZY, 4 = WY
	@sex_chr_countz = (); # O = W, 1 = Z, 2 = Y
	# initialize the arrays
	$genotype_countz[0]=0;
	$genotype_countz[1]=0;
	$genotype_countz[2]=0;
	$genotype_countz[3]=0;
	$genotype_countz[4]=0;
	$sex_chr_countz[0]=0;
	$sex_chr_countz[1]=0;
	$sex_chr_countz[2]=0;
	foreach $num (@mothers) {
	    if ($num == 0) { # the mother genotype is WW
	        $genotype_countz[0]+=1;
			$sex_chr_countz[0]+=2;
	    }
		elsif($num == 1) { # the mother genotype is WZ
        	$genotype_countz[1]+=1;
			$sex_chr_countz[0]+=1;
			$sex_chr_countz[1]+=1;
		}	
	}	
	foreach $num (@fathers) {
	    if ($num == 2) { # the father genotype is ZZ
	        $genotype_countz[2]+=1;
			$sex_chr_countz[1]+=2;
	    }
		elsif($num == 3) { # the father genotype is ZY
        	$genotype_countz[3]+=1;
			$sex_chr_countz[1]+=1;
			$sex_chr_countz[2]+=1;
		}	
		elsif($num == 4) { # the father genotype is WY
        	$genotype_countz[4]+=1;
			$sex_chr_countz[0]+=1;
			$sex_chr_countz[2]+=1;
		}	
	}	
	# now print to outfile
	# generation, genotype counts, sex chr counts
	if((grep(0, @genotype_countz))||(grep(0, @sex_chr_countz))){
		print "The end\n";
		last;
	}
	else{
		print OUTFILE1 $a," @genotype_countz @sex_chr_countz ", ($genotype_countz[0]+$genotype_countz[1])/($genotype_countz[2]+$genotype_countz[3]+$genotype_countz[4]),"\n";
	}
	
}

close OUTFILE1;
```
