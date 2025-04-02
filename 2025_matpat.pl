# Making a matpat output from angsd

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
