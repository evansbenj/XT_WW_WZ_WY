# Mel to Germany assembly

I am trying to figure out why Chr8 has a signal of sex linkage in mel. I mapped the data to the Germany assembly and uysed angsd to analyze sex linkage.

Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2017_mello_GBS/bamz_mapped_to_germany_mello
```
# Figure out which contigs map to chr8 of XT
To do this I mapped each of the Germany assembly contigs to the XT genome using minimap2

Then I extracted the ones that hit Chr8:
```
grep '	Chr8	' XT_to_germany_mel_alignments.paf > XT_to_germany_mel_alignments_Chr8_hitz_only.paf
```

Then I extracted the ones with a match length of at least 1000:
```
awk '$10>999' XT_to_germany_mel_alignments_Chr8_hitz_only.paf > XT_to_germany_mel_alignments_Chr8_hitz_only_matching_gt_1000.paf
```

Here is a script to pull out angsd positions that are on mel contigs that minimap2 mapped to chr8:
```
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a angsd output file and one other file
# that is a list of contigs from the Germany genome that minimap2
# says map a particular XT chromosome (Chr8 in this case)

# execute like this:
# ./Makes_a_mel_Chr8_from_angsdoutput.pl out_Xlaev_additive_F1.lrt0.gz mell_Chr8_contigs.txt mel_Chr8_contigs_angsdout

# where mell_Chr8_contigs.txt are contigs that have at least 1000 bp of matching seqs to XT
# awk '$10>999' XT_to_germany_mel_alignments_Chr8_hitz_only.paf > XT_to_germany_mel_alignments_Chr8_hitz_only_matching_gt_1000.paf

# and mel_Chr8_contigs_angsdout is the output files withg associations for only this chr

# we need to accept any positions on these contigs (and not match specific variants)

my $angsd = $ARGV[0];
my $contigz = $ARGV[1];
my $outputfile = $ARGV[2];
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




#### first open the matsites and patsites and load this into a hash
my @contigz;


unless (open DATAINPUT, $contigz) {
	print "Can not find the matsites file, jackass.\n";
	exit;
}
my $counter=0;
while ( my $line = <DATAINPUT>) {
	chomp($line);
	#@columns=split(/\s+/,$line);
	$contigz[$counter] = $line;
	$counter+=1;
}	
close DATAINPUT;


# now load the angsd output 
if ($angsd =~ /.gz$/) {
	open(DATAINPUT2, "gunzip -c $angsd |") || die "can’t open pipe to $angsd";
}
else {
	open DATAINPUT2, $angsd or die "Could not read from $angsd: $!";
}

$counter=0;
my @data;

while ( my $line = <DATAINPUT2>) {
	@columns=split(/\s+/,$line);
	foreach(@contigz){
		if(($columns[0] ne 'Chromosome') && ($columns[0] == $_)){
			$data[$counter] = $line;
			$counter += 1;
		}
	}	
} # end while
close DATAINPUT2;


# now print the data
foreach(@data){
	print OUTFILE $_;
}
close OUTFILE;
```
