# Mel to Germany assembly

I am trying to figure out why Chr8 has a signal of sex linkage in mel. I mapped the data to the Germany assembly and uysed angsd to analyze sex linkage.

Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2017_mello_GBS/bamz_mapped_to_germany_mello
```
# Figure out which contigs map to chr8 of XT
To do this I mapped each of the Germany assembly contigs to the XT genome using minimap2
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/ben_scripts/2025_minimap_XT_to_mel_output_paf.sh
```
```
#!/bin/sh
#SBATCH --job-name=minmap
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=minmap.%J.out
#SBATCH --error=minmap.%J.err
#SBATCH --account=def-ben

# minimap2 -x asm10 -a --secondary=no -t8 reference.fasta query.fasta >alignments.sam
module load StdEnv/2020 minimap2/2.24

#/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/bin/minimap2/
minimap2 -x asm10 --secondary=no -t8 /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/Germany_genome/Super_NovaXeno_mega_gt200.fasta > XT_to_germany_mel_alignments.paf
```

Then I extracted the ones that hit Chr7 and Chr8:
```
grep '	Chr7	' XT_to_germany_mel_alignments.paf > XT_to_germany_mel_alignments_Chr7_hitz_only.paf
grep '	Chr8	' XT_to_germany_mel_alignments.paf > XT_to_germany_mel_alignments_Chr8_hitz_only.paf
```

Then I extracted the ones within Chr8 with a match length of at least 1000:
```
awk '$10>999' XT_to_germany_mel_alignments_Chr8_hitz_only.paf > XT_to_germany_mel_alignments_Chr8_hitz_only_matching_gt_1000.paf
awk '$10>999' XT_to_germany_mel_alignments_Chr7_hitz_only.paf > XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000.paf
```

Then I required a map quality of at least 60
```
awk '$12>59' XT_to_germany_mel_alignments_Chr8_hitz_only_matching_gt_1000.paf > XT_to_germany_mel_alignments_Chr8_hitz_only_matching_gt_1000_mq60.paf
awk '$12>59' XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000.paf > XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000_mq60.paf
```

I also am going to extract the ones that match Chr7 above and below 20Mb:
```
awk '$9<20000000' XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000_mq60.paf > XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000_mq60_lt_20Mb.paf
awk '$9>20000000' XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000_mq60.paf > XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000_mq60_gt_20Mb.paf
```

Here is a script to pull out angsd positions that are on mel contigs that minimap2 mapped to chr8 (or whatever). One of the input files is a list of chromosomes from the paf files above:
```
cut -f1 XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000_gt_20Mb.paf > XT_to_germany_mel_alignments_Chr7_hitz_only_matching_gt_1000_gt_20Mb_contigs.txt
```
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
