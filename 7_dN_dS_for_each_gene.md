# First get the coordinates of each gene

I wrote a perl script that does this:
```
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in gff file and outputs
# the coordinates of all CDS for each gene

# to execute type ./Get_coordinates_of_CDS_in_each_gene.pl XENTR_10.0_Xenbase_Chr7_lt_30Mb.gff3 output1 


my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

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
		if($temp[6] eq '+'){
			$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[3]=1;	# forward orientation
		}
		elsif($temp[6] eq '-'){
			$gene_hash{$Gene_name."_".$Gene_ID}{$Transcript_ID}{$exon_counter}[3]=-1;	# reverse orientation
		}
		$exon_counter+=1;
	}
} # end while	
close DATAINPUT;	
# OK, now all the CDS are in a hash

foreach my $key (keys %gene_hash){
	foreach my $transcript_id (keys %{$gene_hash{$key}}){
		foreach my $exon (keys %{$gene_hash{$key}{$transcript_id}}){
			print $key," ",$transcript_id," ",$gene_hash{$key}{$transcript_id}{$exon}[0]," ",$gene_hash{$key}{$transcript_id}{$exon}[1]," ",$gene_hash{$key}{$transcript_id}{$exon}[2]," ";
			print $gene_hash{$key}{$transcript_id}{$exon}[1] - $gene_hash{$key}{$transcript_id}{$exon}[0]," ";
			print $gene_hash{$key}{$transcript_id}{$exon}[3],"\n";
		}
	}		
}	
close OUTFILE;
```
