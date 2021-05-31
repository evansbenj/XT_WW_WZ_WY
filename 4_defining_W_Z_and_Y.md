# Defining W, Z, and Y

The three tads we have are siblings with the following genotypes: WW, WY, and WZ. The mother of the tads was the daughter of the father of the tads. The father was WY. The mother was WZ. So the mother should have the same W chromosome as the father.  If there is recombination on the distal tips of the W and Z and/or the distal tips of the W and the Y, then the WW tad might have heterozygosity in this region.

# Workflow

First step is to evaluate the genotypes that are mapped to v10.  Probably the easiest way to do this is to first export them to a tab-delimited file and then write a perl script to generate a synthetic W, Z, and Y chromosome.  This would be interesting for estimation of genetic distances (prediction is that Z and Y are more similar to each other than either is to the W) and also to explore where differences are distributed (prediction is that differences between Z and Y should be highest immediately surrounding the sex-linked region; differences between the Z/Y and the W should be more broadly distributed).

The coordinates of these chromosomes will not match the reference due to insertion deletion events. For this reason, it probably would be wise to generate paml files directly from the tab delimited file for each gene. This would allow me to combine exons for each gene, check for and quantify frameshifts, and then use paml to quantify NS and S changes.

I think I will begin by exporting a tab-delimited file and then using a Perl script to identify sites in the WW individual and subtracting this variation to generate a Z or Y in the WZ and WY individuals. I'd have to develop some rules for resolving heterozygosity in the W. There almost certainly will be heterozygous positions due to repetitive regions.  Probably the most conservative way to deal with this is to insert Ns for anything with W heterozygosity.  

I could also use the Z- and Y-specific kmer assemblies to verify SNPs. I could generate two types of data: one with each chr based on comparisons only between the called genotypes from the WW, WZ, and WY individuals and the other based on this plus verification with the Z- and Y-specific kmer assemblies.

To accomplish the second goal of analyzing NS and S variation in alleles on each chr, I could export the coding regions of each gene (with exons concatenated) as a vcf file, convert to tab, and then use the same script to output paml formatted files for analysis.

Things to check on
* does paml accept Ns
* can I export multiple regions into a concatenated vcf file
