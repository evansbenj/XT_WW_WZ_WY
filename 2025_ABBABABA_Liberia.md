# ABBABABA with Liberia

I want to test for gene flow between the Liberia sample and each of the other trop samples using Patterson's D. I'm going to use subgenomeL of XL as an outgroup. This means mapping the trop data to this subgenome only, and also extracting mapped reads from XL to subgenomeL, and then remapping these reads to the subgenomeL reference and then doing joint genotyping with everyone.

To extract mapped reads from XL on info (as detailed here https://darencard.net/blog/2017-09-07-extract-fastq-bam/).

First make a bam file with only subgenomeL:
```
samtools view -b -L ../XL_subgenomeL.bed SRR3210959_SRR3210971_SRR3210972_sorted.bam > SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL.bam
```
Now extract mapped reads:
```
samtools view -h -f 1 -F 12 SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL.bam > SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped.bam
```
Now sort by coordinate:
```
samtools sort SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped.bam -o SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped_sorted.bam
```
Now export fastqs (combined F and R and singletons):
```
samtools fastq -0 /dev/null SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped_name_sorted.bam > SRR3210959_SRR3210971_SRR3210972_all_reads.fq 
```
These can be aligned to a reference using `bwa mem -p ` which recognizes these combined reads:
```
bwa mem -p ../XENLA_10.1_genome_subgenomeL_only.fa SRR3210959_SRR3210971_SRR3210972_all_reads.fq | samtools view -Shu - | samtools sort - -o SRR3210959_SRR3210971_SRR3210972_extractedfq_to_subgenomeL_ref_sorted.bam
```

I mapped these XL accessions: SRR3210959_SRR3210971_SRR3210972
