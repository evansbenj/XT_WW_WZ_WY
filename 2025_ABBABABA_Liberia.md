# ABBABABA with Liberia

I want to test for gene flow between the Liberia sample and each of the other trop samples using Patterson's D. I'm going to use subgenomeL of XL as an outgroup. This means mapping the trop data to this subgenome only, and also extracting mapped reads from XL to subgenomeL, and then remapping these reads to the subgenomeL reference and then doing joint genotyping with everyone.

To extract mapped reads from XL on info (as detailed here https://darencard.net/blog/2017-09-07-extract-fastq-bam/).

First make a bam file with only subgenomeL:
```
samtools view -b -L ../XL_subgenomeL.bed SRR3210959_SRR3210971_SRR3210972_sorted.bam > SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL.bam
```
Now extract mapped reads, sort, and export fastqs:
```
samtools view -f 1 -F 12 SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL.bam > SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped.bam
```
or, I actually ended up using this because the above line created a weird bam file that seemed corrupt:
```
view -b -F 4 SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL.bam > mapped.bam
```
samtools sort mapped.bam -o mapped_sorted.bam
/usr/local-centos6/bedtools/2.19.1/bin/bamToFastq -i mapped_sorted.bam -fq SRR3210959_SRR3210971_SRR3210972_subgenomeL_mapped.1.fastq -fq2 SRR3210959_SRR3210971_SRR3210972_subgenomeL_mapped.2.fastq
```

I mapped these XL accessions: SRR3210959_SRR3210971_SRR3210972
