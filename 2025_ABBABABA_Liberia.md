# ABBABABA with Liberia

I want to test for gene flow between the Liberia sample and each of the other trop samples using Patterson's D. I'm going to use subgenomeL of XL as an outgroup. This means mapping the trop data to this subgenome only, and also extracting mapped reads from XL to subgenomeL, and then remapping these reads to the subgenomeL reference and then doing joint genotyping with everyone.

To extract mapped reads from XL (as detailed here https://darencard.net/blog/2017-09-07-extract-fastq-bam/):
```
samtools view -u -f 1 -F 12 lib_002.sorted.md.bam > lib_002_map_map.bam
samtools flagstat lib_002.sorted.md.bam
bamToFastq -i lib_002_mapped.sort.bam -fq lib_002_mapped.1.fastq -fq2 lib_002_mapped.2.fastq
```

I mapped these XL accessions: SRR3210959_SRR3210971_SRR3210972
