# Angsd LD analysis
path:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/raw_data/combined/angsd_LD
```
Before I start, I'd like to convert the names of the chromosomes in the Xt genome to match the ones in my concatenated scaffold file

First subset the big chrs
```
samtools view -b -L chrs_only.bed JBL052_XTv10_ref_individ_scaffolds.bam > subsetJBL052_XTv10_ref_individ_scaffolds.bam_bigchrsonly.bam
```


* first step is to filter bam files (http://popgen.dk/angsd/index.php/Filters#Allele_frequencies)
