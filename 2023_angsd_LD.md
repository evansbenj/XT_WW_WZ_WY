# Angsd LD analysis
path:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/raw_data/combined/angsd_LD
```
Before I start, I'd like to convert the names of the chromosomes in the Xt genome to match the ones in my concatenated scaffold file. For some reason when I map the reads to the concat scaf file, the resulting bam is really small - not sure why.  Better to work with the previously mapped data...

First subset the big chrs
```
samtools view -b -L chrs_only.bed JBL052_XTv10_ref_individ_scaffolds.bam > subsetJBL052_XTv10_ref_individ_scaffolds.bam_bigchrsonly.bam
```
Also subset the concatenated scaffolds from the other file I have:
```
samtools view -b -L Scafs.bed JBL052__sorted.bam > subsetJBL052__sorted.bam_Scafsonly.bam
```

# Trial association:
```
angsd -yBin bin_sex.ybin -doAsso 1 -GL 1 -out out -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -bam bam.filelist -P 5
```

# could try:
 `-model  2` with males coded as 1 in the ybin file (0 being the controls, 1 being the cases)
 default is `-model 1` which is additive



* first step is to filter bam files (http://popgen.dk/angsd/index.php/Filters#Allele_frequencies)
