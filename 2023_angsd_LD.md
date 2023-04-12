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

# Memory issues, run by sites
A trial run that was filtering and then doing the association test failed with an out of memory error with 128 Gb assigned.  I'm rerunning now with 256.

In the meantime I am trying to use the sites command to run the association on specific sites that were output from the failed run. This is only a subset of the genome so it will need to be repeated.

First get two columns from output and then sort by chr and then index using angsd:
```
module load angsd
zcat out_longer.lrt0.gz > log # remove last line which was incomplete using emacs
cut -f1,2 log > sites.temp
# remove the header before sorting using emacs or whatever
sort -k1 sites.txt >sorted_sites.txt 
angsd sites index sorted_sites.txt 
```
now one can run the association test on these particular sites:
```
angsd -yBin bin_sex.ybin -doAsso 1 -doMaf 1 -doMajorMinor 1 -GL 1 -sites sorted_sites.txt -out tempty -Pvalue 1 -bam bam.filelist 
```
 the `-Pvalue 1` causes the output to be a Pvalue instead of a chi square value with df = 1

# could try:
 `-model  2` with males coded as 1 in the ybin file (0 being the controls, 1 being the cases)
 default is `-model 1` which is additive



* first step is to filter bam files (http://popgen.dk/angsd/index.php/Filters#Allele_frequencies)


# This worked

Filter to remove non-significant sites for plotting.
Without the "-Pvalue 1' flag, the 6th column is a chisq value with df=1, so lets save only significant values (higher than 7):
```
zcat tempty.lrt0.gz | awk '$6 < 7 { next } { print }'> sig.only
```
