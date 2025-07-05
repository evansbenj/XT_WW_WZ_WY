# Exploring depth on Chr8L

I suspect that the signal on Chr8L is due to some males having atypically high depth and these positions are then being better called as pseudoheterozygous. TO test this I am quantifying depth at the 76 positions that are het in the father and at lest five sons and comparing this to the chromosome-wide depth.

Make bed file for Chr8 and Chr7:
```
Chr7:1-133565930	1	133565930
```
```
Chr8:1-147241510	1	147241509
```

Directory on graham:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2017_mello_GBS/mapped_to_XTv10_concatscaf
```
Get depth per site:
```
samtools depth -b XT_Chr7.bed male_BJE4181_sorted.bam_rg.bam > male_BJE4181_sorted.bam_rg.bam_Chr7.depth
```

Get average depth:
```
awk '{ sum_values += $3 } END { print sum_values/NR }' male_BJE4181_sorted.bam_rg.bam_Chr7.depth
```

Now get depth of the 76 positions with paternal heterozygosity:
```
egrep '(str1|str2|str3)'
```
