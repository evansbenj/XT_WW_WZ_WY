# General genomics
(on info)

I made an input file using the trop data only:
```
python3 /home/ben/2025_genomics_general/genomics_general/VCF_processing/parseVCF.py -i trop_only_allchrs_concat.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 max=100 -o trop_only_allchrs_concat.vcf.gz.geno.gz
```
Then I removed positions with Ns:
```
zcat trop_only_allchrs_concat.vcf.gz.geno.gz | sed '/N\/N/d' > trop_only_allchrs_concat.vcf.gz.geno_a
cat trop_only_allchrs_concat.vcf.gz.geno_a | sed '/N|N/d' > trop_only_allchrs_concat.vcf.gz.geno_b.geno
```
Now calculate Fst (on info2020)
```
python3 /home/ben/2025_genomics_general/genomics_general/popgenWindows.py -w 5000000 -m 100 -g trop_only_allchrs_concat.vcf.gz.geno_b.geno -o trop_only_allchrs_concat.vcf.gz.geno_b.geno_diversity.csv.gz -f phased -T 5 -p F_Ghana_WZ_BJE4687_combined__sorted.bam -p F_IvoryCoast_xen228_combined__sorted.bam -p F_SierraLeone_AMNH17272_combined__sorted.bam -p F_SierraLeone_AMNH17274_combined__sorted.bam -p JBL052_concatscafs_sorted.bam -p M_Ghana_WY_BJE4362_combined__sorted.bam -p M_Ghana_ZY_BJE4360_combined__sorted.bam -p M_Nigeria_EUA0335_combined__sorted.bam  M_SierraLeone_AMNH17271_combined__sorted.bam -p M_SierraLeone_AMNH17273_combined__sorted.bam -p XT10_WZ_no_adapt._sorted.bam -p XT11_WW_trim_no_adapt_scafconcat_sorted.bam -p XT1_ZY_no_adapt._sorted.bam -p XT7_WY_no_adapt__sorted.bam -p all_ROM19161_sorted.bam --popsFile pops1.txt
```
