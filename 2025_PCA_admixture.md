# PCA and Admixture

I am doing the PCA and admixture analysis again starting with only trop samples (excluding the hybrid Nigeria samples). For PCA and Admixture I am going to filter sites with missing data and then thin every 5000 bp:
```
vcftools --gzvcf trop_only_allchrs_concat.vcf.gz --max-missing-count 0 --minQ 30 --recode --recode-INFO-all --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30
```
```
vcftools --vcf trop_only_allchrs_concat_maxmissingcount_0_genoqual30.recode.vcf --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.vcf --thin 5000 --recode
```
