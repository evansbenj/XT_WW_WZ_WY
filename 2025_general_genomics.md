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
