# General genomics
(on info)

I made an input file using the trop data only:
```

```
Then I removed positions with Ns:
```
zcat trop_only_allchrs_concat.vcf.gz.geno.gz | sed '/N\/N/d' > trop_only_allchrs_concat.vcf.gz.geno_a
cat trop_only_allchrs_concat.vcf.gz.geno_a | sed '/N|N/d' > trop_only_allchrs_concat.vcf.gz.geno_b.geno
```
