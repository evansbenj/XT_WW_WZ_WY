# Nuclear / WGS phylogeny 

Directory:
```
/home/ben/2025_XL_v10_Lsubgenome_ref/XT_data_only
```

We are going to make a ML phylogeny from the mapped data.

Using the hardfiltered vcf as input, filter again to retain only positions that have no missing data and that have biallelic sites:
```
vcftools --gzvcf trop_only_allchrs_concat.vcf.gz --max-missing-count 0 --min-alleles 2 --max-alleles 2 --minQ 30 --recode --recode-INFO-all --stdout | gzip -c > trop_only_allchrs_concat_maxmissingcount_0_biallelic_genoqual30.vcf.gz
```

# Model selection
```
/home/ben/IQtree/iqtree-2.4.0-Linux-intel/bin/iqtree2 -s trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.nexus -m MF
```

# Bootstrap tree
```
/home/ben/IQtree/iqtree-2.4.0-Linux-intel/bin/iqtree2 -s trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.nexus -m TVM+F+R3 -b 1000
```
