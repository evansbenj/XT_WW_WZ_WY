# Nuclear / WGS phylogeny 

Directory:
```
/home/ben/2025_XL_v10_Lsubgenome_ref/XT_data_only
```

We are going to make a ML phylogeny from the mapped data.

Using the hardfiltered vcf as input, convert to nexus file; require no missing data (min of 15 samples with data for each site).
```
python3 vcf2phylip/vcf2phylip.py -i trop_only_allchrs_concat.vcf.gz -m 15 -n --output-prefix trop_only_allchrs_concat.vcf.gz_nexus
```

# Model selection
```
/home/ben/IQtree/iqtree-2.4.0-Linux-intel/bin/iqtree2 -s trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.nexus -m MF
```

# Bootstrap tree
```
/home/ben/IQtree/iqtree-2.4.0-Linux-intel/bin/iqtree2 -s trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.nexus -m TVM+F+R3 -b 1000
```

# Or use thinned vcf file:
```
python3 vcf2phylip/vcf2phylip.py -i trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.vcf.gz -m 15 -n --output-prefix trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000_nexus
```
```
/home/ben/IQtree/iqtree-2.4.0-Linux-intel/bin/iqtree2 -s trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.nexus -m TVM+F+R3 -b 1000 --redo-tree
```
