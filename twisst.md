# Twisst
Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/ben_scripts
```
first make nj trees

```
python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/phylo/phyml_sliding_windows.py -T 10 -g ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz --prefix ../raw_data/XTgenomez_Chr7.vcf.gz_SNPsonly_first20mil_XT11nohet.vcf.recode.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz_treez_w50 -w 50 --windType sites --model GTR
```

# but this won't work without an outgroup...
