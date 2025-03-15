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
python3 /home/ben/2025_genomics_general/genomics_general/popgenWindows.py -w 5000000 -m 100 -g trop_only_allchrs_concat.vcf.gz.geno_b.geno.gz -o trop_only_allchrs_concat.vcf.gz.geno_b.geno_diversity.csv.gz -f phased -T 5 -p GH1 -p IC1 -p SL1 -p SL2 -p NG1 -p GH2 -p GH3 -p SL3 -p SL4 -p GH4 -p GH5 -p GH6 -p GH7 -p LB1 --popsFile pops1.txt
```
