# PCA and Admixture

I am doing the PCA and admixture analysis again starting with only trop samples (excluding the hybrid Nigeria samples). I first will work with a vcf that includes only trop samples (including EUA335):

On info in this directory:
```
/home/ben/2025_XL_v10_Lsubgenome_ref/XT_data_only
```
Extract the XT samples:
```
bcftools view -S trop_only.txt combined_Chr9.g.vcf.gz_Chr9_GenotypedSNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz -Oz -o troponly_Chr9.vcf.gz
```
Now merge the chromosomes:
```
bcftools concat troponly_Chr{1..10}.vcf.gz -Oz -o troponly_allChrs.vcf.gz
```


For PCA and Admixture I am going to filter sites with missing data and then thin every 5000 bp:
```
vcftools --gzvcf trop_only_allchrs_concat.vcf.gz --max-missing-count 0 --minQ 30 --recode --recode-INFO-all --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30
```
```
vcftools --vcf trop_only_allchrs_concat_maxmissingcount_0_genoqual30.recode.vcf --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.vcf --thin 5000 --recode
```

# Admixture analysis

on info I concatenated  chrs:
```
bcftools concat all_162_maqs_chr{1..20}_maxmissingcount_0_genoqual30_thin_5000.recode.vcf -Ov -o all_162_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.recode.vcf
```
and then I used plink to make the input files:
```
plink --gzvcf trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.vcf --make-bed --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000 --allow-extra-chr
```

ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
```
awk '{$1="0";print $0}' trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim > trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim.tmp
mv trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim.tmp trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim
```
I then opened 20 screens and did independent analyses (with different seeds) in separate folders:
```
for i in {2..7}
do
/usr/local/admixture/admixture --cv --seed $((1 + $RANDOM % 1000)) ../trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bed $i > log${i}.out
done
```

# AdmixturePlotter (https://github.com/TCLamnidis/AdmixturePlotter)

```
./AdmixturePlotter.R -i ../Plotting/compound.labelled.QperK.txt -p ../poporder.txt
```
