# ABBABABA with Liberia

I want to test for gene flow between the Liberia sample and each of the other trop samples using Patterson's D. I'm going to use subgenomeL of XL as an outgroup. This means mapping the trop data to this subgenome only, and also extracting mapped reads from XL to subgenomeL, and then remapping these reads to the subgenomeL reference and then doing joint genotyping with everyone.

To extract mapped reads from XL on info (as detailed here https://darencard.net/blog/2017-09-07-extract-fastq-bam/).

First make a bam file with only subgenomeL:
```
samtools view -b -L ../XL_subgenomeL.bed SRR3210959_SRR3210971_SRR3210972_sorted.bam > SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL.bam
```
Now extract mapped reads:
```
samtools view -h -f 1 -F 12 SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL.bam > SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped.bam
```
Now sort by coordinate:
```
samtools sort SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped.bam -o SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped_sorted.bam
```
Now export fastqs (combined F and R and singletons):
```
samtools fastq -0 /dev/null SRR3210959_SRR3210971_SRR3210972_sorted_subgenomeL_mapped_name_sorted.bam > SRR3210959_SRR3210971_SRR3210972_all_reads.fq 
```
These can be aligned to a reference using `bwa mem -p ` which recognizes these combined reads:
```
bwa mem -p ../XENLA_10.1_genome_subgenomeL_only.fa SRR3210959_SRR3210971_SRR3210972_all_reads.fq | samtools view -Shu - | samtools sort - -o SRR3210959_SRR3210971_SRR3210972_extractedfq_to_subgenomeL_ref_sorted.bam
```

I mapped these XL accessions: SRR3210959_SRR3210971_SRR3210972

I genotyped and filtered these samples using GATK as previously (e.g. Evans et al. 2022) including these filters:
```
/home/ben/2025_XL_v10_Lsubgenome_ref/XT_fq/gatk-4.6.0.0/gatk --java-options -Xmx8G VariantFiltration -V XT_Lsubgenome_Chr3L_genotyped.vcf.gz\
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O XT_Lsubgenome_Chr3L_genotyped_filtered.vcf.gz
```
```
/home/ben/2025_XL_v10_Lsubgenome_ref/XT_fq/gatk-4.6.0.0/gatk --java-options -Xmx8G SelectVariants \
	        --exclude-filtered \
	        -V XT_Lsubgenome_Chr9_10L_genotyped_filtered.vcf.gz \
	        -O XT_Lsubgenome_Chr9_10L_genotyped_filtered_removed.vcf.gz
```
# General genomics

On info113, concatenate the filtered files:
```
bcftools concat XT_Lsubgenome_Chr{1..9}L_genotyped_filtered_removed.vcf.gz -Ov -o XT_Lsubgenome_allchrsgenotyped_filtered_removed.vcf
```
Then I filtered to remove positions with any missing genotypes and retain only positiojns with up to two variants:
```
vcftools --gzvcf XT_Lsubgenome_allchrsgenotyped_filtered_removed.vcf.gz --max-missing-count 0 --max-alleles 2 --recode --recode-INFO-all --out XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2

mv XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.recode.vcf XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.vcf
```
Compress (only on info113):
```
bgzip -c XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.vcf > XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.vcf.gz
```

On info2020, I converted these files to geno format like this:
```
python3 /home/ben/2025_genomics_general/genomics_general/VCF_processing/parseVCF.py -i XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 max=100 -o XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.geno.gz
```

These files still had uncalled genotypes (N/N and N|N). So, for each chromosome I removed these and then concatenated all of the chromosomes for the final analysis:

```
zcat XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.geno | sed '/N\/N/d' > XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_a.geno

cat XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_a.geno | sed '/N|N/d' > XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_b.geno
```
compress (on info113):
```
bgzip -c XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_b.geno > XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_b.geno.gz
```



# ABABABAtest

I used 5 million bp nonoverlapping windows and required at least 100 informative positions (on info2020):

```
python3 /home/ben/2025_genomics_general/genomics_general/ABBABABAwindows.py -g XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_b.geno.gz -f phased -o XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_b_ABBA_SL_IC_LB_OUT.csv --windType coordinate -w 5000000 -m 100 -s 5000000 -P1 SL -P2 IC -P3 LB -O OUT -T 10 --minData 0.5 --popsFile pops2.txt --writeFailedWindows
```
where pops2.txt is:
```
AMNH17272  SL
AMNH17274  SL
BJE4360  GW
BJE4362  GE
BJE4687  GE
EUA0335  NG
JBL052  NG
M_SierraLeone_AMNH17271  SL
M_SierraLeone_AMNH17273  SL
ROM19161  LB
SRR3210959_SRR3210971_SRR3210972  OUT
XT1  GE
XT10_WZ  GE
XT11_WW  GE
XT7_WY  GE
xen228  IC
```
