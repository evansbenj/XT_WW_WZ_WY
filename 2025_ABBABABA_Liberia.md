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

On info, I converted these files to geno format like this:
```
python3 /home/ben/2025_genomics_general/genomics_general/VCF_processing/parseVCF.py -i XT_Lsubgenome_Chr9_10L_genotyped_filtered_removed.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 max=100 -o XT_Lsubgenome_Chr9_10L_genotyped_filtered_removed.geno.gz
```
ABABABAtest:
```
python3 /home/ben/2025_genomics_general/genomics_general/ABBABABAwindows.py -g XT_Lsubgenome_Chr1L_genotyped_filtered_removed.geno.gz -f phased -o XT_Lsubgenome_Chr1L_ABBA.csv --windType coordinate -w 100000 -m 100 -s 100000 -P1 SL1 -P2 IC1 -P3 LIB -O OUT -T 10 --minData 0.5 --popsFile pops.txt --writeFailedWindows
```
where pops.txt is:
```
AMNH17272 SL1
AMNH17274 SL2
BJE4360 GH1
BJE4362 GH2
BJE4687 GH3
JBL052 NG1
M_SierraLeone_AMNH17271 SL3
M_SierraLeone_AMNH17273 SL4
ROM19161 LIB
SRR3210959_SRR3210971_SRR3210972 OUT
XT1 LB1
XT10_WZ LB2
XT11_WW LB3
XT7_WY LB4
xen228 IC1
```
