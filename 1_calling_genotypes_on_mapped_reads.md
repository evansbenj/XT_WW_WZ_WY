# Calling genotypes on mapped reads

# First add readgroups using picard:
```
#!/bin/sh
#SBATCH --job-name=picard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=2gb
#SBATCH --output=picard.%J.out
#SBATCH --error=picard.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_addreadgroups.sh


module load StdEnv/2020 picard/2.23.3

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=../raw_data/XT10_WZ_trim_sorted.bam O=../raw_data/XT10_WZ_trim_sorted_rg.bam RGID=XT10_WZ RGLB=XT10_WZ RGPL=illum
ina RGPU=XT10_WZ RGSM=XT10_WZ
```
(and the same for XT11_WW and XT7_WY)

# Index the new bam files
```
#!/bin/sh
#SBATCH --job-name=picard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=2gb
#SBATCH --output=picard.%J.out
#SBATCH --error=picard.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch 2020_samtools_index_bam.sh pathtobamfile/bamfile

module load StdEnv/2020 samtools/1.11

samtools index ${1}
```

# Realign indels using GATK
```
#!/bin/sh
#SBATCH --job-name=gatk_indelrealigner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=gatk_indelrealigner.%J.out
#SBATCH --error=gatk_indelrealigner.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch 2020_GATK_indelrealigner.sh pathtoref/ref  pathtobamfile/bamfile
# sbatch 2020_GATK_indelrealigner.sh ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta bamfile

module --force purge
module load nixpkgs/16.09
module load gatk/3.8

java -Xmx2g  -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${1} -I ${2} -o ${2}.intervals

java -Xmx2g  -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${1} -I ${2} -targetIntervals ${2}.intervals -o ${2}_realigned.bam


module load StdEnv/2020 samtools/1.11
samtools index ${2}_realigned.bam

```

# Call genotypes
```
#!/bin/sh
#SBATCH --job-name=haplotypecaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=haplotypecaller.%J.out
#SBATCH --error=haplotypecaller.%J.err
#SBATCH --account=def-ben

# sbatch 2020_chrX_0_gatk_HaplotypeCaller.sh pathtobamfile/bamfilname

module --force purge
module load nixpkgs/16.09
module load gatk/3.8


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta.gz -I ${1} --output-mode EMIT_ALL_CONFI
DENT_SITES --emitRefConfidence GVCF -o ${1}_noBSQR.g.vcf.gz

```

# Combining gvcfs:
```
#!/bin/sh
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32000M
#SBATCH --output=GenotypeGVCFs.%J.out
#SBATCH --error=GenotypeGVCFs.%J.err
#SBATCH --account=def-ben

module load nixpkgs/16.09
module load gatk/3.8


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /home/ben/scratch/2020_XT_WW_WZ_WY/2020_XT_v10_refgen
ome/XENTR_10.0_genome.fasta -V XT10_WZ_trim_sorted_rg.bam_realigned.bam_${1}_noBSQR.g.vcf.gz -V XT11_WW_trim_sorted_
rg.bam_realigned.bam_${1}_noBSQR.g.vcf.gz -V XT7_WY_trim_sorted_rg.bam_realigned.bam_${1}_noBSQR.g.vcf.gz -allSites 
-o XTgenomez_${1}.vcf.gz
```

# Filter
```
#!/bin/sh
#SBATCH --job-name=VarFilt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=VarFilt.%J.out
#SBATCH --error=VarFilt.%J.err
#SBATCH --account=def-ben

module load nixpkgs/16.09
module load gatk/3.8


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10
.0_genome.fasta -V ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR.vcf.gz --filterExpression "QD < 2.0 || FS > 30.0 |
| MQ < 40.0 || ReadPosRankSum < -8.0 || MQRankSum < -10.00" --filterName "lowqual" --genotypeFilterExpression "DP < 5 ||
 DP > 100" --genotypeFilterName "genotypefilter" --setFilteredGtToNocall -o ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_
noBSQR_flagged.vcf.gz
```
```
#!/bin/sh
#SBATCH --job-name=indels
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --output=indels.%J.out
#SBATCH --error=indels.%J.err
#SBATCH --account=def-ben

module load nixpkgs/16.09
module load gatk/3.8


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_
genome.fasta -V ../genotypez/XT_XT11_WW_XT10_WZ_XT7_WY_Chr7_noBSQR_flagged.vcf.gz -o ../genotypez/XT_XT11_WW_XT10_WZ_XT7
_WY_Chr7_noBSQR_filtered.vcf.gz
```

# Below not used
Phased using Beagle in advance of general_genomics:
```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=def-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dept
h_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/beagle.18May20.d20.jar gt=../raw_data/XTgenomez_${1}.vcf.g
z out=../raw_data/XTgenomez_${1}_phased.vcf.gz impute=true
```
