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

# Realign indels using GATK
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
# sbatch 2020_GATK_indelrealigner.sh pathtoref/ref  pathtobamfile/bamfile
# sbatch 2020_GATK_indelrealigner.sh ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta.gz bamfile


module load nixpkgs/16.09
module load gatk/3.8

java -Xmx2g  -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${1} -I ${2} -o ${2}.intervals

java -Xmx2g  -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${1} -I ${2} -targetIntervals ${2}.intervals -o ${2}_realigned.bam


module load StdEnv/2020 samtools/1.11
samtools index ${2}_realigned.bam
```

# call genotypes
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

module load nixpkgs/16.09
module load gatk/3.8


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta.gz -I ${1} --output-mode EMIT_ALL_CONFI
DENT_SITES --emitRefConfidence GVCF -o ${1}_noBSQR.g.vcf.gz

```
