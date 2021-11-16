# Sequencing
NovaSeq 6000 PE150 for all four samples.  3 were run at Genome Quebec and the fourth (XT1_ZY) was run at TCAG

# Working directory
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY
```
# Trimmomatic
```
#!/bin/sh
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=8gb
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_trimmomatic.sh ../raw_data/plate1


module load StdEnv/2020
module load trimmomatic/0.39


#v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*_WW_R*.fastq.gz; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
  	#if [[ $v -eq 1 ]]
      #then # if/then branch
    echo java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-15}_WW_R1.fastq.gz ${file::-15}_WW_R2.fastq.gz ${fi
le::-15}_trim.R1.fq.gz ${file::-15}_trim.R1_single.fq.gz ${file::-15}_trim.R2.fq.gz ${file::-15}_trim.R2_single.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-15}_WW_R1.fastq.gz ${file::-15}_WW_R2.fastq.gz ${fil
e::-15}_trim.R1.fq.gz ${file::-15}_trim.R1_single.fq.gz ${file::-15}_trim.R2.fq.gz ${file::-15}_trim.R2_single.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	#	  v=0
	#else # else branch
  	#	v=1
	#fi
  fi
done 

```

# Mapping reads to XT v10
```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=48:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez sexchr_genotype
# sbatch 2020_align_paired_fq_to_ref.sh ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta.gz ../raw_data/X
T10_WZ_trim_noadapters WZ

module load bwa/0.7.17
module load samtools/1.10


for file in ${2}/*${3}*.R1.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	bwa mem ${1} ${file::-9}.R1.fq.gz ${file::-9}.R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${file::-9}_
sorted.bam
	samtools index ${file::-9}_sorted.bam
  fi
done
```


# Indel realign
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
# sbatch 2020_GATK_indelrealigner.sh pathtoref/ref  pathtobamfile/bamfile_prefix
# sbatch 2020_GATK_indelrealigner.sh ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta bamfile

# first add readgroups with picard
module load StdEnv/2020 picard/2.23.3

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      I=${2}.bam \
      O=${2}_rg.bam \
      RGID=1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=${2}

# index the new bam file
module load StdEnv/2020 samtools/1.12
samtools index ${2}_rg.bam

# now do indel realignment with GATK
module --force purge
module load nixpkgs/16.09
module load gatk/3.8

java -Xmx2g  -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${1} -I ${2}_rg.bam -o ${2}.intervals

java -Xmx2g  -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${1} -I ${2}_rg.bam -targetIntervals ${2}.intervals -
o ${2}_rg_realigned.bam

module load StdEnv/2020 samtools/1.12
# now index the realigned file
samtools index ${2}_rg_realigned.bam
```

# Dedup
```
#!/bin/sh
#SBATCH --job-name=picard_dedup
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=picard_dedup.%J.out
#SBATCH --error=picard_dedup.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch 2020_GATK_indelrealigner.sh pathtobamfile/bamfile_prefix
# sbatch 2020_GATK_indelrealigner.sh bamfile

# first add readgroups with picard
module load StdEnv/2020 picard/2.23.3

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
     REMOVE_DUPLICATES=true \
     I=${1}.bam \
     O=${1}_dedup.bam \
     M=marked_dup_metrics.txt
```
# Haplotypecaller
```
#!/bin/sh
#SBATCH --job-name=haplotypecaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=128:00:00
#SBATCH --mem=32gb
#SBATCH --output=haplotypecaller.%J.out
#SBATCH --error=haplotypecaller.%J.err
#SBATCH --account=def-ben

# sbatch 2020_0_gatk_HaplotypeCaller.sh pathtobamfile/bamfilname chr

module load samtools/1.12
samtools index ${1}.bam

module --force purge
module load nixpkgs/16.09
module load gatk/3.8


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_geno
me.fasta -I ${1}.bam -L ${2} --output_mode EMIT_ALL_CONFIDENT_SITES --emitRefConfidence GVCF -o ${1}_+${2}_noBSQR.g.vcf.gz
```
# GenotypeGVCFs
```
XXX
```

# Below not used

# Genomic regions with no coverage

```
#!/bin/sh
#SBATCH --job-name=bedtools
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=2gb
#SBATCH --output=bedtools.%J.out
#SBATCH --error=bedtools.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_regions_with_no_coverage.sh bamfile prefix
# sbatch 2020_regions_with_no_coverage.sh ../raw_data/XT11_WW_trim_sorted.bam XT_WW

module load StdEnv/2020 bedtools/2.29.2

bedtools genomecov -ibam ${1} -bg | awk '$4 < 1' > ${2}_intervals_with_no_coverage.txt
```


