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
le::-15}_trim.R1.fq.gz ${file::-15}_trim.R1_single.fq.gz ${file::-15}_trim.R2.fq.gz ${file::-15}_trim.R2_single.fq.gz SLID
INGWINDOW:4:15 MINLEN:36 HEADCROP:3
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-15}_WW_R1.fastq.gz ${file::-15}_WW_R2.fastq.gz ${fil
e::-15}_trim.R1.fq.gz ${file::-15}_trim.R1_single.fq.gz ${file::-15}_trim.R2.fq.gz ${file::-15}_trim.R2_single.fq.gz SLIDI
NGWINDOW:4:15 MINLEN:36 HEADCROP:3
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
# sbatch 2020_align_paired_fq_to_ref.sh ~/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta.gz pathtofqfi
lez WW

module load bwa/0.7.17
module load samtools/1.10


for file in ${2}/*${3}*.R1.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	bwa mem ${1} ${file::-9}.R1.fq.gz ${file::-9}.R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${file::-
9}_sorted.bam
	samtools index ${file::-9}_sorted.bam
  fi
done
```


