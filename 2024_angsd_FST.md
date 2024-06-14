# Using angsd for Fst

An advantage of using angsd for Fst is that it works from bam files. We'd like to plot pairwise Fst between each trop individual but excluding the sex-linked region (the first 20Mb of Chr7).

bam files are here:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_bams_XT_lib_mel_cal_readgroups
```

# step 1
```
#!/bin/sh
#SBATCH --job-name=angsd_fst_step1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=32gb
#SBATCH --output=angsd_fst.%J.out
#SBATCH --error=angsd_fst.%J.err
#SBATCH --account=def-ben

# you need a text file called bam.txt that has the name of the bam file in it
# do this for each bam file in a pairwise comparison
# sbatch ../ben_scripts/2024_angsd_fst_step1.sh bamfilename

module load StdEnv/2020 angsd/0.939

# angsd -bam ${1}.txt -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -rf /home/ben/projects/rrg-ben/ben/2022_Liberia/20_bams_XT_lib_mel_cal_readgroups/genome_wo_Chr7SL.bed -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -out ${1}_ -dosaf 1 -gl 1

angsd -b ${1}.txt -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -rf /home/ben/projects/rrg-ben/ben/2022_Liberia/20_bams_XT_lib_mel_cal_readgroups/genome_wo_Chr7SL.bed -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -out ${1}b_ -dosaf 1 -gl 1
```
