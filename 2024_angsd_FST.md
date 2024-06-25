# Using angsd for Fst

An advantage of using angsd for Fst is that it works from bam files. We'd like to plot pairwise Fst between each trop individual but excluding the sex-linked region (the first 20Mb of Chr7).

bam files are here:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_bams_XT_lib_mel_cal_readgroups
```
Discussion here about using angsd for pairwise Fst between individuals:
```
https://github.com/ANGSD/angsd/issues/413
```


# Step 1
```
#!/bin/sh
#SBATCH --job-name=angsd_fst_step1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
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
# Calculate the 2dsfs prior 
From: https://www.popgen.dk/angsd/index.php/Fst
```
../misc/realSFS pop1.saf.idx pop2.saf.idx >pop1.pop2.ml
```

# prepare the fst for easy window analysis etc
```
../misc/realSFS fst index pop1.saf.idx pop2.saf.idx -sfs pop1.pop2.ml -fstout here
```
...or on graham:
```
#!/bin/sh
#SBATCH --job-name=angsd_fst_step2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=32gb
#SBATCH --output=angsd_fst_step2.%J.out
#SBATCH --error=angsd_fst_step2.%J.err
#SBATCH --account=rrg-ben

# you need a text file called bam.txt that has the name of the bam file in it
# do this for each bam file in a pairwise comparison
# sbatch ../ben_scripts/2024_angsd_fst_step1.sh bamfilename

module load StdEnv/2023 angsd/0.940

#realSFS temp1.saf.idx temp2.saf.idx > temp1_temp2.ml
#realSFS fst index temp1.saf.idx temp2.saf.idx -sfs temp1_temp2.ml -fstout pop1_pop2_fsst
realSFS fst index ${1} ${2} -sfs ${3} -fstout ${1}_${2}_fsst
realSFS fst stats ${1}_${2}_fsst.fst.idx -bootstrap 1000 > ${1}_${2}_boot.txt
```
# Get the global estimate
```
../misc/realSFS fst stats here.fst.idx
```
-> FST.Unweight:0.069395 Fst.Weight:0.042349

# Below is not tested that much, but seems to work
```
../misc/realSFS fst stats2 here.fst.idx -win 50000 -step 10000 >slidingwindow
```
