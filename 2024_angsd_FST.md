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
#SBATCH --time=2:00:00
#SBATCH --mem=32gb
#SBATCH --output=angsd_fst.%J.out
#SBATCH --error=angsd_fst.%J.err
#SBATCH --account=def-ben

# you need a text file called bam.txt that has the name of the bam file in it
# do this for each bam file in a pairwise comparison
# sbatch ../ben_scripts/2024_angsd_fst_step1.sh ref bamfilename bedfilename

# /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat_goodnamez.fasta
# /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa

#module load StdEnv/2020 angsd/0.939
module load StdEnv/2023 angsd/0.940

# angsd -bam ${1}.txt -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -rf /home/ben/projects/rrg-ben/ben/2022_Liberia/20_bams_XT_lib_mel_cal_readgroups/genome_wo_Chr7SL.bed -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -out ${1}_ -dosaf 1 -gl 1

# with bed file
# angsd -b ${2}.txt -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -rf ${3} -anc ${1} -out ${2}_ -dosaf 1 -gl 1

# without bed file
angsd -b ${2}.txt -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -anc ${1} -out ${2}_ -dosaf 1 -gl 1
```
# Genome-wide Fst Step2:
```
#!/bin/sh
#SBATCH --job-name=angsd_fst_step2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --mem=256gb
#SBATCH --output=angsd_fst_step2.%J.out
#SBATCH --error=angsd_fst_step2.%J.err
#SBATCH --account=def-ben

# you need a text file called bam.txt that has the name of the bam file in it
# do this for each bam file in a pairwise comparison
# sbatch ../ben_scripts/2024_angsd_fst_step1.sh bamfilename

#module load StdEnv/2020 angsd/0.939
module load StdEnv/2023 angsd/0.940

z_files="./*.saf.idx"
echo $z_files
for i in $z_files; do
  #echo "$i"
    for j in $z_files; do
	#echo "$j"
	if [ "$i" \> "$j" ]; then

	 echo "$i" "$j"
        realSFS $i $j > ${i:2:-9}_${j:2:-8}.mml
      fi
  done
done
```

# Genomewide Fst Step3:
```
#!/bin/sh
#SBATCH --job-name=angsd_fst_step3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --mem=256gb
#SBATCH --output=angsd_fst_step3.%J.out
#SBATCH --error=angsd_fst_step3.%J.err
#SBATCH --account=def-ben

# you need a text file called bam.txt that has the name of the bam file in it
# do this for each bam file in a pairwise comparison
# sbatch ../ben_scripts/2024_angsd_fst_step1.sh bamfilename

module load StdEnv/2020 angsd/0.939

z_files="*.saf.idx"
for i in $z_files; do
  for j in $z_files; do
     # if [ "$i" \> "$j" ]; then
      filename="${i}_${j}.ml"
      echo ${filename}
        realSFS $i $j >${filename}
     # fi
  done
done
```
# Genomewide Fst step4:
```
#!/bin/sh
#SBATCH --job-name=angsd_fst_step4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=32gb
#SBATCH --output=angsd_fst_step4.%J.out
#SBATCH --error=angsd_fst_step4.%J.err
#SBATCH --account=def-ben

# you need a text file called bam.txt that has the name of the bam file in it
# do this for each bam file in a pairwise comparison
# sbatch ../ben_scripts/2024_angsd_fst_step1.sh bamfilename

module load StdEnv/2023 angsd/0.940

#realSFS ${1} ${2} > ${1}_${2}.ml


z_files="*.saf.idx"
for i in $z_files; do
  for j in $z_files; do
     # if [ "$i" \> "$j" ]; then
        #echo "$i" "$j"
        #realSFS $i $j > $i_$j.ml
	realSFS fst index ${i} ${j} -sfs ${i}_${j}.ml -fstout ${i}_${j}_fsst
	realSFS fst stats ${i}_${j}_fsst.fst.idx -bootstrap 1000 > ${i}_${j}_boot.txt
    # fi
  done
done

#realSFS fst index ${1} ${2} -sfs ${1}_${2}.ml -fstout ${1}_${2}_fsst
#realSFS fst stats ${1}_${2}_fsst.fst.idx -bootstrap 1000 > ${1}_${2}_boot.txt
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
