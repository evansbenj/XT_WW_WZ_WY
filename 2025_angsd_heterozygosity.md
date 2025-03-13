# Heterozygosity with angsd

First make a text file with the path to each bam. Then run this sbatch script:
```
#!/bin/sh
#SBATCH --job-name=angsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=24gb
#SBATCH --output=angsd.%J.out
#SBATCH --error=angsd.%J.err
#SBATCH --account=rrg-ben


module load StdEnv/2023 angsd/0.940

#angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta -GL 1 -out ${1}_angsd_out

# angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta -GL 1 -out ${1}_out
angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -GL 1 -out ${1}_out
realSFS ${1}_out.saf.idx -P 24 -fold 1 > ${1}_out.sfs
realSFS saf2theta ${1}_out.saf.idx -sfs ${1}_out.sfs -outname ${1}_out
thetaStat do_stat ${1}_out.thetas.idx -win 50000 -step 10000  -outnames ${1}_theta.thetasWindow.gz
```
