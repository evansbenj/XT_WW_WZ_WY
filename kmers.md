# Kmers with meryl

Working in this directory on graham:
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/ben_scripts
```

I installed meryl from github (https://github.com/marbl/meryl), but first had to load StdEnv/2020

```
module load StdEnv/2020
```

On graham, meryl is here:
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin
```
I made meryl kmer databases like this:

```
#!/bin/sh
#SBATCH --job-name=meryl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl.%J.out
#SBATCH --error=meryl.%J.err
#SBATCH --account=def-ben


# sbatch 2020_meryl_make_kmerdb.sh fastqfile
# ../raw_data/NS.1462.004.IDT_i7_102---IDT_i5_102.XT7_WY_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_102---IDT_i5_102.XT7_WY_R2.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_trim.R1_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_trim.R2_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_WZ_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_114---IDT_i5_114.XT10_WZ_R2.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_trim.R1_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_trim.R2_single.fq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_WW_R1.fastq.gz
# ../raw_data/NS.1462.004.IDT_i7_126---IDT_i5_126.XT11_WW_R2.fastq.gz


/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl count ${1} threads=4 memory=128 k=29 ou
tput ${1}_meryldb.out
```

# Make a union-sum of the kmer-dbs of the forward and reverrse reads:
```
#SBATCH --job-name=meryl                                                                                                        
#SBATCH --nodes=1                                                                                                               
#SBATCH --ntasks-per-node=1                                                                                                     
#SBATCH --time=48:00:00                                                                                                         
#SBATCH --mem=128gb                                                                                                             
#SBATCH --output=meryl.%J.out                                                                                                   
#SBATCH --error=meryl.%J.err                                                                                                    
#SBATCH --account=def-ben                                                                                                       


# sbatch 2020_meryl_union_kmer_dbs.sh db1 db2 out 

# sbatch 2020_meryl_union_kmer_dbs.sh ../raw_data/XT10_WZ_trim.R1.fq.gz_meryldb.out ../raw_data/XT10_WZ_trim.R2.fq.gz_meryldb.out\
# ../raw_data/XT10_WZ_R1R2_meryldb.out

# sbatch 2020_meryl_union_kmer_dbs.sh ../raw_data/XT11_WW_trim.R1.fq.gz_meryldb.out ../raw_data/XT11_WW_trim.R2.fq.gz_meryldb.out\
# ../raw_data/XT11_WW_R1R2_meryldb.out

# sbatch 2020_meryl_union_kmer_dbs.sh ../raw_data/XT7_WY_trim.R1.fq.gz_meryldb.out ../raw_data/XT7_WY_trim.R2.fq.gz_meryldb.out .\
# ./raw_data/XT7_WY_R1R2_meryldb.out

/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY/bin/meryl/build/bin/meryl union-sum ${1} ${2} threads=4 memory=128 k=29 output \
${3}
```
