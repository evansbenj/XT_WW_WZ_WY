# Extract reads that don't map to reference:
```
```
# Assemble them

XT7_WY__1_unaligned_seqs.fq and XT7_WY__2_unaligned_seqs.fq are paired reads that do not match the ref

XT7_WY__M_unaligned_seqs.fq includes singleton forward and reverse reads that do not match the ref

First replace the '/2' of XT7_WY__M_unaligned_seqs.fq with '/1'
```
sed -i 's/\/2/\/1/g' XT7_WY__M_unaligned_seqs.fq
```
Now cat the forward paired reads with all the singletons
```
cat XT7_WY__1_unaligned_seqs.fq XT7_WY__M_unaligned_seqs.fq > XT7_WY_unmappd_paired_and_singletons.fq
```
Now assemble them
```
#!/bin/sh
#SBATCH --job-name=trinity
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --output=trinity.%J.out
#SBATCH --error=trinity.%J.err
#SBATCH --account=rrg-ben

# sbatch 2020_trinity.sh ../unaligned_seqs/XT7_WY/XT7_WY_unmappd_paired_and_singletons.fq ../unaligned_seqs/XT7_WY/XT7_WY__2_unaligned_seqs.fq ../unaligned_seqs/XT7_WY/trinity_output

module load nixpkgs/16.09  gcc/7.3.0 nixpkgs/16.09
module load openmpi/3.1.2
module load samtools/1.9
module load salmon/0.11.3
module load bowtie2/2.3.4.3
module load jellyfish/2.2.6
module load trinity/2.8.4

# get total memory from the assigned node
avail_mem=$(free -h| grep Mem| tr -s ' '| cut -d ' ' -f2)

Trinity --seqType fq --left ${1} --right ${2} --CPU 32 --full_cleanup --max_memory "${avail_mem}" --min_kmer_cov 2 --include_supertranscripts --bflyCalculateCPU  --bflyCPU 10 --bflyHeapSpaceMax 12G --output ${3}```
