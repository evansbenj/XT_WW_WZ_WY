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
module load StdEnv/2020 trinity/2.11.0
module load samtools jellyfish 
module load gcc/9.3.0  openmpi/4.0.3 salmon/1.3.0

Trinity --seqType fq --left XT7_WY_unmappd_paired_and_singletons.fq --right XT7_WY__2_unaligned_seqs.fq --no_normalize_reads --max_memory 10G
```
