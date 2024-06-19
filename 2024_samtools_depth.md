# Getting depth of coverage for mapped reads

Script is here:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/ben_scripts
```

```
#!/bin/sh
#SBATCH --job-name=samtools_depth
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=samtools_depth.%J.out
#SBATCH --error=samtools_depth.%J.err
#SBATCH --account=rrg-ben

module load StdEnv/2023  gcc/12.3 samtools/1.20
samtools coverage -q 20 --plot-depth ${1} -o ${1}_coverage.txt
```
