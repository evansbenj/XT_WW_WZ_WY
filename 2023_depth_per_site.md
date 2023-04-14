# Depth per site

Get depth from all bam files into one file:
```
#!/bin/sh
#SBATCH --job-name=samtools_depthpersite
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=8:00:00
#SBATCH --mem=2gb
#SBATCH --output=samtools_depthpersite.%J.out
#SBATCH --error=samtools_depthpersite.%J.err
#SBATCH --account=def-ben

# sbatch 2023_samtools_depth_from_bam.sh bamfile

module load samtools/1.10

samtools depth -b Chr7.bed -f ./list_of_bams.txt -o Chr7_depthpersite.txt
```
where Chr7.bed is this:
```
Chr7:1-133565930	1	133565930
```
and list_of_bams.txt is this:
```
F_Ghana_WZ_BJE4687_combined__sorted.bam_rg.bam
M_Nigeria_EUA0335_combined__sorted.bam_rg.bam
F_IvoryCoast_xen228_combined__sorted.bam_rg.bam
M_SierraLeone_AMNH17271_combined__sorted.bam_rg.bam
F_Nigeria_EUA0331_combined__sorted.bam_rg.bam
M_SierraLeone_AMNH17273_combined__sorted.bam_rg.bam
F_Nigeria_EUA0333_combined__sorted.bam_rg.bam
ROM19161__sorted.bam_rg.bam
F_SierraLeone_AMNH17272_combined__sorted.bam_rg.bam
XT10_WZ_no_adapt._sorted.bam_rg.bam
F_SierraLeone_AMNH17274_combined__sorted.bam_rg.bam
XT11_WW_trim_no_adapt_scafconcat_sorted.bam_rg.bam
M_Ghana_WY_BJE4362_combined__sorted.bam_rg.bam
XT1_ZY_no_adapt._sorted.bam_rg.bam
M_Ghana_ZY_BJE4360_combined__sorted.bam_rg.bam
XT7_WY_no_adapt__sorted.bam_rg.bam
M_Nigeria_EUA0334_combined__sorted.bam_rg.bam
```

# Subsample the depth per site file
Using `awk`:
```
awk '$2 < 8000000 { next } { print }' new_Chr7_depthpersite.txt | awk '$2 > 9000000 { next } { print }' > Chr7_depthpersite_8MB_to_9MB.txt
```
