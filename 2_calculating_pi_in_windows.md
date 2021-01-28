# Calculating pi in windows for each genome

I'm going to try two approaches - vcftools and general_genomics. The former works from a vcf file, which I now have.  The latter works from a geno file which I am making like this (2020_make_geno_from_vcf.sh):
```
#!/bin/sh
#SBATCH --job-name=makegeno
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=2gb
#SBATCH --output=makegeno.%J.out
#SBATCH --error=makegeno.%J.err
#SBATCH --account=def-ben

# sbatch 2020_make_geno_from_vcf.sh path_and_name_of_vcf.gz_file

python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dep
th_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i ${1} 
| gzip > ${1}.geno.gz
```
