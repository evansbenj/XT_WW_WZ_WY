# Liberia and others
working directory:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/raw_data/combined
```
and
```
/home/ben/projects/rrg-ben/ben/2020_XT_WW_WZ_WY
```
and mellotrop from Germany
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/Germany_genome/raw_data
```
and XT from ref genome
```
/home/ben/projects/rrg-ben/ben/2020_XT_v10_raw_data/XT-v10_rawdata/
```

In total we have 19 WGS samples:
```
Fieldnumber	Museum_or_barcode	species	Country	Locality	Sex	Sex chromosome genotype if known
AMNH17271	2644.055	Xenopus tropicalis	Sierra Leone	near Freetown	Male	
AMNH17272	unknown	Xenopus tropicalis	Sierra Leone	near Freetown	Female	
AMNH17273	unknown	Xenopus tropicalis	Sierra Leone	near Freetown	Male	
AMNH17274	unknown	Xenopus tropicalis	Sierra Leone	near Freetown	Female	
xen228	no specimen	Xenopus tropicalis	Ivory Coast	Adiopodoume	Female	
BJE4362	3D6.1D5980178C	Xenopus tropicalis	Ghana	Ametozofe (Ghana East)	Male	WY
BJE4687	3D6.1D5980179E	Xenopus tropicalis	Ghana	Ametozofe (Ghana East)	Female	WZ
BJE4360	alive	Xenopus tropicalis	Ghana	Ankasa West	male	ZY
EUA0331	Z23682	Xenopus tropicalis	Nigeria	Benin	Female	
EUA0333	Z23684	Xenopus tropicalis	Nigeria	Benin	Female	
EUA0334	Z23685	Xenopus tropicalis	Nigeria	Benin	Male	
EUA0335	Z23686	Xenopus tropicalis	Nigeria	Benin	Male	
ROM19161	ROM19161	Xenopus tropicalis	Liberia		Female	
XTR_1_ZY	tadpole	Xenopus tropicalis	Ghana		probM	ZY
XT7_WY	tadpole	Xenopus tropicalis	Ghana		probM	WY
XT11_WW	tadpole	Xenopus tropicalis	Ghana		probF	WW
XT10_WZ	tadpole	Xenopus tropicalis	Ghana		probF	WZ
JBL052 adult Xenopus tropicalis lab F probWZ
```

# XT genome with concatenated scaffolds:
```
/home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta
```

# Align and genotyping
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts/2020_align_paired_fq_to_ref.sh
```
```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=64:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2018_Austin_XB_genome/Austin_genome/Xbo.v1.fa.gz patht
ofqfilez
# or for XL genome use:
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa.gz pathtofq
filez

module load bwa/0.7.17
module load samtools/1.10


for file in ${2}/*trim_R1.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	#${file::-9}
	echo bwa mem ${1} ${file::-13}trim_R1.fq.gz ${file::-13}trim_R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - 
-o ${file::-9}_sorted.bam
	bwa mem ${1} ${file::-13}trim_R1.fq.gz ${file::-13}trim_R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${
file::-13}_sorted.bam
	samtools index ${file::-13}_sorted.bam
  fi
done

```

# haplotype caller
```
/home/ben/projects/rrg-ben/ben/2021_M_f_aurea/ben_scripts/2021_HaplotypeCaller.sh
```
```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --output=HaplotypeCaller.%J.out
#SBATCH --error=HaplotypeCaller.%J.err
#SBATCH --account=def-ben


# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 

# execute like this:
# sbatch 2021_HaplotypeCaller.sh /home/ben/projects/rrg-ben/ben/2021_rhemac_v10/rheMac10.fa.sa path chr

module load nixpkgs/16.09 gatk/4.1.0.0

for file in ${2}*_sorted.bam_rg.bam
do
    gatk --java-options -Xmx8G HaplotypeCaller  -I ${file} -R ${1} -L ${3} -O ${file}_${3}.g.vcf -ERC GVCF
done
```
