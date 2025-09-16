# PCA and Admixture

I am doing the PCA and admixture analysis again starting with only trop samples (excluding the hybrid Nigeria samples). I first will work with a vcf that includes only trop samples (including EUA335):

On info in this directory:
```
/home/ben/2025_XL_v10_Lsubgenome_ref/XT_data_only
```
Extract the XT samples:
```
bcftools view -S trop_only.txt combined_Chr9.g.vcf.gz_Chr9_GenotypedSNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz -Oz -o troponly_Chr9.vcf.gz
```
Now merge the chromosomes:
```
bcftools concat troponly_Chr{1..10}.vcf.gz -Oz -o troponly_allChrs.vcf.gz
```


For PCA and Admixture I am going to filter sites with missing data and then thin every 5000 bp:
```
vcftools --gzvcf trop_only_allchrs_concat.vcf.gz --max-missing-count 0 --minQ 30 --recode --recode-INFO-all --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30
```
```
vcftools --vcf trop_only_allchrs_concat_maxmissingcount_0_genoqual30.recode.vcf --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.vcf --thin 5000 --recode
```

For the PCA of the XT data mapped to the XLsubgenomeL, I also thinned to include only every 500bp (on info):
directory:
```
/home/ben/2025_XL_v10_Lsubgenome_ref/XT_subgenomeL
```
```
vcftools --gzvcf XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.vcf.gz --minQ 30 --out XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2_thin_5000.vcf --thin 5000 --recode
```

for subgenome L the same file used for PCA is this one:
```
/home/ben/2025_XL_v10_Lsubgenome_ref/individual_vcfz/troponly_forreal_XT_Lsubgenome_allchrsgenotyped_filtered_removed_nomissing_2.vcf.gz_thin_5000.vcf.recode.vcf.gz
```

# Admixture analysis

on info I concatenated  chrs:
```
bcftools concat all_162_maqs_chr{1..20}_maxmissingcount_0_genoqual30_thin_5000.recode.vcf -Ov -o all_162_maqs_allautsomal_chrs_maxmissingcount_0_genoqual30_thin_5000.recode.vcf
```
and then I used plink to make the input files:
```
plink --gzvcf trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.vcf --make-bed --out trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000 --allow-extra-chr
```

ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
```
awk '{$1="0";print $0}' trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim > trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim.tmp
mv trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim.tmp trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bim
```
I then opened 20 screens and did independent analyses (with different seeds) in separate folders:
```
for i in {2..7}
do
/usr/local/admixture/admixture --cv --seed $((1 + $RANDOM % 1000)) ../trop_only_allchrs_concat_maxmissingcount_0_genoqual30_thin_5000.bed $i > log${i}.out
done
```
# Rename files
I used a perl script to rename files for Admixture plotter:
```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This script will rename directories for admixtureplotter

# Execute it in a directory that has 2 folders in it called
# admix1...admix20


my $y= 7; # this is the number of ancestral populations
my $x;
my $z;
my $command;


#my @directorynamez = ("admix1/","admix2/","admix3/","admix4/","admix5/","admix6/",
#						"admix7/","admix8/","admix9/","admix10/",
#						"admix11/","admix12/","admix13/","admix14/","admix15/","admix16/",
#						"admix17/","admix18/","admix19/","admix20/",);

# make directories
for ($x = 2 ; $x <= $y ; $x++ ) {
	$command = "mkdir ".$x;
	system($command);
	$command = "mkdir ".$x."/Logs";
	system($command);
	for ($z = 1 ; $z <= 20 ; $z++ ) {
		$command = "mkdir ".$x."/".$z;
		system($command);
	}	
}


for ($z = 1 ; $z <= 20 ; $z++ ) {
	for ($x = 2 ; $x <= $y ; $x++ ) {
		$command = "mv admix".$z."/*".$x.".P ".$x."/".$z."/.";
		system($command);
		$command = "mv admix".$z."/*".$x.".Q ".$x."/".$z."/.";
		system($command);
		$command = "mv admix".$z."/log".$x.".out ".$x."/Logs/".$x."_".$z.".log";
		system($command);
	}	
}
```


# AdmixturePlotter (https://github.com/TCLamnidis/AdmixturePlotter)
You need to download the `.bed` file that was used to run admixture and put it in the same directory as the renamed output files

Then run this command (after copying these scripts and files to the same directory:
```
poporder.txt
allsamples.ind
colorlist.txt
CompileData.sh
```

```
./CompileData.sh
```
```
./AdmixturePlotter.R -i ../Plotting/compound.labelled.QperK.txt -p ./poporder.txt -c ./colorlist.txt
```
