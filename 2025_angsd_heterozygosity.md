# Heterozygosity with angsd

directory:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/20_bams_XT_lib_mel_cal_readgroups
```
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
thetaStat do_stat ${1}_out.thetas.idx -win 50000 -step 50000  -outnames ${1}_theta.thetasWindow.gz
```

# Plotting
```R
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2021_XT_sex_linked_markers/2025_pi_windows")
library(ggplot2)
library(readr)
library(dplyr)
dir <- "/Users/Shared/Previously Relocated Items/Security/projects/2021_XT_sex_linked_markers/2025_pi_windows"
list.files(dir)

# load the data 

AMNH17271 <- read.table("AMNH17271_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(AMNH17271) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
AMNH17271$sample <- "SL1"
AMNH17271$color <- "steelblue"

AMNH17272 <- read.table("AMNH17272_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(AMNH17272) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
AMNH17272$sample <- "SL2"
AMNH17272$color <- "steelblue"

AMNH17273 <- read.table("AMNH17273_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(AMNH17273) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
AMNH17273$sample <- "SL3"
AMNH17273$color <- "steelblue"

AMNH17274 <- read.table("AMNH17274_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(AMNH17274) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
AMNH17274$sample <- "SL4"
AMNH17274$color <- "steelblue"

lib <- read.table("lib_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(lib) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
lib$sample <- "LB"
lib$color <- "orange"

ic <- read.table("ic_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(ic) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
ic$sample <- "IC"
ic$color <- "purple"



jbl <- read.table("jbl_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(jbl) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
jbl$sample <- "NG1"
jbl$color <- "green"

EUA0331 <- read.table("EUA0331_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(EUA0331) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
EUA0331$sample <- "NG2"
EUA0331$color <- "green"

EUA0333 <- read.table("EUA0333_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(EUA0333) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
EUA0333$sample <- "NG3"
EUA0333$color <- "green"

EUA0334 <- read.table("EUA0334_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(EUA0334) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
EUA0334$sample <- "NG4"
EUA0334$color <- "green"

EUA0335 <- read.table("EUA0335_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(EUA0335) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
EUA0335$sample <- "NG5"
EUA0335$color <- "green"

BJE4360 <- read.table("BJE4360_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(BJE4360) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
BJE4360$sample <- "GW"
BJE4360$color <- "red"

BJE4362 <- read.table("BJE4362_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(BJE4362) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
BJE4362$sample <- "GE1"
BJE4362$color <- "red"

BJE4687 <- read.table("BJE4687_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(BJE4687) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
BJE4687$sample <- "GE2"
BJE4687$color <- "red"


XT1 <- read.table("XT1_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(XT1) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
XT1$sample <- "GE3"
XT1$color <- "red"

XT7 <- read.table("XT7_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(XT7) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
XT7$sample <- "GE4"
XT7$color <- "red"

XT10 <- read.table("XT10_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(XT10) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
XT10$sample <- "GE5"
XT10$color <- "red"

XT11 <- read.table("XT11_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(XT11) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
XT11$sample <- "GE6"
XT11$color <- "red"


cal <- read.table("cal_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(cal) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
cal$sample <- "Cal"
cal$color <- "gray"

mel <- read.table("mel_theta.thetasWindow.gz.pestPG", header = T, skip = 1)
colnames(mel) <- c("details","Chr","midpoint","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
mel$sample <- "Mel"
mel$color <- "black"




my_df <- rbind(AMNH17271,AMNH17272,AMNH17273,AMNH17274,
               lib,ic,BJE4360,BJE4362,BJE4687,XT1,XT7,XT10,XT11,jbl,EUA0331,EUA0333,EUA0334,EUA0335,cal,mel)
#my_df <- rbind(SL1,SL2,SL3,SL4,LB,IC,GW,GE1,GE2,GE3,GE4,GE5,GE6,NG1,NG2,NG3,NG4,NG5,Cal,Mel)

my_df$pi_per_site <- my_df$tP/my_df$nSites


#my_df$color <- "red"
#my_df$color[my_df$color == "mal" ] <- "blue"


#my_df$sample <- factor(my_df$sample, levels = c('AMNH17271','AMNH17272','AMNH17273','AMNH17274',
#                            'lib','ic','BJE4360','BJE4362','BJE4687',
#                            'XT1','XT7','XT10','XT11','jbl','EUA0335','EUA0331',
#                            'EUA0333','EUA0334','cal','mel'), order = T)

my_df$sample <- factor(my_df$sample, levels = c('SL1','SL2','SL3','SL4',
                            'LB','IC','GW','GE1','GE2',
                            'GE3','GE4','GE5','GE6','NG1','NG2','NG3',
                            'NG4','NG5','Cal','Mel'), order = T)

library(Hmisc)
pdf("./XT_pi_boxplot.pdf",w=8, h=2.0, version="1.4", bg="transparent")
  p<-ggplot(my_df %>% arrange(sample), aes(x=sample, y=pi_per_site, fill = color, color = color)) + 
    geom_violin(width = 1.5) +
    stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
                 colour = "black", size=0.15) +
    #geom_boxplot() +
    # color the stuff the way I want
    scale_fill_manual(values = c("black" = "black", "gray"="gray", "orange" = "orange",
                                  "red" = "red","green" = "green", "purple" = "purple",
                                  "steelblue" = "steelblue")
                      )+
    scale_color_manual(values = c("black" = "black", "gray"="gray", "orange" = "orange",
                                 "red" = "red","green" = "green", "purple" = "purple",
                                 "steelblue" = "steelblue")
    )+
    xlab("Sample") + ylab(expression(paste(italic(pi)))) +
    # get rid of gray background
    theme_classic() +
    theme(legend.position="none") #+
    p
dev.off()
```
