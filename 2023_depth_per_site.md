# Depth per site

directory:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/raw_data/combined/bam_bad_chrnames/depth_per_site
```

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

module load StdEnv/2023  gcc/12.3 samtools/1.20

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

# Plot
```R
library (ggplot2)
library(tidyverse)
library(reshape2) # this facilitates the overlay plot
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2021_XT_sex_linked_markers/2022_WGS/depth_per_site")
# read in the data 

depth_per_site <- read.table("new_Chr7_depthpersite_0MB_to_12MB.txt", header = T)
is.data.frame(depth_per_site)
table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site)
#head(depth_per_site)

depth_per_site[, c(3:19)] <- sapply(depth_per_site[, c(3:19)], as.numeric)

# Standardize the depth by the mean depth per sample
x <- mean(na.omit(depth_per_site$F_Ghana_WZ_BJE4687));x
# get rid of outliers
depth_per_site$F_Ghana_WZ_BJE4687 <- sapply(depth_per_site$F_Ghana_WZ_BJE4687,
                                            function(r) ifelse( (is.numeric(r))&&(r<x*3),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$F_Ghana_WZ_BJE4687));x
# standardize
depth_per_site$F_Ghana_WZ_BJE4687_standardized <- sapply(depth_per_site$F_Ghana_WZ_BJE4687,
                                                         function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$M_Nigeria_EUA0335));x
# get rid of outliers
depth_per_site$M_Nigeria_EUA0335 <- sapply(depth_per_site$M_Nigeria_EUA0335, 
                                            function(r) ifelse( (is.numeric(r))&&(r<x*3),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$M_Nigeria_EUA0335));x
# standardize
depth_per_site$M_Nigeria_EUA0335_standardized <- sapply(depth_per_site$M_Nigeria_EUA0335,
                                                        function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$F_IvoryCoast_xen228));x
# get rid of outliers
depth_per_site$F_IvoryCoast_xen228 <- sapply(depth_per_site$F_IvoryCoast_xen228, 
                                             function(r) ifelse( (is.numeric(r))&&(r<x*3),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$F_IvoryCoast_xen228));x
# standardize
depth_per_site$F_IvoryCoast_xen228_standardized <- sapply(depth_per_site$F_IvoryCoast_xen228,
                                                          function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$M_SierraLeone_AMNH17271));x
# get rid of outliers
depth_per_site$M_SierraLeone_AMNH17271 <- sapply(depth_per_site$M_SierraLeone_AMNH17271, 
                                                 function(r) ifelse( (is.numeric(r))&&(r<x*3),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$M_SierraLeone_AMNH17271));x
# standardize
depth_per_site$M_SierraLeone_AMNH17271_standardized <- sapply(depth_per_site$M_SierraLeone_AMNH17271, 
                                                              function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$F_Nigeria_EUA0331));x
# get rid of outliers
depth_per_site$F_Nigeria_EUA0331 <- sapply(depth_per_site$F_Nigeria_EUA0331, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$F_Nigeria_EUA0331));x
# standardize
depth_per_site$F_Nigeria_EUA0331_standardized <- sapply(depth_per_site$F_Nigeria_EUA0331, 
                                                        function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$M_SierraLeone_AMNH17273));x
# get rid of outliers
depth_per_site$M_SierraLeone_AMNH17273 <- sapply(depth_per_site$M_SierraLeone_AMNH17273, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$M_SierraLeone_AMNH17273));x
# standardize
depth_per_site$M_SierraLeone_AMNH17273_standardized <- sapply(depth_per_site$M_SierraLeone_AMNH17273, 
                                                              function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$F_Nigeria_EUA0333));x
# get rid of outliers
depth_per_site$F_Nigeria_EUA0333 <- sapply(depth_per_site$F_Nigeria_EUA0333, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$F_Nigeria_EUA0333));x
# standardize
depth_per_site$F_Nigeria_EUA0333_standardized <- sapply(depth_per_site$F_Nigeria_EUA0333, 
                                                        function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$ROM19161));x
# get rid of outliers
depth_per_site$ROM19161 <- sapply(depth_per_site$ROM19161, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$ROM19161));x
# standardize
depth_per_site$ROM19161_standardized <- sapply(depth_per_site$ROM19161, 
                                               function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$F_SierraLeone_AMNH17272));x
# get rid of outliers
depth_per_site$F_SierraLeone_AMNH17272 <- sapply(depth_per_site$F_SierraLeone_AMNH17272, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$F_SierraLeone_AMNH17272));x
# standardize
depth_per_site$F_SierraLeone_AMNH17272_standardized <- sapply(depth_per_site$F_SierraLeone_AMNH17272, 
                                                              function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$XT10_WZ));x
# get rid of outliers
depth_per_site$XT10_WZ <- sapply(depth_per_site$XT10_WZ, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$XT10_WZ));x
# standardize
depth_per_site$XT10_WZ_standardized <- sapply(depth_per_site$XT10_WZ, 
                                              function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$F_SierraLeone_AMNH17274));x
# get rid of outliers
depth_per_site$F_SierraLeone_AMNH17274 <- sapply(depth_per_site$F_SierraLeone_AMNH17274, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$F_SierraLeone_AMNH17274));x
# standardize
depth_per_site$F_SierraLeone_AMNH17274_standardized <- sapply(depth_per_site$F_SierraLeone_AMNH17274, 
                                                              function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$XT11_WW));x
# get rid of outliers
depth_per_site$XT11_WW <- sapply(depth_per_site$XT11_WW, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$XT11_WW));x
# standardize
depth_per_site$XT11_WW_standardized <- sapply(depth_per_site$XT11_WW, 
                                              function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$M_Ghana_WY_BJE4362));x
# get rid of outliers
depth_per_site$M_Ghana_WY_BJE4362 <- sapply(depth_per_site$M_Ghana_WY_BJE4362, 
                                 function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$M_Ghana_WY_BJE4362));x
# standardize
depth_per_site$M_Ghana_WY_BJE4362_standardized <- sapply(depth_per_site$M_Ghana_WY_BJE4362, 
                                                         function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$XT1_ZY));x
# get rid of outliers
depth_per_site$XT1_ZY <- sapply(depth_per_site$XT1_ZY, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$XT1_ZY));x
# standardize
depth_per_site$XT1_ZY_standardized <- sapply(depth_per_site$XT1_ZY, 
                                             function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$M_Ghana_ZY_BJE4360));x
# get rid of outliers
depth_per_site$M_Ghana_ZY_BJE4360 <- sapply(depth_per_site$M_Ghana_ZY_BJE4360, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$M_Ghana_ZY_BJE4360));x
# standardize
depth_per_site$M_Ghana_ZY_BJE4360_standardized <- sapply(depth_per_site$M_Ghana_ZY_BJE4360, 
                                                         function(r) ifelse ((is.numeric(r)),r/x,NA))


x <- mean(na.omit(depth_per_site$XT7_WY));x
# get rid of outliers
depth_per_site$XT7_WY <- sapply(depth_per_site$XT7_WY, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$XT7_WY));x
# standardize
depth_per_site$XT7_WY_standardized <- sapply(depth_per_site$XT7_WY, 
                                             function(r) ifelse ((is.numeric(r)),r/x,NA))

x <- mean(na.omit(depth_per_site$M_Nigeria_EUA0334));x
# get rid of outliers
depth_per_site$M_Nigeria_EUA0334 <- sapply(depth_per_site$M_Nigeria_EUA0334, 
                                           function(q) ifelse( (is.numeric(q))&&(q<x*3),q,NA))
# recalculate mean
x <- mean(na.omit(depth_per_site$M_Nigeria_EUA0334));x
# standardize
depth_per_site$M_Nigeria_EUA0334_standardized <- sapply(depth_per_site$M_Nigeria_EUA0334, 
                                                        function(r) ifelse ((is.numeric(r)),r/x,NA))



# calculate the F-M difference in mean depth 
depth_per_site$mean_diff_F_minus_M_notads <- rowMeans(depth_per_site[ , c(20,22,24,26,27,28,30)], na.rm=TRUE) -
                            rowMeans(depth_per_site[ , c(21,23,25,32,34,36)], na.rm=TRUE)



png(filename = "depthdiff_per_site_8.58_8.62Mb.png",w=1800, h=800,units = "px", bg="transparent")
p<-ggplot(depth_per_site, aes(x=POS, y=mean_diff_F_minus_M_notads)) + 
  geom_line() +
  scale_x_continuous(name="Position on Chr7", limits = c(8580000,8620000)) +
  scale_y_continuous(name="(F-M) Mean Standardized Depth Difference") + #, limits = c(-50,50)) +
  theme_classic()
p
dev.off()

# these are positions that are deleted in females
head(depth_per_site[order(depth_per_site$mean_diff_F_minus_M_notads),c("POS","mean_diff_F_minus_M_notads") ],n=20)
# these are positions that are deleted in males
head(depth_per_site[order(-depth_per_site$mean_diff_F_minus_M_notads),c("POS","mean_diff_F_minus_M_notads") ],n=20)


# get rid of the first column which is the Chr name
new <- depth_per_site[,c(-1,-20:-37)];new

# make new df with only the standardized depths
new <- depth_per_site[,c(-1,-3:-19)];new
# make a new table with column headers as a factor 
df <- gather(new, "sample","depth", -POS);head(df)
dim(df)
unique(df$sample)

# make a new colum with color for each sex
df$color <- "red"
df$color[df$sample =="M_Nigeria_EUA0335_standardized"] <- "blue"
df$color[df$sample =="M_SierraLeone_AMNH17271_standardized"] <- "blue"
df$color[df$sample =="M_SierraLeone_AMNH17273_standardized"] <- "blue"
df$color[df$sample =="M_Ghana_WY_BJE4362_standardized"] <- "blue"
df$color[df$sample =="XT1_ZY_standardized"] <- "blue"
df$color[df$sample =="M_Ghana_ZY_BJE4360_standardized"] <- "blue"
df$color[df$sample =="XT7_WY_standardized"] <- "blue"
df$color[df$sample =="M_Nigeria_EUA0334_standardized"] <- "blue"

df$sample <- factor(df$sample,
                    levels = c("XT11_WW_standardized","XT10_WZ_standardized",
                               "ROM19161_standardized","F_IvoryCoast_xen228_standardized","F_Ghana_WZ_BJE4687_standardized",
                               "F_Nigeria_EUA0331_standardized",
                               "F_Nigeria_EUA0333_standardized",
                               "F_SierraLeone_AMNH17272_standardized","F_SierraLeone_AMNH17274_standardized",
                                "M_Nigeria_EUA0335_standardized","M_SierraLeone_AMNH17271_standardized",
                               "M_SierraLeone_AMNH17273_standardized","M_Ghana_WY_BJE4362_standardized",
                               "M_Ghana_ZY_BJE4360_standardized","M_Nigeria_EUA0334_standardized",
                               "XT1_ZY_standardized","XT7_WY_standardized"), ordered = T)

png(filename = "depth_per_site_standardized.png",w=1200, h=2400,units = "px", bg="transparent")
p<-ggplot(df %>% arrange(sample), aes(x=as.numeric(POS)/1000000, y=as.numeric(depth), color = color)) + 
  geom_line() +
  scale_x_continuous(name="Position on Chr7 (Mb)") + #, limits = c(8.047746,8.054998)) +
  scale_y_continuous(name="Depth") + #, limits = c(0,100)) +
  facet_wrap(~ sample, scales = "free_y", ncol = 1, strip.position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside")
p
dev.off()

newwww <-df %>% arrange(color);View(newwww)





png(filename = "log_depth_per_site_ratio_6_to_8Mb.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(depth_per_site, aes(x=position, y=ratio)) + 
  #geom_point(aes(color = factor(variable))) +
  #scale_color_manual(breaks = c("WW_minus_WZ", "WW_minus_WY", "WZ_minus_WY"), values=c("red", "blue", "gray")) +
  geom_point(size=0.25, alpha = 0.2) +
  #geom_line(aes(colour = "red"), linetype = 1) +
  #scale_y_continuous(name="Coverage ratio", limits=c(-5,5)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr7", limits=c(6000000, 8000000)) +
  # get rid of gray background
  theme_bw() 
p
dev.off()





png(filename = "log_depth_per_site_ratio_8_to_10Mb_.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(depth_per_site, aes(x=position, y=ratio)) + 
  #geom_point(aes(color = factor(variable))) +
  #scale_color_manual(breaks = c("WW_minus_WZ", "WW_minus_WY", "WZ_minus_WY"), values=c("red", "blue", "gray")) +
  geom_point(size=0.25, alpha = 0.2) +
  #geom_line(aes(colour = "red"), linetype = 1) +
  #scale_y_continuous(name="Coverage ratio", limits=c(-5,5)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr7", limits=c(8000000, 10000000),breaks=seq(8000000,10000000,100000)) +
  # get rid of gray background
  theme_bw() 
p
dev.off()




# here is a function to calculate moving averages
# https://stackoverflow.com/questions/743812/calculating-moving-average
moving_fun <- function(x, w, FUN, ...) {
  # x: a double vector
  # w: the length of the window, i.e., the section of the vector selected to apply FUN
  # FUN: a function that takes a vector and return a summarize value, e.g., mean, sum, etc.
  # Given a double type vector apply a FUN over a moving window from left to the right, 
  #    when a window boundary is not a legal section, i.e. lower_bound and i (upper bound) 
  #    are not contained in the length of the vector, return a NA_real_
  if (w < 1) {
    stop("The length of the window 'w' must be greater than 0")
  }
  output <- x
  for (i in 1:length(x)) {
    # plus 1 because the index is inclusive with the upper_bound 'i'
    lower_bound <- i - w + 1
    if (lower_bound < 1) {
      output[i] <- NA_real_
    } else {
      output[i] <- FUN(x[lower_bound:i, ...])
    }
  }
  output
}


# example
# v <- seq(1:10)

# compute a MA(2)
# moving_fun(v, 2, mean)


# make three new columns
depth_per_site_subset$WW_minus_WZ <- depth_per_site_subset$XT11_WW - depth_per_site_subset$XT10_WZ
depth_per_site_subset$WW_minus_WY <- depth_per_site_subset$XT11_WW - depth_per_site_subset$XT7_WY
depth_per_site_subset$WZ_minus_WY <- depth_per_site_subset$XT10_WZ - depth_per_site_subset$XT7_WY

# calculate the moving average
WW_minus_WZ_moving_average <- moving_fun(depth_per_site_subset[,6], 1000, mean)
WW_minus_WY_moving_average <- moving_fun(depth_per_site_subset[,7], 1000, mean)
WZ_minus_WY_moving_average <- moving_fun(depth_per_site_subset[,8], 1000, mean)

depth_per_site_subset <- cbind(depth_per_site_subset,
                               WW_minus_WZ_moving_average,
                               WW_minus_WY_moving_average,
                               WZ_minus_WY_moving_average)

# this creates a new dataframe for the overlay plot
# the column names will be "position", "variable", and "value"
# and "variable" will be the grouping variable in the overlay plot
# the values of "variable" will be "WW_minus_WZ","WW_minus_WY", and "WZ_minus_WY"
overlay_plot_Data <- melt(depth_per_site_subset[,c(2,9:11)], id="position")
head(overlay_plot_Data)
# subset the data so plotting is quicker

# overlay_plot_Data_subset <- overlay_plot_Data[
#          (overlay_plot_Data$position>10000000 & 
#             overlay_plot_Data$position<12000000),]
# make a plot 
# Open a png file

png(filename = "depth_per_site.png",w=1200, h=800,units = "px", bg="transparent")
  p<-ggplot(overlay_plot_Data, aes(x=position, y=value, col = variable)) + 
    #geom_point(aes(color = factor(variable))) +
    #scale_color_manual(breaks = c("WW_minus_WZ", "WW_minus_WY", "WZ_minus_WY"), values=c("red", "blue", "gray")) +
    geom_point(size=0.25, alpha = 0.2) +
    #geom_line(aes(colour = "red"), linetype = 1) +
    scale_y_continuous(name="Coverage difference", limits=c(-100,100)) +
    # log transform y-axis
    scale_x_continuous(name="Position on Chr7", limits=c(10000000, 12000000)) +
    # get rid of gray background
    theme_bw() 
    # Get rid of the legend
    #theme(legend.position = "none")
  p + facet_wrap( ~ variable, nrow=3)
dev.off()
  
```
