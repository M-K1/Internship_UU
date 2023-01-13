#!/usr/bin/Rscript


#Load libraries
library(cowplot)
library(facetscales)
library(ggplot2)
library(dplyr)
library(ggdark)
library(ggpubr)

#Import data

kmer_mappings_fx10545 <- read.csv("D:/matthijs/WUR/Internship/Output/AnthocyaninContent/fx10.545/data_for_manhattan_kmer.csv")


#Convert chromosome names to actual chromosome number
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022518.1'] <- 1
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022519.1'] <- 2
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022520.1'] <- 3
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022521.1'] <- 4
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022522.1'] <- 5
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022523.1'] <- 6
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022524.1'] <- 7
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022525.1'] <- 8
kmer_mappings_fx10545[kmer_mappings_fx10545 == 'CM022526.1'] <- 9

kmer_mappings_fx10545[kmer_mappings_fx10545 == 'NC_056623.1'] <- 1

kmer_mappings_fx10545$P <- -log10(as.numeric(kmer_mappings_fx10545$P))

#only kmers in chromosome 5
kmer_mappings_fx10545_chr5 <- subset(kmer_mappings_fx10545, CHR==5)

#add column with the positions in Mbp
kmer_mappings_fx10545_chr5$BP_Mbp <- as.numeric(kmer_mappings_fx10545_chr5$BP)/1000000

#Set intervals for first plot
minmax_kmer_pos_fx10545 <- range(kmer_mappings_fx10545_chr5$BP)
minmax_kmer_pos_fx10545_rounded <- round(minmax_kmer_pos_fx10545, -3)
intervals_fx10545 <- seq(from = minmax_kmer_pos_fx10545_rounded[1]-100000,to = (minmax_kmer_pos_fx10545_rounded[2]+100000), by = 100000)
intervals_fx10545 <- intervals_fx10545/1000000

hist(kmer_mappings_fx10545_chr5$BP_Mbp,
     main = "Histogram of fx10545 (bin size 100000)",
     xlab = "Position on chromosome 5 (Mbp)",
     breaks = intervals_fx10545)

#Using ggplot
hist_fx10545 <- ggplot(kmer_mappings_fx10545_chr5, aes(x=BP_Mbp)) +
  geom_histogram(binwidth = 0.01, fill="black") +
  # ggtitle("Histogram of fx10545, chromosome 5, 84 to 94 Mbp (bin size 10000)") +
  ylab("count") + xlab("Position on chromosome 5 (Mbp)") +
  # geom_vline(xintercept=6.080001, linetype='dotted', col = 'red',size=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(84, 94)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3000)) +
  # expand_limits(y = 6000) +
  theme_bw()

#histogram with new intervals (smaller than 100 mbp and larger than 80 mbp)
kmer_mappings_fx10545_chr5_100mbp <- subset(kmer_mappings_fx10545_chr5, BP_Mbp < 93 & BP_Mbp > 84)
intervals_fx10545 <- seq(from = 83000000-100000,to = 93000000+100000, by = 10000)
intervals_fx10545 <- intervals_fx10545/1000000

hist(kmer_mappings_fx10545_chr5_100mbp$BP_Mbp,
     main = "Histogram of fx10545 close up (bin size 10000)",
     xlab = "Position on chromosome 5 (Mbp)",
     breaks = intervals_fx10545)

#even closer
kmer_mappings_fx10545_chr5_closer <- subset(kmer_mappings_fx10545_chr5, BP_Mbp < 86.3 & BP_Mbp > 85.7)
intervals_fx10545 <- seq(from = 85800000-100000,to = 86200000+100000, by = 10000)
intervals_fx10545 <- intervals_fx10545/1000000

hist(kmer_mappings_fx10545_chr5_closer$BP_Mbp,
     main = "Histogram of fx10545 close up (bin size 10000)",
     xlab = "Position on chromosome 5 (Mbp)",
     breaks = intervals_fx10545)

#jitter plot baby
jitter_fx10_545 <- ggplot(kmer_mappings_fx10545_chr5_closer,aes(BP_Mbp,P))+
  geom_jitter(height = 0.5)+
  xlab("Position (Mbp)") + ylab("-log10(p)")+
  scale_x_continuous(expand = c(0.01, 0)) +
  expand_limits(y = 17) +
  # ggtitle("Jitter plot of significant kmers between 85.8 and 86.2 Mbp of fx10.545 on chromosome 5") +
  theme_bw() 
  
mean(kmer_mappings_fx10545_chr5_closer[["P"]])

#Other peak
kmer_mappings_fx10545_chr5_other <- subset(kmer_mappings_fx10545_chr5, BP_Mbp < 94 & BP_Mbp > 92)
jitter_fx10_545_peak_2 <- ggplot(kmer_mappings_fx10545_chr5_other,aes(BP_Mbp,P))+
  geom_jitter(height = 0.5)+
  xlab("Position (Mbp)") + ylab("-log10(p)")+
  scale_x_continuous(expand = c(0.01, 0)) +
  expand_limits(y = 17) +
  theme_bw()

mean(kmer_mappings_fx10545_chr5_other[["P"]])

range(kmer_mappings_fx10545_chr5_closer$BP_Mbp)

point_fx10_545 <- ggplot(kmer_mappings_fx10545_chr5_closer,aes(BP_Mbp,P))+
  geom_point()+
  xlab("Position (Mbp)") + ylab("-log10(p)")+
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_vline(xintercept=85.956110, linetype='dotted', col = 'red',size=1) +
  # ggtitle("Jitter plot of significant kmers between 85.8 and 86.2 Mbp of fx10.545 on chromosome 5") +
  theme_bw()

#Combine plots into one plot using gpubr

fx10545 <- ggarrange(mhpl, hist_fx10545, ggarrange(jitter_fx10_545, jitter_fx10_545_peak_2, labels = c("C", "D"), ncol = 2 ), nrow = 3, labels = c("A", "B"), heights = 1)
fx10545
