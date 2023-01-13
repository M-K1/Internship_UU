#!/usr/bin/Rscript


#Load libraries
library(cowplot)
library(facetscales)
library(ggplot2)
library(ggdark)
library(ggpubr)

#Import data

kmer_mappings <- read.csv("D:/matthijs/WUR/Internship/Output/introgression_6010001/data_for_manhattan_kmer.csv")

#Convert chromosome names to actual chromosome number
kmer_mappings[kmer_mappings == 'CM022518.1'] <- 1
kmer_mappings[kmer_mappings == 'CM022519.1'] <- 2
kmer_mappings[kmer_mappings == 'CM022520.1'] <- 3
kmer_mappings[kmer_mappings == 'CM022521.1'] <- 4
kmer_mappings[kmer_mappings == 'CM022522.1'] <- 5
kmer_mappings[kmer_mappings == 'CM022523.1'] <- 6
kmer_mappings[kmer_mappings == 'CM022524.1'] <- 7
kmer_mappings[kmer_mappings == 'CM022525.1'] <- 8
kmer_mappings[kmer_mappings == 'CM022526.1'] <- 9

#only kmers in chromosome 1
kmer_mappings_chr1 <- subset(kmer_mappings, CHR==1)

minmax_kmer_pos <- range(kmer_mappings_chr1$BP)
minmax_kmer_pos_rounded <- round(minmax_kmer_pos, -3)
intervals <- seq(from = minmax_kmer_pos_rounded[1]-10000,to = (minmax_kmer_pos_rounded[2]+10000), by = 10000)
kmer_mappings_chr1$BP_Mbp <- as.numeric(kmer_mappings_chr1$BP)/1000000


#Create histogram
hist(kmer_mappings_chr1$BP,
main = "Histogram of introgression 6010001 (bin size 10000)",
xlab = "Position on chromosome 1",
breaks = intervals)

#Using ggplot
hist_601 <- ggplot(kmer_mappings_chr1, aes(x=BP_Mbp)) +
  geom_histogram(binwidth = 0.01, fill="black") +
  # ggtitle("Histogram of introgression 6080001, 4 to 8 Mbp (bin size 10000)") +
  ylab("count") + xlab("Position on chromosome 1 (Mbp)") +
  geom_vline(xintercept=6.010001, linetype='dotted', col = 'red',size=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

#plot only between region of interest with kmer position and p-value
kmer_mappings_interest <- subset(kmer_mappings, BP < 5954833 & BP > 5933483)
kmer_mappings_interest$P <- -log10(as.numeric(kmer_mappings_interest$P))
kmer_mappings_interest$BP_Mbp <- as.numeric(kmer_mappings_interest$BP)/1000000

blocks_601 <- ggplot(kmer_mappings_interest,aes(BP_Mbp,P))+
  geom_point()+
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Position (Mbp)") + ylab("-log10(p)")+
  expand_limits(y = max(kmer_mappings_interest$P) + 1) +
  # ggtitle("Jitter plot of significant kmers between 5.8 and 6.2 Mbp") +
  geom_vline(xintercept=5.933483, linetype='dashed', col = 'red') +
  geom_vline( xintercept = 5.954833, col="purple", linetype="dashed")+
  geom_vline( xintercept = 5.982821, col="blue", linetype="dashed")+
  geom_vline( xintercept = 6.007020, col="green", linetype="dashed")+
  theme_bw()
blocks_601

#Frequency of kmer mapping
n_occur <- data.frame(table(kmer_mappings$SEQ))
n_occur[n_occur$Freq > 1,]
count(n_occur[n_occur$Freq > 1,])
