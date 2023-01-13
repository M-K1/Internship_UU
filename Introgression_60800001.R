#!/usr/bin/Rscript


#Load libraries
library(cowplot)
library(facetscales)
library(ggplot2)
library(dplyr)
library(ggdark)
library(ggpubr)

#Import data

kmer_mappings_6080001 <- read.csv("D:/matthijs/WUR/Internship/Output/introgression_6080001/data_for_manhattan_kmer.csv")

#function for label for blocks
add_x_break <- function(plot, xval) {
  
  p2 <- ggplot_build(plot)
  breaks <- p2$layout$panel_params[[1]]$x$breaks
  breaks <- breaks[!is.na(breaks)]
  
  plot +
    geom_vline(xintercept = xval) +
    scale_x_continuous(breaks = sort(c(xval, breaks)))
}

#Convert chromosome names to actual chromosome number
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022518.1'] <- 1
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022519.1'] <- 2
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022520.1'] <- 3
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022521.1'] <- 4
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022522.1'] <- 5
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022523.1'] <- 6
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022524.1'] <- 7
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022525.1'] <- 8
kmer_mappings_6080001[kmer_mappings_6080001 == 'CM022526.1'] <- 9

kmer_mappings_6080001$P <- -log10(as.numeric(kmer_mappings_6080001$P))

#only kmers in chromosome 1
kmer_mappings_6080001_chr1 <- subset(kmer_mappings_6080001, CHR==1)

#add column with the positions in Mbp
kmer_mappings_6080001_chr1$BP_Mbp <- as.numeric(kmer_mappings_6080001_chr1$BP)/1000000

#Set intervals for first plot
minmax_kmer_pos_6080001 <- range(kmer_mappings_6080001_chr1$BP)
minmax_kmer_pos_6080001_rounded <- round(minmax_kmer_pos_6080001, -3)
intervals_6080001 <- seq(from = minmax_kmer_pos_6080001_rounded[1]-100000,to = (minmax_kmer_pos_6080001_rounded[2]+100000), by = 100000)
intervals_6080001 <- intervals_6080001/1000000


#Create histogram
hist(kmer_mappings_6080001_chr1$BP_Mbp,
     main = "Histogram of introgression 6080001 (bin size 100000)",
     xlab = "Position on chromosome 1 (Mbp",
     breaks = intervals_6080001)

#histogram with new intervals (max 25 Mbp since that is our region of interest)
kmer_mappings_6080001_chr1_25mbp <- subset(kmer_mappings_6080001_chr1, BP_Mbp < 26)
intervals_6080001 <- seq(from = minmax_kmer_pos_6080001_rounded[1]-100000,to = 25000000+100000, by = 100000)
intervals_6080001 <- intervals_6080001/1000000

hist(kmer_mappings_6080001_chr1_25mbp$BP_Mbp,
     main = "Histogram of introgression 6080001 close up (bin size 100000)",
     xlab = "Position on chromosome 1 (Mbp)",
     breaks = intervals_6080001)

#even closer
kmer_mappings_6080001_chr1_5_8mbp <- subset(kmer_mappings_6080001_chr1, BP_Mbp > 5 & BP_Mbp < 8)
intervals_6080001 <- seq(from = 5000000-100000,to = 8000000+100000, by = 10000)
intervals_6080001 <- intervals_6080001/1000000

hist(kmer_mappings_6080001_chr1_5_8mbp$BP_Mbp,
     main = "Histogram of introgression 6080001 5 to 8 Mbp (bin size 10000)",
     xlab = "Position on chromosome 1 (Mbp)",
     breaks = intervals_6080001,
     xlim = c(intervals_6080001[1],intervals_6080001[length(intervals_6080001)])
     )

#Using ggplot
hist_608 <- ggplot(kmer_mappings_6080001_chr1_5_8mbp, aes(x=BP_Mbp)) +
  geom_histogram(binwidth = 0.01, fill="black") +
  # ggtitle("Histogram of introgression 6080001, 5 to 8 Mbp (bin size 10000)") +
  ylab("count") + xlab("Position on chromosome 1 (Mbp)") +
  geom_vline(xintercept=6.080001, linetype='dotted', col = 'red',size=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 6000) +
  theme_bw()

#plot only between region of interest with kmer position and p-value
kmer_mappings_interest_6080001 <- subset(kmer_mappings_6080001_chr1, BP < 7500000 & BP > 5300000)

ggplot(kmer_mappings_interest_6080001,aes(BP_Mbp,P))+
  geom_jitter(height = 0.5)+
  xlab("Position (Mbp)") + ylab("-log10(p)")+
  ggtitle("Jitter plot of significant kmers between 5.25 and 7.5 Mbp of 6080001 (jitter height 0.5)")


#with block lines
blocks_608 <- ggplot(kmer_mappings_interest_6080001,aes(BP_Mbp,P))+
  geom_jitter(height = 0.5, fill = "black", color = "black")+
  xlab("Position (Mbp)") + ylab("-log10(p)")+
  scale_x_continuous(limits = c(5.25, 7.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = max(kmer_mappings_interest_6080001$P) + 10) +
  geom_vline( xintercept = 5.340318, col="yellow", linetype="dashed")+
  geom_vline( xintercept = 6.069797, col="blue", linetype="dashed")+
  geom_vline( xintercept = 6.346654, col="red", linetype="dashed")+
  geom_vline( xintercept = 7.445193, col="green", linetype="dashed")+
  geom_vline( xintercept = 6.904801, col="orange", linetype="dashed")+
  geom_vline( xintercept = 7.143962, col="brown", linetype="dashed")+
  geom_vline( xintercept = 7.382407, col="deeppink", linetype="dashed")+
  geom_vline( xintercept = 6.539012, col="chartreuse3", linetype="dashed")+
  geom_vline( xintercept = 6.736287, col="cyan2", linetype="dashed")+
  geom_vline( xintercept = 5.861800, col="purple", linetype="dashed")+
  geom_vline( xintercept = 5.409798, col="aquamarine4", linetype="dashed")+
  geom_vline( xintercept = 6.272627, col="violet", linetype="dashed")+
  # ggtitle("Jitter plot of significant k-mers between 5.25 and 7.5 Mbp of 6080001") +
  theme_bw()


#Subsample for quicker results
kmer_mappings_interest_6080001_subsample <- sample_frac(kmer_mappings_interest_6080001, size = 0.1)


ggplot(kmer_mappings_interest_6080001_subsample,aes(BP_Mbp,P))+
  geom_jitter(height = 0.5, fill = "black", color = "black")+
  xlab("Position (Mbp)") + ylab("-log10(p)")+
  scale_x_continuous(limits = c(6.5, 6.8)) +
  scale_y_continuous(expand = c(0, 5)) +
  geom_vline( xintercept = 5.340318, col="yellow", linetype="dashed")+
  geom_vline( xintercept = 6.069797, col="blue", linetype="dashed")+
  geom_vline( xintercept = 6.346654, col="red", linetype="dashed")+
  geom_vline( xintercept = 7.445193, col="green", linetype="dashed")+
  geom_vline( xintercept = 6.904801, col="orange", linetype="dashed")+
  geom_vline( xintercept = 7.143962, col="brown", linetype="dashed")+
  geom_vline( xintercept = 7.382407, col="deeppink", linetype="dashed")+
  geom_vline( xintercept = 6.539012, col="chartreuse3", linetype="dashed")+
  geom_vline( xintercept = 6.736287, col="cyan2", linetype="dashed")+
  geom_vline( xintercept = 5.861800, col="purple", linetype="dashed")+
  geom_vline( xintercept = 5.409798, col="aquamarine4", linetype="dashed")+
  geom_vline( xintercept = 6.272627, col="violet", linetype="dashed")+
  ggtitle("Jitter plot of significant kmers (subsample) between 5.25 and 7.5 Mbp of 6080001")

#Frequency of k-mers mapping
n_occur_2 <- data.frame(table(kmer_mappings_6080001$SEQ, kmer_mappings_6080001$BP))
colnames(n_occur_2) <- c("SEQ", "Freq")
n_occur_2[n_occur_2$Freq > 1,]
count(n_occur_2[n_occur_2$Freq > 2,])

map_freq <- table(n_occur_2$Freq)

n_occur_regofint <-data.frame(table(kmer_mappings_interest_6080001$SEQ))
n_occur_regofint[n_occur_regofint$Freq > 1,]
count(n_occur_regofint[n_occur_regofint$Freq > 1,])

seq_and_pos <- data.frame(kmer_mappings_6080001$SEQ, kmer_mappings_6080001$BP, kmer_mappings_6080001$CHR)
colnames(seq_and_pos) <- c("SEQ", "BP", "CHR")

pos_n_occur_2 <- merge(seq_and_pos, n_occur_2, by = "SEQ")
View(pos_n_occur_2)
