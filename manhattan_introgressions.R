#!/usr/bin/Rscript


#Load libraries
library(cowplot)
library(facetscales)
library(ggplot2)
library(ggpubr)

#Import data
args = commandArgs(trailingOnly=TRUE)
intro_601.res <- read.csv("D:/matthijs/WUR/Internship/Output/introgression_6010001/data_for_manhattan_kmer.csv")
intro_608.res <- read.csv("D:/matthijs/WUR/Internship/Output/introgression_6080001/data_for_manhattan_kmer.csv")
threshold_601 <- read.table("D:/matthijs/WUR/Internship/Output/introgression_6010001/threshold_5per")
threshold_608 <- read.table("D:/matthijs/WUR/Internship/Output/introgression_6080001/threshold_5per")

#Convert chromosome names to actual chromosome number
intro_601.res[intro_601.res == 'CM022518.1'] <- 1
intro_601.res[intro_601.res == 'CM022519.1'] <- 2
intro_601.res[intro_601.res == 'CM022520.1'] <- 3
intro_601.res[intro_601.res == 'CM022521.1'] <- 4
intro_601.res[intro_601.res == 'CM022522.1'] <- 5
intro_601.res[intro_601.res == 'CM022523.1'] <- 6
intro_601.res[intro_601.res == 'CM022524.1'] <- 7
intro_601.res[intro_601.res == 'CM022525.1'] <- 8
intro_601.res[intro_601.res == 'CM022526.1'] <- 9

intro_608.res[intro_608.res == 'CM022518.1'] <- 1
intro_608.res[intro_608.res == 'CM022519.1'] <- 2
intro_608.res[intro_608.res == 'CM022520.1'] <- 3
intro_608.res[intro_608.res == 'CM022521.1'] <- 4
intro_608.res[intro_608.res == 'CM022522.1'] <- 5
intro_608.res[intro_608.res == 'CM022523.1'] <- 6
intro_608.res[intro_608.res == 'CM022524.1'] <- 7
intro_608.res[intro_608.res == 'CM022525.1'] <- 8
intro_608.res[intro_608.res == 'CM022526.1'] <- 9


#Hard coded scale
scales_x <- list(
  `1` = scale_x_continuous(limits = c(0, 214.8)),
  `2` = scale_x_continuous(limits = c(0, 217.1)),
  `3` = scale_x_continuous(limits = c(0, 257.8)),
  `4` = scale_x_continuous(limits = c(0, 337.4)),
  `5` = scale_x_continuous(limits = c(0, 339.6)),
  `6` = scale_x_continuous(limits = c(0, 193.1)),
  `7` = scale_x_continuous(limits = c(0, 195.5)),
  `8` = scale_x_continuous(limits = c(0, 309.6)),
  `9` = scale_x_continuous(limits = c(0, 203.9))
)

## Manhattan plot 601

intro_601.res$BP <- as.numeric(intro_601.res$BP)/1000000 #Position from basepairs to mega basepairs

intro_601.res <- rbind(intro_601.res, list('seq_chr_1', '1', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_2', '2', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_3', '3', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_4', '4', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_5', '5', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_6', '6', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_7', '7', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_8', '8', '1', '1'))
intro_601.res <- rbind(intro_601.res, list('seq_chr_9', '9', '1', '1'))

intro_601.res$CHR <- as.numeric(intro_601.res$CHR) #To make sure it is not a string
intro_601.res$BP <- as.numeric(intro_601.res$BP) #Make numeric
intro_601.res$P <- -log10(as.numeric(intro_601.res$P)) #If needed transform to -log10

tr_601 <- threshold_601$V1


mhpl_601 <- ggplot(intro_601.res,aes(BP,P))+
  geom_point(alpha=0.7)+
  facet_grid_sc(~CHR,scales=list(x = scales_x,y="free"),space="free")+
  xlab("Position (Mbp)") + ylab("-log10(p)") +
  # ggtitle(args[1])+
  theme_cowplot()+
  theme(panel.border = element_rect(colour = "black",linetype = "solid"),
        panel.spacing = unit(0.1,"cm"), axis.text.x = element_text(angle = 75, hjust = 1, size = 10))+
  scale_x_continuous(breaks=c(0, 50, 100, 200, 300))+
  geom_hline(yintercept=tr_601, linetype='dotted', col = 'red',size=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),  limits = c(2, NA))


## Manhattan plot 608

intro_608.res$BP <- as.numeric(intro_608.res$BP)/1000000 #Position from basepairs to mega basepairs

intro_608.res <- rbind(intro_608.res, list('seq_chr_1', '1', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_2', '2', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_3', '3', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_4', '4', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_5', '5', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_6', '6', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_7', '7', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_8', '8', '1', '1'))
intro_608.res <- rbind(intro_608.res, list('seq_chr_9', '9', '1', '1'))

intro_608.res$CHR <- as.numeric(intro_608.res$CHR) #To make sure it is not a string
intro_608.res$BP <- as.numeric(intro_608.res$BP) #Make numeric
intro_608.res$P <- -log10(as.numeric(intro_608.res$P)) #If needed transform to -log10

tr_608 <- threshold_601$V1

mhpl_608 <- ggplot(intro_608.res,aes(BP,P))+
  geom_point(alpha=0.7)+
  facet_grid_sc(~CHR,scales=list(x = scales_x,y="free"),space="free")+
  xlab("Position (Mbp)") + ylab("-log10(p)") +
  # ggtitle(args[1])+
  theme_cowplot()+
  theme(panel.border = element_rect(colour = "black",linetype = "solid"),
        panel.spacing = unit(0.1,"cm"), axis.text.x = element_text(angle = 75, hjust = 1, size = 10))+
  scale_x_continuous(breaks=c(0, 50, 100, 200, 300))+
  geom_hline(yintercept=tr_608, linetype='dotted', col = 'red',size=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),  limits = c(2, NA))

#Combine all plots
combine_introgression <- ggarrange(mhpl_601, mhpl_608, 
      ggarrange(hist_601, hist_608, blocks_601, blocks_608, labels = c("C", "D", "E", "F"), nrow = 2, ncol = 2 ),  
      nrow = 3, labels = c("A", "B"))


combine_introgression_test <- ggarrange(mhpl_601, mhpl_601, 
                                   ggarrange(mhpl_601, mhpl_601, mhpl_601, mhpl_601, labels = c("C", "D", "E", "F"), nrow = 2, ncol = 2 ),  
                                   nrow = 3, labels = c("A", "B"))

#Save that shit
ggsave("intro_combined.png", plot = combine_introgression, dpi = 300, device = "png", bg='white', units = "cm", width = 40, height = 40 )
