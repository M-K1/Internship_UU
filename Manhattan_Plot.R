#!/usr/bin/Rscript


#Load libraries
library(cowplot)
library(facetscales)
library(ggplot2)

#Import data
args = commandArgs(trailingOnly=TRUE)
gwas.res <- read.csv(args[1])
threshold_file <- read.table(args[3])

#Convert chromosome names to actual chromosome number
gwas.res[gwas.res == 'NC_056623.1'] <- 1
gwas.res[gwas.res == 'NC_056624.1'] <- 2
gwas.res[gwas.res == 'NC_056625.1'] <- 3
gwas.res[gwas.res == 'NC_056626.1'] <- 4
gwas.res[gwas.res == 'NC_056627.1'] <- 5
gwas.res[gwas.res == 'NC_056628.1'] <- 6
gwas.res[gwas.res == 'NC_056629.1'] <- 7
gwas.res[gwas.res == 'NC_056630.1'] <- 8
gwas.res[gwas.res == 'NC_056631.1'] <- 9


gwas.res[gwas.res == 'CM022518.1'] <- 1
gwas.res[gwas.res == 'CM022519.1'] <- 2
gwas.res[gwas.res == 'CM022520.1'] <- 3
gwas.res[gwas.res == 'CM022521.1'] <- 4
gwas.res[gwas.res == 'CM022522.1'] <- 5
gwas.res[gwas.res == 'CM022523.1'] <- 6
gwas.res[gwas.res == 'CM022524.1'] <- 7
gwas.res[gwas.res == 'CM022525.1'] <- 8
gwas.res[gwas.res == 'CM022526.1'] <- 9

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



## Manhattan plot

gwas.res$BP <- as.numeric(gwas.res$BP)/1000000 #Position from basepairs to mega basepairs

gwas.res <- rbind(gwas.res, list('seq_chr_1', '1', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_2', '2', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_3', '3', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_4', '4', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_5', '5', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_6', '6', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_7', '7', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_8', '8', '1', '1'))
gwas.res <- rbind(gwas.res, list('seq_chr_9', '9', '1', '1'))

gwas.res$CHR <- as.numeric(gwas.res$CHR) #To make sure it is not a string
gwas.res$BP <- as.numeric(gwas.res$BP) #Make numeric
gwas.res$P <- -log10(as.numeric(gwas.res$P)) #If needed transform to -log10

threshold <- threshold_file$V1

mhpl <- ggplot(gwas.res,aes(BP,P))+
  geom_point(alpha=0.7)+
  facet_grid_sc(~CHR,scales=list(x = scales_x,y="free"),space="free")+
  xlab("Position (Mbp)") + ylab("-log10(p)") +
  ggtitle(args[4])+
  theme_cowplot()+
  theme(panel.border = element_rect(colour = "black",linetype = "solid"),
        panel.spacing = unit(0.1,"cm"), axis.text.x = element_text(angle = 75, hjust = 1, size = 10))+
  scale_x_continuous(breaks=c(0, 50, 100, 200, 300))+
  geom_hline(yintercept=threshold, linetype='dotted', col = 'red',size=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),  limits = c(2, NA))


#Create pdf of manhattan plot
ggsave(args[2], dpi = 300, device = "png", bg='white')

