#20180105
#OBJECTIVE: Analysis of mRNAseq Time Course data
#Gene level expression - summary QC plotting, data statistics

install.packages("ggplot2")

setwd("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC")
GeneTPM <- read.csv("Galaxy Output/Genes TPM Table_Samples 1-24.csv",header = T)
head(GeneTPM)

install.packages("beanplot")
library("beanplot")
beanplot(... = GeneTPM)

x <- rnorm(22)
x
typeof(x)
beanplot(x)

GeneTPM[2:24,2:20]
beanplot(GeneTPM[2:24]) #plot works, but takes a few minutes to render, very raw, too broad y-scale

tags <- colnames(GeneTPM)
tags <- (1:24)
tags
beanplot(GeneTPM[2:25], names = tags, cut = 1000)

boxplot(GeneTPM[2:25], names = tags, outline = F, xlab = "Sample", ylab = "TPM", col = "Grey")
