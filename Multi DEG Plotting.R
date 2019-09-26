# DESeq2 post analysis plots of expression values with a multipage PDF catalog
# https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf

# The point
# I get asked all the time, “hey, could you make a boxplot of the expression values of SomeGene for me?” The plotCounts() function in DESeq will do this, but only one gene at a time. This shows you how to do this for any arbitrary number of genes in your data (or all of them!). You can plot a small number in a single plot with facets, or you could plot a large number in a text-searchable multi-page PDF.
# 
# Setup
# Load the packages you’ll use.
# 
# # CRAN

# install.packages("rlang")

library(dplyr)
# install.packages("tidyr")
library(tidyr)
library(ggplot2)
# theme_set(theme_bw(base_size=14) + theme(strip.background = element_blank()))
# 
# # Bioconductor
library(DESeq2)

# Create a DESeqDataSet and run the DESeq pipeline on it.
# Inspect the colData.
colData(ddsTC)

# Get the results for the dex treated vs untreated samples and arrange by p-value.

# res <- results(dds, tidy=TRUE, contrast=c("dex", "trt", "untrt")) %>%
#   arrange(padj, pvalue) %>%
#   tbl_df()
# res
resTC

# Now, define some genes of interest. These must be the names of the genes as in the count data. 
# That is, your genes of interest must be  %in% names(dds). These are the genes you’ll plot. You could do as many as you like.
# Since you’ve arranged by p-value, let’s take the top few results (9 in this case, so we make a square faceted plot).

# Define the genes of interest.
goi <- rownames(resTC[1:30,])
stopifnot(all(goi %in% names(ddsTC)))
goi


# Join & tidy
# Create a tidy/transposed count matrix. Here’s a line-by-line explanation:
  
#   Here, you’ll take the counts() of the dds object, normalizing it, without outlier replacement.
#   You’ll add a half count to it, because the next thing you’ll do is log2() it, and you don’t want any -Infs. 
#   So now you have a log-transformed normalized count matrix. Now, transpose that matrix so the sample names are the row.names.
# You’ll now merge that thing with the colData(dds), where the sample names are also the row.names. 
# You’ll now have a really wide data.frame with a single row for each sample, followed by all the colData,
# followed by a column for each of the genes of interest.
# Gather that up again so you now have one row per sample per gene, with all the sample info, and that gene’s expression.




tcounts <- t(log2((counts(ddsTC[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

# Take a look.

tcounts %>% 
  select(Row.names, Strain, HrsPHS, gene, expression) %>% 
  # head %>% 
  knitr::kable()



# Single faceted plot
# Now, create a single faceted plot. You don’t have to add a fill aesthetic - I only did so as an example, 
# using the fake variable I made up.

ggplot(tcounts, aes(HrsPHS, expression, fill=Strain)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="Top Results")

###$### My version: smooth lines instead of boxplots
ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = Strain, color = Strain)) + 
  geom_smooth(se = F) +
  facet_wrap(~gene, scales="free_y") + 
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="Top Results") +
  scale_x_continuous(labels = c("0","3","6","12","30"))


# ggplot(WBGene00004930, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
#   geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00004930 sod-1")

###$$### My version: multi genes - no faceting
ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  geom_smooth(se = F) +
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="Top Results") +
  scale_x_continuous(labels = c("0","3","6","12","30"))

################################################################
#### creating multi DEG plot of upreg(3,6,12,30) - 478 genes ###
################################################################
goi <- read.csv("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/478 common elements in Upreg 3 hrs, 6 hrs, 12 hrs and 30 hrs.csv", header = F, colClasses = "character")
# goi <- as.list(goi)
goi <- c(goi[1:478,1])
goi

tcounts <- t(log2((counts(ddsTC[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

# tcounts <- t(((counts(ddsTC[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
#   merge(colData(ddsTC), ., by="row.names") %>%
#   gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

head(tcounts)

ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  # geom_smooth(se = F, size = 1.5, alpha = 0.1) +
  geom_line(stat = "smooth", method = "loess", alpha = 0.3, size = 1.5) +
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="Upregulated (3,6,12,30) - 478 genes") +
  scale_x_continuous(labels = c("0","3","6","12","30"))



which.max(tcounts$expression)
tcounts[4255,]
(sort(tcounts$expression,decreasing = T))

##########################################################################
# 420 common elements in "Downreg 3 hrs", "6 hrs", "12 hrs" and "30 hrs":
##########################################################################
goi <- read.csv("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/420 common elements in Downreg 3 hrs, 6 hrs, 12 hrs and 30 hrs.csv", header = F, colClasses = "character")
goi[1,] <- "WBGene00020200"
goi <- c(goi[1:420,1])
goi

tcounts <- t(log2((counts(ddsTC[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

head(tcounts)

ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = 0.3, size = 1.5) +
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="Downregulated (3,6,12,30) - 420 genes") +
  scale_x_continuous(labels = c("0","3","6","12","30"))


##########################################################################
### 225 Most Variable Genes Contributing to Principle Component 1 ########
##########################################################################
pca.contrib
pca.contrib.data[pca.contrib.data$contrib>0.2,]
count(pca.contrib.data[pca.contrib.data$contrib>0.5,])

goi <- rownames(pca.contrib.data[pca.contrib.data$contrib>0.5,])
goi

tcounts <- t(log2((counts(ddsTC[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

head(tcounts)

ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = 0.7, size = 1) +
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="19 Most Variable Genes Contributing to PC1") +
  scale_x_continuous(labels = c("0","3","6","12","30"))

# write.csv(tcounts, file = "PC1 top 225 genes count table.csv")



############################################################################################################
######## Multi page PDF of DEG plots #######################################################################
############################################################################################################

# Or in a loop, create a multi-page PDF (text searchable) with one page per gene.

# pdf("multi-ggplot2-catalog.pdf")
# for (i in goi) {
#   p <- ggplot(filter(tcounts, gene==i), aes(dex, expression, fill=fake)) + geom_boxplot() + ggtitle(i)
#   print(p)
# }
# dev.off()
goi[1:10]
resTC["WBGene00000002",c(2,6)]
names(ddsTC)

pdf("multi-ggplot2-catalog.pdf")
for (i in names(ddsTC)) {
  g <- plotCounts(ddsTC, i, intgroup = c("Strain","HrsPHS"), returnData = T)
  p <- ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "loess") + ggtitle(i) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
    xlab("Hours PHS")
  print(p)
}
dev.off()



################################################################################################################################
################################################################################################################################
################################# Finalizing Counts Plots for publication figures ##############################################

library(doBy)
library(scales)

# Get the average FC for a subset of genes
#FC.avg <- summaryBy(FC ~ HrsPHS, data = FC.long[1:4000,], FUN = c(length,mean,sd))
#FC.avg


# elt-2 WBGene00001250
g <- plotCounts(ddsTC, "WBGene00001250", intgroup = c("Strain","HrsPHS"), returnData = T)

s <- summaryBy(count ~ Strain:HrsPHS, data = g)


ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = s, mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("elt-2 WBGene00001250") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))


# elt-7 WBGene00015981
g <- plotCounts(ddsTC, "WBGene00015981", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("elt-7 WBGene00015981") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# pho-1 WBGene00004020 !!!nope. not good example!!!!!
g <- plotCounts(ddsTC, "WBGene00004020", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("pho-1 WBGene00004020") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# ifb-2	WBGene00002054
g <- plotCounts(ddsTC, "WBGene00002054", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("ifb-2	WBGene00002054") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# end-1 WBGene00001310
g <- plotCounts(ddsTC, "WBGene00001310", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("end-1 WBGene00001310") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# end-3	WBGene00001311
g <- plotCounts(ddsTC, "WBGene00001311", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("end-3 WBGene00001311") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# skn-1 WBGene00004804
g <- plotCounts(ddsTC, "WBGene00004804", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("skn-1 WBGene00004804") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# myo-2 WBGene00003514
g <- plotCounts(ddsTC, "WBGene00003514", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("myo-2 WBGene00003514") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# ceh-22 WBGene00000445
g <- plotCounts(ddsTC, "WBGene00000445", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("WBGene00000445") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# ifa-1 WBGene00002050
g <- plotCounts(ddsTC, "WBGene00002050", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("ifa-1 WBGene00002050") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# ule-3 WBGene00008378
g <- plotCounts(ddsTC, "WBGene00008378", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("ule-3 WBGene00008378") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# pes-23 WBGene00003987
g <- plotCounts(ddsTC, "WBGene00003987", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("pes-23 WBGene00003987") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# rsef-1 WBGene00016344
g <- plotCounts(ddsTC, "WBGene00016344", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("rsef-1 WBGene00016344") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

# fkh-6 WBGene00001438
g <- plotCounts(ddsTC, "WBGene00001438", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(g, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  theme_classic() +
  geom_point(size = 2) + 
  geom_line(data = summaryBy(count ~ Strain:HrsPHS, data = g), mapping = aes(y = count.mean), size = 1.5) +
  scale_x_continuous(labels = c("0","3","6","12","30")) +
  scale_y_continuous(labels = scientific) +
  ggtitle("fkh-6 WBGene00001438") +
  xlab("Hours PHS") +
  ylab("TPM") +
  theme(axis.text = element_text(size = 12))

