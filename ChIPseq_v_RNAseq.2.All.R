######## Repeating analysis in "ChIPseq_v_RNAseq.R" using the entire list of 217 transcription factors ##########

library(dplyr)
library(tidyr)
library(ggplot2)

head(genome.Summit2Gene, 100)
str(genome.Summit2Gene)
summary(genome.Summit2Gene)
summary(genome.Summit2Gene$Stage)
length(genome.Summit2Gene$Chr) # 444,136 rows

write.csv(genome.Summit2Gene, "genome.Summit2Gene.csv")

### What is the "background" frequency of binding sites for each TF in the ChIP-seq data subset? ###
# Plotted ChIP-seq summits within gene promoter regions
ggplot(data = genome.Summit2Gene, aes(x = TF)) + 
  geom_bar() + 
  theme_bw() +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,30000))

# Counted number of ChIP-seq summits within gene promoter regions
length(levels(genome.Summit2Gene$TF)) # 217
TF.peaks.counts <- count(genome.Summit2Gene, vars = genome.Summit2Gene$TF)
TF.peaks.counts[1:20,]

write.csv(TF.peaks.counts, "TF.peaks.counts_217TF.csv")

# Vector with the names for all 217 Transcription Factors
TF.vector.All <- as.vector(TF.peaks.counts$vars)
TF.vector.All

# Break up the genome wide list of all TF:Gene mappings into separate data frames for each TF, all contained within one List object #

TF.list.All <- as.list(TF.vector.All)
TF.list.All

for (tf in 1:length(TF.list.All)) {
  TF.list.All[[tf]] <- genome.Summit2Gene[genome.Summit2Gene$TF==TF.vector.All[tf],]
}

head(TF.list.All)
## rename lists by TF name
names(TF.list.All) <- TF.vector.All
str(TF.list.All)
TF.list.All$`EOR-1`

### Need the background number of genes which are bound by each TF ###

# TF.unique.genes

length(TF.list.All[[10]]$Gene) # 3123 = number of TF binding sites
length(unique(TF.list.All[[10]]$Gene)) # 3040 = number of GENES with TF binding sites

# starting data.frame container for summarized TF::gene binding data
TF.unique.genes <- data.frame(TF.vector.All)
# looping through the list of TF binding sites and counting the number of unique genes for each
for (i in 1:length(TF.unique.genes$TF.vector.All)) {
  TF.unique.genes$num.genes[i] <- length(unique(TF.list.All[[i]]$Gene))
}

head(TF.unique.genes)

## Create a data frame for a given DEG set (Up.6), where each TF has a row, and columns for:
## (1) number of TF binding sites "num.bs"
## (2) number of genes w/ >= 1 ChIP peak "num.genes"
## (3) percent of the background sites present in the DEG set "pc.bs"
## (4) percent of the DEG set genes which have >= 1 ChIP peak "pc.genes"

Up.3.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Up.3.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Up.3,]
  Up.3.df$num.bs[i] <- length(overlap$Gene)
  Up.3.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.3.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.3.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.3)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Up.6.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Up.6.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Up.6,]
  Up.6.df$num.bs[i] <- length(overlap$Gene)
  Up.6.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.6.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.6.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.6)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Up.12.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Up.12.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Up.12,]
  Up.12.df$num.bs[i] <- length(overlap$Gene)
  Up.12.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.12.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.12.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.12)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Up.30.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Up.30.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Up.30,]
  Up.30.df$num.bs[i] <- length(overlap$Gene)
  Up.30.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.30.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.30.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.30)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.3.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Down.3.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Down.3,]
  Down.3.df$num.bs[i] <- length(overlap$Gene)
  Down.3.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.3.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.3.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.3)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.6.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Down.6.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Down.6,]
  Down.6.df$num.bs[i] <- length(overlap$Gene)
  Down.6.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.6.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.6.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.6)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.12.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Down.12.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Down.12,]
  Down.12.df$num.bs[i] <- length(overlap$Gene)
  Down.12.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.12.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.12.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.12)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.30.df <- data.frame(TF = TF.vector.All)

for (i in 1:length(Down.30.df$TF)) {
  tmp <- TF.list.All[[i]]
  overlap <- tmp[tmp$Gene %in% Down.30,]
  Down.30.df$num.bs[i] <- length(overlap$Gene)
  Down.30.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.30.df$pc.bs[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.30.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.30)
}


##############################################################################################################################

############################## D O T  P L O T S ##############################################################################

##############################################################################################################################

###### Bar graphs & Dotplots for additional DEG sets ######

#### Up.3

### PLOT: % of DEGs with TF binding peaks
ggplot(Up.3.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Upregulated 3 hrsPHS (", length(Up.3), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.3))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Up.3.df[Up.3.df$pc.genes>0.15,], mapping = aes(x = TF, y = pc.genes*100,
                                                      size = Up.3.df$num.bs[Up.3.df$pc.genes>0.15],
                                                      colour = Up.3.df$pc.bs[Up.3.df$pc.genes>0.15]*100)) + 
  geom_point() +
  ggtitle(paste0("Upregulated 3 hrsPHS (", length(Up.3), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(15,30)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Up.3.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#### Up.12

### PLOT: % of DEGs with TF binding peaks
ggplot(Up.12.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Upregulated 12 hrsPHS (", length(Up.12), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.4))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Up.12.df[Up.12.df$pc.genes>0.15,], mapping = aes(x = TF, y = pc.genes*100,
                                                        size = Up.12.df$num.bs[Up.12.df$pc.genes>0.15],
                                                        colour = Up.12.df$pc.bs[Up.12.df$pc.genes>0.15]*100)) + 
  geom_point() +
  ggtitle(paste0("Upregulated 12 hrsPHS (", length(Up.12), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(15,35)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Up.12.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")


#########################################################################################################################################

################################### S T A T S ###########################################################################################

#########################################################################################################################################

##### Statistical analysis of correlation significance using Fisher's Exact tests #####

### for loop to compute p-value for all TFs in a DEG set ###

length(TF.vector.All)

# # # # # # # # # # # # # # # #
# # # Up.3 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Up.3.df$num.genes[TF]
  v <- length(Up.3)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.3.df$pval[TF] <- ct.ft$p.value
}

Up.3.df

# # # # # # # # # # # # # # # #
# # # Up.6 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Up.6.df$num.genes[TF]
  v <- length(Up.6)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.6.df$pval[TF] <- ct.ft$p.value
}

Up.6.df

# # # # # # # # # # # # # # # #
# # # Up.12 fisher tests  # # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Up.12.df$num.genes[TF]
  v <- length(Up.12)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.12.df$pval[TF] <- ct.ft$p.value
}

Up.12.df

# # # # # # # # # # # # # # # #
# # # Up.30 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Up.30.df$num.genes[TF]
  v <- length(Up.30)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.30.df$pval[TF] <- ct.ft$p.value
}

Up.30.df

# # # # # # # # # # # # # # # #
# # # Down.3 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Down.3.df$num.genes[TF]
  v <- length(Down.3)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.3.df$pval[TF] <- ct.ft$p.value
}

Down.3.df

# # # # # # # # # # # # # # # #
# # # Down.6 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Down.6.df$num.genes[TF]
  v <- length(Down.6)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.6.df$pval[TF] <- ct.ft$p.value
}

Down.6.df

# # # # # # # # # # # # # # # #
# # # Down.12 fisher tests  # # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Down.12.df$num.genes[TF]
  v <- length(Down.12)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.12.df$pval[TF] <- ct.ft$p.value
}

Down.12.df

# # # # # # # # # # # # # # # #
# # # Down.30 fisher tests  # #
# # # # # # # # # # # # # # # #

for (TF in 1:217) {
  a <- Down.30.df$num.genes[TF]
  v <- length(Down.30)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.unique.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.30.df$pval[TF] <- ct.ft$p.value
}

Down.30.df

#############################################################################################

####################### P - V A L U E  H E A T  M A P #######################################

#############################################################################################


## Need to create a matrix of data for heatmap. Can also use `as.matrix()` on a data frame.
library(gplots)
# help(heatmap.2)

name.heatmap.cols <- c("Up 3","Up 6","Up 12","Up 30","Down 3","Down 6","Down 12","Down 30")


pvals.matrix <- data.frame(Up.3.df$pval, Up.6.df$pval, Up.12.df$pval, Up.30.df$pval,
                           Down.3.df$pval, Down.6.df$pval, Down.12.df$pval, Down.30.df$pval)
rownames(pvals.matrix) <- TF.vector.All
colnames(pvals.matrix) <- name.heatmap.cols
pvals.matrix <- as.matrix(pvals.matrix)
is.matrix(pvals.matrix)
pvals.matrix[1:10,1:6]
heatmap.2(pvals.matrix)

write.csv(pvals.matrix, "TF.All217.DEGs.pvalues.csv")

pvals.matrix.log <- -log(pvals.matrix)

pvals.matrix.log[is.infinite(pvals.matrix.log)] = 666
is.infinite(pvals.matrix.log)
heatmap.2(pvals.matrix.log)


### Full heatmap of p-values for all 217 TFs
heatmap.2(pvals.matrix.log, trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 4)

# saving image via code
tiff(filename = "217TF.pval.heatmap.05binNo30hrColscale.tiff", width = 6, height = 24, units = "in", res = 300)
dev.off()

### Full heatmap of p-values for all 217 TFs | Version 2
heatmap.2(pvals.matrix.log, trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "col", key = T, keysize = 1.1, density.info = "none", colsep = 4)


# binning large p-values
p <- pvals.matrix.log > 100
pvals.matrix.log.bin <- pvals.matrix.log
pvals.matrix.log.bin[p] = 100

### Full heatmap of p-values for all 217 TFs | Version 3
heatmap.2(pvals.matrix.log.bin, trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 4)

### Full heatmap of p-values for all 217 TFs | Version 4 (no 30 hrsPHS)
heatmap.2(pvals.matrix.log.bin[,c(1,2,3,5,6,7)], trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 3)

### Full heatmap of p-values for all 217 TFs | Version 5 (no 30 hrsPHS)
heatmap.2(pvals.matrix.log.bin[,c(1,2,3,5,6,7)], trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "col", key = T, keysize = 1.1, density.info = "none", colsep = 3)



## trimming TFs with larger pvalues (down to 60)
rowSums(pvals.matrix.log[,c(1,2,3,5,6,7)])
plot(rowSums(pvals.matrix.log[,c(1,2,3,5,6,7)]))
z <- rowSums(pvals.matrix.log[,c(1,2,3,5,6,7)]) > 200
pvals.matrix.log.trim <- pvals.matrix.log[z,]

## trimmed heatmap
heatmap.2(pvals.matrix.log.trim, trace = "none", dendrogram = "row", Colv = F, col = bluered(20),
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 4)

## trimmed heatmap | Version 2 (no 30 hrsPHS)
heatmap.2(pvals.matrix.log.trim[,c(1,2,3,5,6,7)], trace = "none", dendrogram = "row", Colv = F, col = bluered(20),
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 3)

##### Trimming V2 #####

## trimming TFs with larger pvalues (down to 60)
s <- rowSums(pvals.matrix.log[,c(1,2,3)])
plot(s)

z <- rowSums(pvals.matrix.log[,c(1,2,3,5,6,7)]) > 200
pvals.matrix.log.trim <- pvals.matrix.log[s>40,]

## trimmed heatmap
heatmap.2(pvals.matrix.log.trim, trace = "none", dendrogram = "row", Colv = F, col = bluered(20),
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 4)

## trimmed heatmap | Version 2 (no 30 hrsPHS)
heatmap.2(pvals.matrix.log.trim[,c(1,2,3,5,6,7)], trace = "none", dendrogram = "row", Colv = F, col = bluered(20),
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 3)

## trimmed heatmap | Version 3 (binned pval(100); no 30 hrsPHS)
heatmap.2(pvals.matrix.log.bin[s>40,c(1,2,3,5,6,7)], trace = "none", dendrogram = "row", Colv = F, col = bluered(20),
          scale = "col", key = T, keysize = 1.1, density.info = "none", colsep = 3)

## trimmed heatmap | Version 4 (full pval range (up to 666); no 30 hrsPHS; column scaled; TFs with highest p-vals (Up3/6/12))
heatmap.2(pvals.matrix.log[s>50,c(1,2,3,5,6,7)], trace = "none", dendrogram = "row", Colv = F, col = bluered(20),
          scale = "col", key = T, keysize = 1.1, density.info = "none", colsep = 3)




