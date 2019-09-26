##### ELT-2 ChIP-seq associated genes | Plotting groups by DEGs | Enrichment analysis #####

# Import list of of ELT-2 bound gene (promoters) and DEG lists
ELT2.DEG.Venn <- read.delim(file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Rstudio/20190128 RNAseq vs ChIPseq/genes.boundby.ELT2.L3_DEG.Venn.txt",
                            header = T, stringsAsFactors = F, sep = "\t", strip.white = T)
head(ELT2.DEG.Venn)
tail(ELT2.DEG.Venn)
length(ELT2.DEG.Venn$ELT.2.Up.6)

ELT2.DEG.List <- as.list(ELT2.DEG.Venn)
# tail(ELT2.DEG.List[[1]])
# names(ELT2.DEG.List)
ELT2.DEG.List <- lapply(ELT2.DEG.List, function(x) x[startsWith(x, prefix = "WB")])
str(ELT2.DEG.List)

##########################################################################################################################

################################### ELT-2 G E N E  E X P R E S S I O N  P L O T S ########################################

##########################################################################################################################
library(dplyr)
library(tidyr)
library(doBy)
library(ggplot2)
TissuePlot(x = ELT2.DEG.List$ELT.2.Up.3)
### Loop to create ggplot data for each ELT-2 gene list
# Container for TissuePlot data
ELT2.DEG.plots <- vector("list",length = length(ELT2.DEG.List))
# Looping TissuePlot through list of tissue expressed genes
for (i in 1:length(ELT2.DEG.plots)) {
  tp <- TissuePlot(ELT2.DEG.List[[i]])
  tp2 <- tp + ggtitle(paste0(as.character(names(ELT2.DEG.List[i]))," (",length(ELT2.DEG.List[[i]])," genes)"))
  ELT2.DEG.plots[[i]] <- tp2
}
# Test plot first tissue
ELT2.DEG.plots[[6]]
### Exporting image files for all tissues ###
for (i in 1:length(ELT2.DEG.plots)) {
  ggsave(paste0("ELT2.DEG.plot.",i,".tiff"), plot = ELT2.DEG.plots[[i]], width = 4, height = 4)
}


### Plot Fold Change for ELT-2 bound DEGs ###
# Container for TissuePlot data
ELT2.DEG.plots.FC <- vector("list",length = length(ELT2.DEG.List))
# Looping FoldChangeTissuePlot through list of tissue expressed genes
for (i in 1:length(ELT2.DEG.plots.FC)) {
  
  ELT2.DEG.plots.FC[[i]] <- FoldChangeTissuePlot(ELT2.DEG.List[[i]], title = names(ELT2.DEG.List[i]))
  
}
# Test plot first tissue
ELT2.DEG.plots.FC[[6]]
### Exporting image files for all tissues ###
for (i in 1:length(ELT2.DEG.plots.FC)) {
  ggsave(paste0("ELT2.DEG.FCplot.",i,".FC.tiff"), plot = ELT2.DEG.plots.FC[[i]], width = 4, height = 4)
}

##########################################################################################################################

################################### ELT-2 E N R I C H M E N T  A N A L Y S I S ###########################################

##########################################################################################################################

library("clusterProfiler")
library("org.Ce.eg.db")
str(ELT2.DEG.List)

# MF enrichment for ELT-2 TFBS DEGs
ELT2.enrichMF <- compareCluster(geneClusters = ELT2.DEG.List,
                                  fun = "enrichGO",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE",
                                  ont = "MF")
ELT2.enrichMF # 291 obs.
View(as.data.frame(ELT2.enrichMF))
dotplot(ELT2.enrichMF)
ELT2.enrichMF_2 <- simplify(ELT2.enrichMF)
ELT2.enrichMF_2 # 173 obs.
dotplot(ELT2.enrichMF_2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, size = 10))

ggsave("ELT-2.dotplot.enrichMF_2.tiff", device = "tiff", width = 12, height = 7)

# BP enrichment
ELT2.enrichBP <- compareCluster(geneClusters = ELT2.DEG.List,
                                  fun = "enrichGO",
                                  ont = "BP",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
ELT2.enrichBP # 1540 obs. !?!
head(as.data.frame(ELT2.enrichBP))
dotplot(ELT2.enrichBP)
ELT2.enrichBP_2 <- simplify(ELT2.enrichBP)
ELT2.enrichBP_2 # 478 obs.
dotplot(ELT2.enrichBP_2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, size = 10))

ggsave("ELT-2.dotplot.enrichBP_2.tiff", device = "tiff", width = 9, height = 7)

# CC enrichment
ELT2.enrichCC <- compareCluster(geneClusters = ELT2.DEG.List,
                                  fun = "enrichGO",
                                  ont = "CC",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
ELT2.enrichCC # 579 obs.
head(as.data.frame(ELT2.enrichCC))
dotplot(ELT2.enrichCC)
ELT2.enrichCC_2 <- simplify(ELT2.enrichCC)
ELT2.enrichCC_2 # 157 obs.
dotplot(ELT2.enrichCC_2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, size = 10))

ggsave("ELT-2.dotplot.enrichCC_2.tiff", device = "tiff", width = 7, height = 7)

