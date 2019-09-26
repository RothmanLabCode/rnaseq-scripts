####################################################################################################
##### Modification of "Cluster-Venn.R" script for publication clustering GO enrichment figures #####
####################################################################################################

##### Venn Diagram Gene Set Overlaps as Clusters for Analysis #####
library("clusterProfiler")
library("org.Ce.eg.db")

# Import lists of DEGs from Venny
library("readxl")
VennUp <- read_excel("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/Upreg Venn Clusters.xlsx")
View(VennUp)
VennDown <- read_excel("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/Downreg Venn Clusters.xlsx")
View(VennDown)

VennUp.list <- as.list(VennUp)
VennUp.list <- lapply(VennUp.list, function(x) x[!is.na(x)])
str(VennUp.list)
VennUp.list$Up.3.6.30

VennDown.list <- as.list(VennDown)
VennDown.list <- lapply(VennDown.list, function(x) x[!is.na(x)])
str(VennDown.list)

# Upreg cluster BP enrichment
VennUp.enrichBP <- compareCluster(geneClusters = VennUp.list,
                                  fun = "enrichGO",
                                  ont = "BP",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
# Upreg Venn clusters enriched Biological Process results
head(as.data.frame(VennUp.enrichBP))
dotplot(VennUp.enrichBP)
VennUp.enrichBP #705 obs.

# I want to reduce GO term redundancy using 'simplify()' and decrease number of terms in dotplots.
VennUp.enrichBP.S <- simplify(VennUp.enrichBP)
VennUp.enrichBP.S #297 obs.
?dotplot
dotplot(VennUp.enrichBP.S)
dotplot(VennUp.enrichBP.S, showCategory = 3, font.size = 11)

###########################################
### Now to do the same for Downreg DEGs ###
###########################################

# Downreg cluster BPenrichment
VennDown.enrichBP <- compareCluster(VennDown.list,
                                    "enrichGO",
                                    ont = "BP",
                                    OrgDb = org.Ce.eg.db,
                                    keyType = "WORMBASE")
# Downreg Venn clusters enriched Biological Process results
head(as.data.frame(VennDown.enrichBP))
dotplot(VennDown.enrichBP)
VennDown.enrichBP #1369 obs.

# reduce redundancy and plotted terms
VennDown.enrichBP.S <- simplify(VennDown.enrichBP)
VennDown.enrichBP.S #412 obs.
dotplot(VennDown.enrichBP.S)
dotplot(VennDown.enrichBP.S, showCategory = 3, font.size = 11)





#######################################################################################################
#######################################################################################################
##### 2019.09.09 new code: tissue plot for ALL 20K genes ##############################################

all.genes <- TissuePlot(names(ddsTC)) + ggtitle("All Genes (20094)")
all.genes



tcounts <- scale(t(counts(ddsTC, normalized = T, replaced = F) + 0.5)) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(ddsTC)+1):ncol(.))
# 416 rows corresponding to 13 genes have NaN expression values. Need to remove for summaryBy to work
### (why are there NaN genes? On wormbase little info, may not be real genes, only predicted but never transcribed?!)

tcounts[!is.nan(tcounts$expression),]

tcounts.sum <- summaryBy(expression ~ Strain + HrsPHS, data = tcounts[!is.nan(tcounts$expression),], FUN = c(length,mean,sd))

p <- ggplot(tcounts, aes(x = HrsPHS, y = expression, group = interaction(gene, Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = .5) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_grid(Strain ~ .) +
  labs(x="Hours PHS", 
       y="Expression (scaled normalized TPM)", 
       fill="Strain") +
  geom_line(data = tcounts.sum, aes(x = HrsPHS, y = expression.mean, group = Strain), color = "black") +
  ggtitle("All Genes (20094)")

p

ggsave("All Genes.tissueplot.tiff", plot = p, width = 4, height = 4)
library(svglite)
ggsave("All Genes.tissueplot.svg", plot = p, device = "svg", width = 4, height = 4)

##########################################################################################################
##########################################################################################################


