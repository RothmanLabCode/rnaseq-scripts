##### Venn Diagram Gene Set Overlaps (of UP and Downregulated gene sets) as Clusters for Analysis #####
library("clusterProfiler")
library("org.Ce.eg.db")
# Import lists of DEGs from Venny
library("readxl")
UpDown <- read_excel("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/Up.Down.Venn.Clusters.xlsx")
View(UpDown)

UpDown.list <- as.list(UpDown)
UpDown.list <- lapply(UpDown.list, function(x) x[!is.na(x)])
str(UpDown.list)

# MF enrichment for Upreg & Downreg Venn clusters
UpDown.enrichMF <- compareCluster(geneClusters = UpDown.list,
                                  fun = "enrichGO",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE",
                                  ont = "MF")
UpDown.enrichMF # 108 enriched terms
as.data.frame(UpDown.enrichMF)
dotplot(UpDown.enrichMF)
UpDown.enrichMF_2 <- simplify(UpDown.enrichMF)
UpDown.enrichMF_2 # 47 enriched terms
dotplot(UpDown.enrichMF_2)

# BP enrichment
UpDown.enrichBP <- compareCluster(geneClusters = UpDown.list,
                                  fun = "enrichGO",
                                  ont = "BP",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
UpDown.enrichBP # 556 enriched terms !?!
head(as.data.frame(UpDown.enrichBP))
dotplot(UpDown.enrichBP)
UpDown.enrichBP_2 <- simplify(UpDown.enrichBP)
UpDown.enrichBP_2 # 233 enriched terms
dotplot(UpDown.enrichBP_2)

# CC enrichment
UpDown.enrichCC <- compareCluster(geneClusters = UpDown.list,
                                  fun = "enrichGO",
                                  ont = "CC",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
UpDown.enrichCC # 205 enriched terms
head(as.data.frame(UpDown.enrichCC))
dotplot(UpDown.enrichCC)
UpDown.enrichCC_2 <- simplify(UpDown.enrichCC)
UpDown.enrichCC_2 # 88 enriched terms
dotplot(UpDown.enrichCC_2)

#### Reactome Pathway Enrichment
## Requires converting gene IDs to entrez....
library("ReactomePA")
keytypes(org.Ce.eg.db)
UpDown.bitr <- lapply(UpDown.list, bitr, fromType = "WORMBASE",
                      toType = c("ENTREZID","ALIAS","UNIPROT"),
                      OrgDb = org.Ce.eg.db)
UpDown.bitr$`Up.3.6-Down.12`
# Create empty lists for EntrezIDs
UpDown.eg <- vector("list",length(UpDown.list))
# Add EntrezIDs to Upreg lists
for (i in 1:length(UpDown.eg)) {
  UpDown.eg[[i]] <- UpDown.bitr[[i]]$ENTREZID
}
# Names for EntrezID lists
names(UpDown.eg) <- names(UpDown.list)
### VennCluster Lists with EntrezIDs ###
str(UpDown.eg)
# Reactome Pathway Analysis - RPA
UpDown.RPA <- compareCluster(geneClusters = UpDown.eg,
                                       fun = "enrichPathway",
                                       organism = "celegans")
UpDown.RPA # 69 enriched pathways
head(as.data.frame(UpDown.RPA))
dotplot(UpDown.RPA)

### Exporting csv results files for UpDown clusters enrichment ###
write.csv(as.data.frame(UpDown.enrichBP_2), file = "UpDown.enrichBP_2 20180626.csv")
write.csv(as.data.frame(UpDown.enrichCC_2), file = "UpDown.enrichCC_2 20180626.csv")
write.csv(as.data.frame(UpDown.enrichMF_2), file = "UpDown.enrichMF_2 20180626.csv")
write.csv(as.data.frame(UpDown.RPA), file = "UpDown.RPA 20180626.csv")
