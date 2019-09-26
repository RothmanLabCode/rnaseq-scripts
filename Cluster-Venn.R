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

# GO enrichment for Upreg Venn clusters
VennUp.enrichGO <- compareCluster(geneClusters = VennUp.list,
                                  fun = "enrichGO",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
as.data.frame(VennUp.enrichGO)
dotplot(VennUp.enrichGO)

# GO enrichment for Downreg Venn clusters
VennDown.enrichGO <- compareCluster(geneClusters = VennDown.list,
                                    fun = "enrichGO",
                                    OrgDb = org.Ce.eg.db,
                                    keyType = "WORMBASE")
as.data.frame(VennDown.enrichGO)
dotplot(VennDown.enrichGO)



# groupGO returns the number of identified genes belonging to EVERY GO term at the specified level, for each cluster
# no statistical values or "enrichment" returned
# result is very large...
VennUp.groupGO <- compareCluster(geneClusters = VennUp.list,
                                 fun = "groupGO",
                                 ont = "BP",
                                 level = 2,
                                 OrgDb = org.Ce.eg.db,
                                 keytype = "WORMBASE")
as.data.frame(VennUp.groupGO)[,1:5]
length(as.data.frame(VennUp.groupGO)[,1])
# BP level 3: 8,775 rows
# BP level 2: 495 rows
dotplot(VennUp.groupGO)
# Not particularly informative, okay to get an overview, but GO terms with greatest GeneRatios are fairly consistent across clusters.

### KEGG pathway enrichment needs particular IDs to work... need to convert first...
# VennUp.enrichKEGG <- compareCluster(geneClusters = VennUp.list,
#                                     fun = "enrichKEGG",
#                                     organism = "cel",
#                                     # OrgDb = org.Ce.eg.db,
#                                     keyType = "ENSEMBL")

##### Reactome Pathway Analysis / Enrichment #####
# biocLite("ReactomePA")
library("ReactomePA")
# This package is designed for reactome pathway-based analysis. Reactome is an open-source, 
# open access, manually curated and peer-reviewed pathway database.
# Currently ReactomePA supports several model organisms, including ‘celegans’, ‘fly’, ‘human’, ‘mouse’, ‘rat’, ‘yeast’ and ‘zebrafish’.

# The input gene ID should be Entrez gene ID. 
keytypes(org.Ce.eg.db)
VennUp.bitr <- lapply(VennUp.list, bitr, fromType = "WORMBASE",
                      toType = c("ENTREZID","ALIAS","UNIPROT"),
                      OrgDb = org.Ce.eg.db)
# VennUp.bitr$Up.3
VennDown.bitr <- lapply(VennDown.list, bitr, fromType = "WORMBASE",
                        toType = c("ENTREZID","ALIAS","UNIPROT"),
                        OrgDb = org.Ce.eg.db)
# VennDown.bitr$Down.3

# Create empty lists for EntrezIDs
VennUp.eg <- vector("list",length(VennUp.list))
VennDown.eg <- vector("list",length(VennDown.list))
# Add EntrezIDs to Upreg lists
for (i in 1:length(VennUp.eg)) {
  VennUp.eg[[i]] <- VennUp.bitr[[i]]$ENTREZID
}
# Add EntrezIDs to Downreg lists
for (i in 1:length(VennDown.eg)) {
  VennDown.eg[[i]] <- VennDown.bitr[[i]]$ENTREZID
}
# Names for EntrezID lists
names(VennUp.eg) <- names(VennUp.list)
names(VennDown.eg) <- names(VennDown.list)
#
#
### VennCluster Lists with EntrezIDs ###
str(VennUp.eg)
str(VennDown.eg)
#
#
# Reactome Pathway Enrichment using compareCluster
## Upreg enrichPathway
VennUp.enrichPathway <- compareCluster(geneClusters = VennUp.eg,
                                  fun = "enrichPathway",
                                  organism = "celegans")
head(as.data.frame(VennUp.enrichPathway))
## Downreg enrichPathway
VennDown.enrichPathway <- compareCluster(geneClusters = VennDown.eg,
                                         fun = "enrichPathway",
                                         organism = "celegans")
head(as.data.frame(VennDown.enrichPathway))

dotplot(VennUp.enrichPathway, showCategory=NULL)
dotplot(VennDown.enrichPathway, showCategory=10)

### Specific Venn List Pathway Enrichment
VennUp.6.enrichPathway <- enrichPathway(VennUp.eg[["Up.6"]],
                                        organism = "celegans",
                                        readable = T)
dotplot(VennUp.6.enrichPathway)
as.data.frame(VennUp.6.enrichPathway)
emapplot(VennUp.6.enrichPathway, showCategory = 89)
cnetplot(VennUp.6.enrichPathway, showCategory = 30, circular = F, colorEdge = F)
Up.6.Pathway <- simplify(VennUp.6.enrichPathway)


##############################################################################################################################
##### Gene set enrichment analysis (GSEA) ? #####
# The GSEA algorithm using all genes sorted by numerical values (e.g. expression value, fold change etc.).
# see https://github.com/GuangchuangYu/DOSE/wiki/how-to-prepare-your-own-geneList for explanation of preparing your data.
## Copying results data of intrest
res6hr.gseDF <- as.data.frame(res6hr$log2FoldChange)
## Corresponding WB Gene names
res6hr.gseDF$WBGene <- rownames(res6hr)
# Best to rename data columns
names(res6hr.gseDF) = c("FC","WBGene")
# Translating to Entrez IDs
res6hr.bitr <- bitr(res6hr.gseDF$WBGene,
                       fromType = "WORMBASE",
                       toType = "ENTREZID",
                       OrgDb = org.Ce.eg.db)
# Combining results FC data with Entrez IDs
res6hr.gseDF <- merge(res6hr.bitr, res6hr.gseDF, by.x = "WORMBASE", by.y = "WBGene")
### Creating a geneList in format accepted by GSEA
res6hr.gseList <- res6hr.gseDF$FC
names(res6hr.gseList) = res6hr.gseDF$ENTREZID
res6hr.gseList = sort(res6hr.gseList, decreasing = T)
# Splitting into Upreg/Downreg gene lists
res6hr.gseList.UP <- res6hr.gseList[res6hr.gseList > 0]
res6hr.gseList.DOWN <- res6hr.gseList[res6hr.gseList < 0]
res6hr.gseList.DOWN = abs(res6hr.gseList.DOWN)
res6hr.gseList.DOWN = sort(res6hr.gseList.DOWN, decreasing = T)
# gsePathway analysis #
res6hr.gsePathway.DOWN <- gsePathway(res6hr.gseList.DOWN,
                                organism = "celegans",
                                nPerm = 2000,
                                minGSSize = 5,
                                maxGSSize = 1000,
                                pvalueCutoff = 0.5)
(as.data.frame(res6hr.gsePathway.DOWN))
gseaplot(res6hr.gsePathway.DOWN, geneSetID = res6hr.gsePathway.DOWN@result$ID[1])
# CONCLUSIONS: not very good results... not sure why gsea doesn't find pathways, but 'enrichPathway' does find significant pathways...
######################################################################################################################################

##### Additional enrichGO profiles for Venn Clusters #####
# Upreg cluster BP enrichment
VennUp.enrichBP <- compareCluster(geneClusters = VennUp.list,
                                  fun = "enrichGO",
                                  ont = "BP",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
# Upreg Venn clusters enriched Biological Process results
head(as.data.frame(VennUp.enrichBP))
dotplot(VennUp.enrichBP)

# Downreg cluster BPenrichment
VennDown.enrichBP <- compareCluster(VennDown.list,
                                    "enrichGO",
                                    ont = "BP",
                                    OrgDb = org.Ce.eg.db,
                                    keyType = "WORMBASE")
# Downreg Venn clusters enriched Biological Process results
head(as.data.frame(VennDown.enrichBP))
dotplot(VennDown.enrichBP)

# # # # # # # # # # # # #

# Upreg cluster CC enrichment
VennUp.enrichCC <- compareCluster(geneClusters = VennUp.list,
                                  fun = "enrichGO",
                                  ont = "CC",
                                  OrgDb = org.Ce.eg.db,
                                  keyType = "WORMBASE")
# Upreg Venn clusters enriched Cellular Component results
head(as.data.frame(VennUp.enrichCC))
dotplot(VennUp.enrichCC)

# Downreg cluster CC enrichment
VennDown.enrichCC <- compareCluster(VennDown.list,
                                    "enrichGO",
                                    ont = "CC",
                                    OrgDb = org.Ce.eg.db,
                                    keyType = "WORMBASE")
# Downreg Venn clusters enriched Cellular Component results
head(as.data.frame(VennDown.enrichCC))
dotplot(VennDown.enrichCC)

# # # # # # # # # # # # # # # # # # #
### Exporting Enrichment results ###
# MF
write.csv(as.data.frame(VennUp.enrichGO),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennUp.enrichMF 20180515.csv")

write.csv(as.data.frame(VennDown.enrichGO),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennDown.enrichMF 20180515.csv")
# BP
write.csv(as.data.frame(VennUp.enrichBP),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennUp.enrichBP 20180515.csv")

write.csv(as.data.frame(VennDown.enrichBP),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennDown.enrichBP 20180515.csv")
# CC
write.csv(as.data.frame(VennUp.enrichCC),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennUp.enrichCC 20180515.csv")

write.csv(as.data.frame(VennDown.enrichCC),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennDown.enrichCC 20180515.csv")
# Pathway
write.csv(as.data.frame(VennUp.enrichPathway),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennUp.enrichPathway 20180515.csv")

write.csv(as.data.frame(VennDown.enrichPathway),
          file = "~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Cluster-Venn/VennDown.enrichPathway 20180515.csv")

#######################################################################################################################################

##### Sub-analysis of specific Venn-clusters and gene sets #####

## Up.6.12.30 enrichCC -> histone acetyl transferase complex?
Up.6.12.30.enrichCC <- enrichGO(gene = VennUp.list$Up.6.12.30,
                                OrgDb = org.Ce.eg.db,
                                keyType = "WORMBASE",
                                ont = "CC",
                                readable = TRUE)
as.data.frame(Up.6.12.30.enrichCC)
dotplot(Up.6.12.30.enrichCC)
cnetplot(Up.6.12.30.enrichCC, showCategory = 15)

# Down 6.12.30 BP -> electron transport chain?
Down.6.12.30.BP <- enrichGO(VennDown.list$Down.6.12.30,
                            OrgDb = org.Ce.eg.db,
                            keyType = "WORMBASE",
                            ont = "BP",
                            readable = TRUE)
dotplot(Down.6.12.30.BP)
cnetplot(Down.6.12.30.BP, showCategory = 83,
         colorEdge = F, circular = F) # long render time
emapplot(Down.6.12.30.BP)
as.data.frame(Down.6.12.30.BP)$ID

##### Up.3.6.12.30 BP -> Immune response? #####
Up.3.6.12.30.BP <- enrichGO(gene = VennUp.list$Up.3.6.12.30,
                            OrgDb = org.Ce.eg.db,
                            keyType = "WORMBASE",
                            ont = "BP",
                            readable = T)
head(as.data.frame(Up.3.6.12.30.BP))
length(as.data.frame(Up.3.6.12.30.BP)$ID)
dotplot(Up.3.6.12.30.BP, showCategory = 25)
emapplot(Up.3.6.12.30.BP)
goplot(Up.3.6.12.30.BP)
cnetplot(Up.3.6.12.30.BP, showCategory = 25)
Up.3.6.12.30.BP_2 <- simplify(Up.3.6.12.30.BP)
cnetplot(Up.3.6.12.30.BP_2, showCategory = 16)
heatplot(Up.3.6.12.30.BP_2)

# get gene names of intrest ("defense response")
# convert enrichGO result to data frame
Up.3.6.12.30.BP_df <- as.data.frame(Up.3.6.12.30.BP)
# extract string containing genes of intrest from data frame - for GO term "defense response"
Up.3.6.12.30.defense_response <- Up.3.6.12.30.BP_df[Up.3.6.12.30.BP_df$Description=="defense response","geneID"]
# convert text string to column of gene Symbol names
genesofintrest <- as.data.frame(strsplit(Up.3.6.12.30.defense_response, split = "/"), col.names = "SYMBOL")
# get WBGene names back
genesofintrest <- bitr(genesofintrest$SYMBOL, fromType = "SYMBOL", toType = "WORMBASE", OrgDb = org.Ce.eg.db)
stopifnot(all(genesofintrest$WORMBASE %in% names(ddsTC)))
# log normalized gene expression
tcounts <- t(log2((counts(ddsTC[genesofintrest$WORMBASE, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(genesofintrest$WORMBASE)+1):ncol(.))
# scaled gene expression
tcounts <- (scale(t(counts(ddsTC[genesofintrest$WORMBASE], normalized=TRUE, replaced=FALSE)))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(genesofintrest$WORMBASE)+1):ncol(.))

# Plotting of either log2 normalized, or scaled, timecourse gene expression
ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = 0.7, size = 1) +
  labs(x="Hours PHS", 
       # y="Expression (log normalized counts)", 
       y="Expression (scaled normalized counts)", 
       fill="Strain", 
       title="GO:0000981 Up.12.MF") +
  scale_x_continuous(labels = c("0","3","6","12","30"))

########## FUNCTION TO PLOT SCALED GENE EXPRESSION FOR GENE LIST FROM ENRICH RESULT ##########
EnrichedGeneSetPlot <- function(enrichResults, GOterm, tit){
  df <- as.data.frame(enrichResults)
  df.GO <- df[df$ID==GOterm,"geneID"]
  genesofintrest <- as.data.frame(strsplit(df.GO, split = "/"), col.names = "SYMBOL")
  genesofintrest <- bitr(genesofintrest$SYMBOL, fromType = "SYMBOL", toType = "WORMBASE", OrgDb = org.Ce.eg.db)
  stopifnot(all(genesofintrest$WORMBASE %in% names(ddsTC)))
# scaled gene expression
  tcounts <- (scale(t(counts(ddsTC[genesofintrest$WORMBASE], normalized=TRUE, replaced=FALSE)))) %>%
    merge(colData(ddsTC), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(genesofintrest$WORMBASE)+1):ncol(.))
  # ggplot of scaled timecourse gene expression
 p <- ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
    geom_line(stat = "smooth", method = "loess", alpha = 0.7, size = 1) +
    labs(x="Hours PHS", 
         y="Expression (scaled normalized counts)", 
         fill="Strain", 
         title=paste(GOterm, tit, sep = " ")) +
    scale_x_continuous(labels = c("0","3","6","12","30"))
 return(p)
 return(genesofintrest)
}

########## FUNCTION TO RETURN LIST OF GENES FROM ENRICH RESULT ############
EnrichedGeneSetList <- function(enrichResults, GOterm){
  df <- as.data.frame(enrichResults)
  df.GO <- df[df$ID==GOterm,"geneID"]
  genesofintrest <- as.data.frame(strsplit(df.GO, split = "/"), col.names = "SYMBOL")
  genesofintrest <- bitr(genesofintrest$SYMBOL, fromType = "SYMBOL", toType = "WORMBASE", OrgDb = org.Ce.eg.db)
  stopifnot(all(genesofintrest$WORMBASE %in% names(ddsTC)))
  return(genesofintrest)
}

### test of EnrichedGeneSetPlot function ###
View(as.data.frame(Up.12.MF))
# GO:0005509 calcium ion binding
EnrichedGeneSetPlot(Up.12.MF, "GO:0005509", tit = "Calcium Ion Binding (Up.12)")

##### Upreg 12 hrs MF -> sequence-specific DNA binding? #####
Up.12.MF <- enrichGO(VennUp.list$Up.12,
                            OrgDb = org.Ce.eg.db,
                            keyType = "WORMBASE",
                            ont = "MF",
                            readable = TRUE)
Up.12.MF #23 enriched terms found
emapplot(Up.12.MF)
heatplot(Up.12.MF)
Up.12.MF_2 <- simplify(Up.12.MF, cutoff = 0.5) # reduced to 9 enriched terms
heatplot(Up.12.MF_2)
cnetplot(Up.12.MF_2, showCategory = 9)

# Gene expression plot of "RNA polymerase II TF activity..."
Up.12.MF_df <- as.data.frame(Up.12.MF)
Up.12.MF_df.pol2 <- Up.12.MF_df[Up.12.MF_df$ID=="GO:0000981","geneID"]
genesofintrest <- as.data.frame(strsplit(Up.12.MF_df.pol2, split = "/"), col.names = "SYMBOL")
genesofintrest <- bitr(genesofintrest$SYMBOL, fromType = "SYMBOL", toType = "WORMBASE", OrgDb = org.Ce.eg.db)
stopifnot(all(genesofintrest$WORMBASE %in% names(ddsTC)))
write.csv(genesofintrest, file = "Up.12.MF RNA Pol II genes.csv", quote = F)
genesofintrest

EnrichedGeneSetPlot(enrichResults = Up.12.MF_2, GOterm = "GO:0000981", tit = "RNA pol II TF activity (Up 12 MF)")


##### UPreg 12 hrs BP -> neurogenesis? #####
Up.12.BP <- enrichGO(VennUp.list$Up.12,
                     OrgDb = org.Ce.eg.db,
                     keyType = "WORMBASE",
                     ont = "BP",
                     readable = TRUE)
Up.12.BP # 148 enriched terms!
emapplot(Up.12.BP, showCategory = 148)
dotplot(Up.12.BP)
Up.12.BP_2 <- simplify(Up.12.BP) # reduced to 47 enriched terms
emapplot(Up.12.BP_2)
cnetplot(Up.12.BP_2, showCategory = 47) #, circular = T, colorEdge = T)
heatplot(Up.12.BP_2)
View(as.data.frame(Up.12.BP_2))
EnrichedGeneSetPlot(enrichResults = Up.12.BP_2, GOterm = "GO:0007399", tit = "nervous system dev. (Up 12)")
FoldChangeTissuePlot(genesofintrest$WORMBASE, title = "Nervous System Development (Up 12)")

##### Upreg 6 hrs CC -> Proteasome? #####
Up.6.CC <- enrichGO(VennUp.list$Up.6,
                    OrgDb = org.Ce.eg.db,
                    keyType = "WORMBASE",
                    ont = "CC",
                    readable = T)
Up.6.CC #73 enriched terms
emapplot(Up.6.CC, showCategory = 73)
Up.6.CC_2 <- simplify(Up.6.CC) #17 enriched terms
dotplot(Up.6.CC_2, showCategory = 20)
cnetplot(Up.6.CC_2, showCategory = 17)
heatplot(Up.6.CC_2)
# get gene names of intrest ("proteasome complex")
# convert enrichGO result to data frame
Up.6.CC_df <- as.data.frame(Up.6.CC_2)
# extract string containing genes of intrest from data frame - for GO term "defense response"
Up.6.CC.proteasome <- Up.6.CC_df[Up.6.CC_df$Description=="proteasome complex","geneID"]
# convert text string to column of gene Symbol names
genesofintrest <- as.data.frame(strsplit(Up.6.CC.proteasome[1], split = "/"), col.names = "SYMBOL")
# get WBGene names back
genesofintrest <- bitr(genesofintrest$SYMBOL, fromType = "SYMBOL", toType = "WORMBASE", OrgDb = org.Ce.eg.db)
stopifnot(all(genesofintrest$WORMBASE %in% names(ddsTC)))


##### Downreg 6.12.30 Biological Process -> oxphos genes #####
Down.6.12.30.BP <- enrichGO(VennDown.list$Down.6.12.30,
                            OrgDb = org.Ce.eg.db,
                            keyType = "WORMBASE",
                            ont = "BP",
                            readable = TRUE)
Down.6.12.30.BP #83 enriched terms
Down.6.12.30.BP_2 <- simplify(Down.6.12.30.BP)
Down.6.12.30.BP_2 # 20 enriched terms
cnetplot(Down.6.12.30.BP_2)#, showCategory = 20)
heatplot(Down.6.12.30.BP_2)
View(as.data.frame(Down.6.12.30.BP_2))
EnrichedGeneSetPlot(enrichResults = Down.6.12.30.BP_2,
                    GOterm = "GO:0022900",
                    tit = "electron transport chain")
goi <- EnrichedGeneSetList(enrichResults = Down.6.12.30.BP_2, GOterm = "GO:0022900")
FoldChangeTissuePlot(goi$WORMBASE, title = "electron transport chain (18 genes)")


##### Downreg 6.12.30 Biological Process -> muscle contraction #####
EnrichedGeneSetPlot(enrichResults = Down.6.12.30.BP_2,
                    GOterm = "GO:0006936",
                    tit = "muscle contraction")
goi <- EnrichedGeneSetList(enrichResults = Down.6.12.30.BP_2, GOterm = "GO:0006936")
FoldChangeTissuePlot(goi$WORMBASE, title = "muscle contraction (10 genes)")


##### Down.12.30 BP Gonad development #####
# EnrichedGeneSetPlot(enrichResults = VennDown.enrichBP,
#                     GOterm = "GO:0008406",
#                     tit = "gonad development")

VennDown.enrichBP_df <- as.data.frame(VennDown.enrichBP)
Down.12.30.gonad.dev <- VennDown.enrichBP_df[VennDown.enrichBP_df$ID=="GO:0008406","geneID"]
genesofintrest <- as.data.frame(strsplit(Down.12.30.gonad.dev, split = "/"), col.names = "WBGene", stringsAsFactors = F)
# TissuePlot(genesofintrest)
# scaled gene expression
tcounts <- (scale(t(counts(ddsTC[genesofintrest$WBGene], normalized=TRUE, replaced=FALSE)))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(genesofintrest$WBGene)+1):ncol(.))
# ggplot of scaled timecourse gene expression
ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = 0.7, size = 1) +
  labs(x="Hours PHS", 
       y="Expression (scaled normalized counts)", 
       fill="Strain", 
       title="GO:0008406 gonad development") +
  scale_x_continuous(labels = c("0","3","6","12","30"))
