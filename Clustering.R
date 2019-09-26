# Clustering

# We can furthermore cluster significant genes by their profiles. We extract a matrix of the
# shrunken log2 fold changes using the coef function:
betas <- coef(ddsTC)
colnames(betas)
head(betas)

# We can now plot the log2 fold changes in a heatmap (figure below).
mat <- betas[, -c(1,2)]
mat <- mat[complete.cases(mat),] #need to remove any rows with NA
summary(mat)

# http://www.2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html
# The first measure is using the sum of squared error (SSE). SSE is defined as the sum of the
# squared distance between each member of a cluster and its cluster centroid. We repeatedly test
# and increasing number of clusters and evaluate the SSE. As we increase the number of clusters
# the distance between any point and it’s centroid will be smaller since the cluster itself is smaller.
# At a certain number of clusters number however, the SSE will not significantly decrease with each new
# addition of a cluster. This is the elbow and suggests a suitable number of clusters:
clust.data <- mat[,5:8]
head(clust.data, 15)
wss <- (nrow(clust.data)-1)*sum(apply(clust.data,2,var))
wss
for (i in 2:40) wss[i] <- sum(kmeans(clust.data, centers=i)$withinss)
plot(1:40, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


# Determining Cluster Number with NbClust package
# chooseCRANmirror()
# source("https://cloud.r-project.org/")
# install.packages("NbClust")
library("NbClust")

# NbClust(data = clust.data, min.nc = 3, max.nc = 25, method = "kmeans", index = "all") #Failed

head(clust.data)
abs(clust.data[1,3:4])

head(rowSums2(abs(clust.data)))
plot(rowSums2(abs(clust.data)))

logFCsums <- rowSums2(abs(clust.data))
summary(logFCsums)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1923  2.2981  4.1526  5.8914  7.7382 44.9533 
length(logFCsums[logFCsums>8.655])

# logFCsums <- logFCsums[logFCsums>8.655]
# plot(logFCsums)

### Trimming the data for genes with top 4,000 greatest absolute amount of fold changes
logFCsums <- logFCsums > 8.655
clust.data.FCtrim <- clust.data[logFCsums,]
head(clust.data.FCtrim)

NbClust(clust.data.FCtrim, min.nc = 3, max.nc = 10, method = "kmeans", index = "all")
# ******************************************************************* 
#   * Among all indices:                                                
#   * 4 proposed 3 as the best number of clusters 
# * 11 proposed 4 as the best number of clusters 
# * 3 proposed 5 as the best number of clusters 
# * 2 proposed 7 as the best number of clusters 
# * 1 proposed 8 as the best number of clusters 
# * 1 proposed 9 as the best number of clusters 
# * 1 proposed 10 as the best number of clusters 
# 
# ***** Conclusion *****                            
#   
#   * According to the majority rule, the best number of clusters is  4 
# 
# 
# ******************************************************************* 
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
pheatmap(clust.data.FCtrim, cluster_col=FALSE, show_rownames = F)
pheatmap(clust.data.FCtrim, cluster_col=FALSE, show_rownames = F, kmeans_k = 12)


wss <- (nrow(clust.data.FCtrim)-1)*sum(apply(clust.data.FCtrim,2,var))
wss
for (i in 2:40) wss[i] <- sum(kmeans(clust.data.FCtrim, centers=i, iter.max = 20)$withinss)
plot(1:40, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


library("cluster")
# set.seed(13)
gap <- clusGap(clust.data.FCtrim, kmeans, 100, B = 100, verbose = interactive())
plot(gap, main = "Gap statistic")
abline(v=which.max(gap$Tab[,3]), lty = 2)
# par(mfrow = c(1,1))
gap$Tab
# cluster number: 15, 17, 18, 20, 22, 23


# # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################################### APclustering #########################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # #

# install.packages("apcluster")
library("apcluster")
vignette("apcluster") #user manual
help(apcluster)
help("APResult")

apresults1 <- apcluster(negDistMat(r=2), clust.data.FCtrim)
apresults1
# Number of samples     =  4000 
# Number of iterations  =  243 
# Input preference      =  -84.76268 
# Sum of similarities   =  -12642.07 
# Sum of preferences    =  -7374.353 
# Net similarity        =  -20016.43 
# Number of clusters    =  87 
apresults1@sim
plot(apresults1, clust.data.FCtrim)
heatmap(apresults1)
# Default results look okay, but a fairly high cluster number, and the heatmap seems to show far fewer clusters, maybe 4-6?

# Merging clusters obtained from affinity propogation using Exemplar-based agglomerative clustering
apresults1agg <- aggExCluster(x = apresults1)
apresults1agg
plot(apresults1agg) #Dendrogram of the 87 clusters
# plot(apresults1agg, showSamples=TRUE, nodePar=list(pch=NA, lab.cex=0.4))
plot(apresults1agg, clust.data.FCtrim, k=12)
apresults1aggCut20 <- cutree(apresults1agg, k=20)
apresults1aggCut20
heatmap(apresults1aggCut20)

#Multi ggplot of Exemplars for 20 cluster results
goi <- names(apresults1aggCut20@exemplars)
goi

tcounts <- t(log2((counts(ddsTC[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

head(tcounts[tcounts$Strain=="JR3642",])
tcounts <- tcounts[tcounts$Strain=="JR3642",]

ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = 1, size = 1.3) +
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="20 Exemplars") +
  scale_x_continuous(labels = c("0","3","6","12","30"))


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

#Reclustering with minimum input preference
apcSimMat <- negDistMat(clust.data.FCtrim, r=2)
apcSimMat[1:9,1:9]

apresults2 <- apcluster(s = apcSimMat, q = 0, details = TRUE)
apresults2
# Number of samples     =  4000 
# Number of iterations  =  142 
# Input preference      =  -2028.373 
# Sum of similarities   =  -35049.88 
# Sum of preferences    =  -24340.47 
# Net similarity        =  -59390.36 
# Number of clusters    =  12 
plot(apresults2, clust.data.FCtrim)
#plot legend
plot.new()
par(oma=rep(0.1,4),mar=rep(0.1,4))
legend("center", legend = paste("Cluster",1:length(apresults2)),col=rainbow(length(apresults2)),pch=19,cex=1.5)

plot(apresults2)
heatmap(apresults2, apcSimMat)

# Merging clusters obtained from affinity propogation using Exemplar-based agglomerative clustering
apresults2agg <- aggExCluster(s = apcSimMat, x = apresults2)
apresults2agg
plot(apresults2agg) #Dendrogram of the 12 clusters
heatmap(apresults2agg, apcSimMat, cexRow=0, cexCol=0, legend="col", col=rev(brewer.pal(9,"YlOrRd")))
# display.brewer.all()
rev(brewer.pal(9,"YlOrRd"))

#Multi ggplot of Exemplars for 12 cluster results
goi <- names(apresults2@exemplars)
goi

tcounts <- t(log2((counts(ddsTC[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

head(tcounts[tcounts$Strain=="JR3642",],20)
# tcounts <- tcounts[tcounts$Strain=="JR3642",]

ggplot(tcounts, aes(x = as.numeric(HrsPHS), y = expression, group = interaction(gene,Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = 1, size = 1.3) +
  labs(x="Hours PHS", 
       y="Expression (log normalized counts)", 
       fill="Strain", 
       title="20 Exemplars") +
  scale_x_continuous(labels = c("0","3","6","12","30"))


#################################################################
##### Plotting Cluster Profiles for apresults2 (12 clusters) ####
#
# Below code notes: make sure scaling method is correct. Loop clusters, plot all.
apres2Exemplars <- names(apresults2@exemplars)
apres2Exemplars[4]

exemplarsCounts <- (scale(t(counts(ddsTC[apres2Exemplars], normalized=TRUE, replaced=FALSE)))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(apres2Exemplars)+1):ncol(.))
head(exemplarsCounts)
# Select count data for single exemplar gene
exemplarsCounts[exemplarsCounts$gene==apres2Exemplars[1],]


# Need to have centered, scaled gene expression
counts(ddsTC["WBGene00006555",], normalized=T)

scale(counts(ddsTC[apresults2@clusters[[2]]], normalized=T))

###########################################
##### Cluster Plotting Final Pipeline #####
goi <- names(apresults2@clusters[[12]])

tcounts <- (scale(t(counts(ddsTC[goi], normalized=TRUE, replaced=FALSE)))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

# head(tcounts)
# head(tcounts[tcounts$Strain=="JR3642",],20)
# tcounts <- tcounts[tcounts$Strain=="JR3642",]

ggplot(tcounts, aes(x = HrsPHS, y = expression, group = interaction(gene, Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = .7) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_grid(Strain ~ .) +
  labs(x="Hours PHS", 
       y="Expression (scaled normalized counts)", 
       fill="Strain", 
       title=paste("Cluster 12 (",length(goi)," genes)", sep = "")) +
  geom_smooth(se=F, data=exemplarsCounts[exemplarsCounts$gene==apres2Exemplars[12],], color="black")

#c(rep("grey",5),rep("black",5)), size=1.5)
#  scale_x_continuous(labels = c("0","3","6","12","30"))
# paste("Cluster 1 (",length(goi)," genes)", sep = "")

###################################################################################################
##### Looping and multi pane plotting #############################################################
apres2Exemplars #gene names - from above code
exemplarsCounts #count data for exemplars - from above code

# attempting to loop first 4 clusters
?par
par(mfrow = c(2,2))

for (i in 1:4) {
  goiLoop <- names(apresults2@clusters[[i]])
  
  tcountsLoop <- (scale(t(counts(ddsTC[goiLoop], normalized=TRUE, replaced=FALSE)))) %>%
    merge(colData(ddsTC), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goiLoop)+1):ncol(.))
  
  p <- ggplot(tcountsLoop, aes(x = HrsPHS, y = expression, group = interaction(gene, Strain), color = Strain)) + 
    geom_line(stat = "smooth", method = "loess", alpha = .7) +
    theme_minimal() +
    labs(x="Hours PHS", 
         y="Expression (scaled normalized counts)", 
         fill="Strain", 
         title= paste("Cluster ", i," [",length(goiLoop)," genes]",sep = "")) +
    geom_smooth(se=F, data=exemplarsCounts[exemplarsCounts$gene==apres2Exemplars[i],], 
                color=c(rep("grey",5),rep("black",5)), size=1.5)
  print(p)

  }
### DOESN'T REALLY WORK AT ALL........
###### Alternative approach: make one giant dataframe with all genes and cluster designations, then plot facets to show all clusters
# ... but does that work with the exemplars too?...


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

##### Cluster Profiler #####
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
# Going to need external annotation data
biocLite("org.Ce.eg.db")
library("clusterProfiler")
library("org.Ce.eg.db")
?org.Ce.eg()
?org.Ce.eg.db
org.Ce.eg_dbInfo()

columns(org.Ce.eg.db) # shows which kinds of data can be returned from org.Ce.eg.db
keytypes(org.Ce.eg.db) # in this case shows the same as columns()

gname <- org.Ce.egGENENAME
mapped_genes <- mappedkeys(gname)
tail(mapped_genes)
glist <- as.list(gname[mapped_genes])  
glist[1:5]  
glist[30000:30010]  

# AnnotationDbi list of Wormbase genes - Bimap
WBGnames <- org.Ce.egWORMBASE
mapped_WBgenes <- mappedkeys(WBGnames)
WBGnames.list <- as.list(WBGnames[mapped_WBgenes])
WBGnames.list[1:5]
head(WBGnames.list)
tail(WBGnames.list)
?Bimap
class(WBGnames)
show(WBGnames) # minimal infor about Bimap object
summary(WBGnames) # bit more info about Bimap object - L/Rkeynames

### see all the AnnDbi mappings for a given gene ###
bitr("WBGene00000695", fromType = "WORMBASE", toType = keytypes(org.Ce.eg.db), OrgDb = org.Ce.eg.db)

# Bimap-direction
?`Bimap-direction`
direction(WBGnames)
Lkeys(WBGnames)
Rkeys(WBGnames)
Llength(WBGnames) # 46734
Rlength(WBGnames) # 20303
mappedLkeys(WBGnames)
count.mappedLkeys(WBGnames)

### example data for compareCluster
data("gcSample")
lapply(gcSample, head)
str(gcSample)
gcSample$X1
str(gcSample$X1)
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))
dotplot(ck)

# trying to get the gene names from the clusters into a format that clusterProfiler can use...
apresults2@clusters[[11]]
c11 <- labels(apresults2[[11]])
c12 <- labels(apresults2[[12]])
?list
# names(c11) <- "c11"
my.clusters <- list(c11,c12)
my.clusters
head(c11)
str(c11)
summary(c11)
class(c11)

?APResult

# compareCluster(apresults2@clusters, fun = "enrichKEGG")

# create empty lists for the cluster gene names
my.clusters <- vector("list",(length(apresults2)))
my.clusters
# populate the slots in the empty list with each cluster's genes
for (i in 1:length(apresults2)) {
  my.clusters[[i]] <- labels(apresults2[[i]])
}
# workable list of clusters!
my.clusters
# export cluster gene names to a txt file, with a new line for each cluster
lapply(my.clusters, write, "APclustering_Results_20180503.txt", append=TRUE, ncolumns=1000)

# ComClust.groupGO <- compareCluster(my.clusters, fun = "groupGO", OrgDb="org.Ce.eg.db")


############################################################################################################################
############################################################################################################################
################################# GO analysis using ClusterProfiler ########################################################


## Test script from http://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
# sample_test <- enrichGO(sample_gene, OrgDb=Cgriseus, pvalueCutoff=1, qvalueCutoff=1)
# head(summary(sample_test))

my.enrichGO <- enrichGO(my.clusters, OrgDb = "org.Ce.eg.db", pvalueCutoff = 0.1)  # ERRORS: need to convert gene names?!?
# --> No gene can be mapped....
# --> Expected input gene ID: 13184808,176800,176860,176622,174606,176622
# --> return NULL...

### Test GO enrichment on smallest gene cluster (c11)
?bitr # Biological Id TranslatoR
length(my.clusters[[11]]) # 96
keytypes(org.Ce.eg.db) # I can find the ID types that I want
# Getting a data frame with various IDs in columns
c11eg <- bitr(my.clusters[[11]], fromType = "WORMBASE", toType = "ENTREZID", OrgDb = "org.Ce.eg.db")
c11eg # has 97 entries, 1 more than c11 due to multiple mappings of EntrezIDs to a WBGene
# Performing GO enrichment for the list of EntrezIDs corresponding to cluster 11
c11eg.enrichGO <- enrichGO(c11eg$ENTREZID, OrgDb = "org.Ce.eg.db")
# info on what the output includes
?`enrichResult-class`
c11eg.enrichGO # structure format, not very readable
as.data.frame(c11eg.enrichGO) # makes enrichment output in a readable table
# translate Entrez geneID into human readable gene names
c11eg.enrichGO.readable <- setReadable(c11eg.enrichGO, OrgDb = org.Ce.eg.db)
as.data.frame(c11eg.enrichGO.readable)
#
##
###
# Visualization of enrichGO results for single test cluster 11
# GO terms vs. GeneRatio with dot size ~ number of genes and dot color ~ p.adjust
dotplot(c11eg.enrichGO.readable)
# Map of interconnectivity of genes and enriched GO terms
cnetplot(c11eg.enrichGO.readable) # possible to color code by fold change!? or p-value..?
# Chart of enriched GO terms and parents in tree-like format
par(cex = .4) # text plotted is too small to read, presetting cex inversely proportional to plotted text size
plotGOgraph(c11eg.enrichGO.readable)
# new plot function?
emapplot(c11eg.enrichGO.readable)
?emapplot
###
##
#
##### Translating IDs in my.clusters to Entrez for compareCluster analysis #####
my.clusters[[1]]
# my.clusters.eg <- bitr(my.clusters, fromType = "WORMBASE", toType = "ENTREZID", OrgDb = "org.Ce.eg.db")
my.clusters.eg <- lapply(my.clusters, bitr, fromType = "WORMBASE", toType = "ENTREZID", OrgDb = org.Ce.eg.db)
my.clusters.eg
summary(my.clusters.eg)
my.clusters.eg[[1]]$ENTREZID
length(my.clusters.eg[[1]]$ENTREZID) # 616

my.clusters.eg_only <- vector("list",(length(my.clusters.eg)))
my.clusters.eg_only
# populate the slots in the empty list with each cluster's genes
for (i in 1:length(my.clusters.eg_only)) {
  my.clusters.eg_only[[i]] <- my.clusters.eg[[i]]$ENTREZID
}
# workable list of clusters!
my.clusters.eg_only
# List of clusters needs to be named for compareCluster
paste0("C",1:12)
names(my.clusters.eg_only) <- paste0("C",1:12)
# Running GO enrichment for all clusters
my.enrichGO <- compareCluster(my.clusters.eg_only, fun = "enrichGO", OrgDb = org.Ce.eg.db)
head(as.data.frame(my.enrichGO))
summary(my.enrichGO)
str(my.enrichGO)
dotplot(my.enrichGO.2, showCategory=NULL)
unique(my.enrichGO.2@compareClusterResult$ID)
my.enrichGO.2 <- simplify(my.enrichGO, cutoff=0.7, by="p.adjust", select_fun=min)

###########################################################################################################################
##### Formula Interface for compareCluster #####
#
## What are the enriched GO terms by timepoint and up/downreg?
# created a .csv in excel
my.df <- read.csv("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Results and Analysis/Clustering/WBID_HrsPHS_DE.csv")
head(my.df)
summary(my.df)
class(my.df)
my.df$ï..WBID
columns(org.Ce.eg.db)
# find mappings between WBID and EntrezID
my.df.eg <- bitr(my.df$ï..WBID, fromType="WORMBASE", toType="ENTREZID",OrgDb=org.Ce.eg.db)
head(my.df.eg)
names(my.df)[1] <- "WORMBASE"
names(my.df.eg)[1]
# merge data frames to add EntrezIDs
my.df.merge <- merge(my.df, my.df.eg, by="WORMBASE", sort=F)
head(my.df.merge)
summary(my.df.merge)
# run enrichment according to experimental conditions
my.df.enrichGO <- compareCluster(ENTREZID~HrsPHS+DE, data = my.df.merge, fun = "enrichGO", OrgDb=org.Ce.eg.db)
as.data.frame(my.df.enrichGO)
dotplot(my.df.enrichGO, showCategory=10)

dotplot(my.df.enrichGO.2, x=~HrsPHS) + ggplot2::facet_grid(~DE) + ggplot2::scale_x_discrete(limits=c("3","6","12","30"))

write.csv(as.data.frame(my.df.enrichGO), file = "myEnrichGO_HrsPHS_DE.csv")
#use simplify to remove redundant GO terms
# https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
my.df.enrichGO.3 <- simplify(my.df.enrichGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(my.df.enrichGO.3, showCategory=10)
str(my.df.enrichGO.2) #246 GO terms?
str(my.df.enrichGO) #584 GO terms?
str(my.df.enrichGO.3) #136 GO terms

##### Enrich KEGG #####
my.clusters
my.clusters.eg
my.clusters.eg_only
search_kegg_organism("cel",by="kegg_code")
keytypes(org.Ce.eg.db)

my.clusters.IDs <- lapply(my.clusters, bitr, fromType = "WORMBASE", toType = c("ENTREZID","REFSEQ","ALIAS","UNIPROT"), OrgDb = org.Ce.eg.db)
my.clusters.IDs[[1]]

my.clusters.uniprot_only <- vector("list",(length(my.clusters.IDs)))
my.clusters.uniprot_only
# populate the slots in the empty list with each cluster's genes
for (i in 1:length(my.clusters.uniprot_only)) {
  my.clusters.uniprot_only[[i]] <- my.clusters.IDs[[i]]$UNIPROT
}
my.clusters.uniprot_only[[11]]

my.clusters.kegg <- lapply(my.clusters.uniprot_only, bitr_kegg, fromType = "uniprot", toType = "kegg", organism = "cel")
# Warning messages:
#   1: In FUN(X[[i]], ...) : 28.17% of input gene IDs are fail to map...
# 2: In FUN(X[[i]], ...) : 28.81% of input gene IDs are fail to map...
# 3: In FUN(X[[i]], ...) : 21.26% of input gene IDs are fail to map...
# 4: In FUN(X[[i]], ...) : 14.49% of input gene IDs are fail to map...
# 5: In FUN(X[[i]], ...) : 14.65% of input gene IDs are fail to map...
# 6: In FUN(X[[i]], ...) : 16.31% of input gene IDs are fail to map...
# 7: In FUN(X[[i]], ...) : 43.97% of input gene IDs are fail to map...
# 8: In FUN(X[[i]], ...) : 29.92% of input gene IDs are fail to map...
# 9: In FUN(X[[i]], ...) : 21.85% of input gene IDs are fail to map...
# 10: In FUN(X[[i]], ...) : 11.5% of input gene IDs are fail to map...
# 11: In FUN(X[[i]], ...) : 14.62% of input gene IDs are fail to map...
# 12: In FUN(X[[i]], ...) : 20.5% of input gene IDs are fail to map...
my.clusters.kegg
my.clusters.kegg_only <- vector("list",(length(my.clusters.IDs)))
my.clusters.kegg_only
# populate the slots in the empty list with each cluster's genes
for (i in 1:length(my.clusters.kegg_only)) {
  my.clusters.kegg_only[[i]] <- my.clusters.kegg[[i]]$kegg
}
my.clusters.kegg_only
# List of clusters needs to be named for compareCluster
names(my.clusters.kegg_only) <- paste0("C",1:12)
### KEGG enrichment for gene clusters ###
my.enrichKEGG <- compareCluster(my.clusters.kegg_only, fun = "enrichKEGG", organism="cel")
as.data.frame(my.enrichKEGG)
dotplot(my.enrichKEGG)
## there are hardly any surviving genes used in the enrichment, but some potentially interesting hits
## cluster10 in particular will be used for additional analysis
c10.kegg <- enrichKEGG(my.clusters.kegg_only[[10]], organism = "cel")
as.data.frame(c10.kegg)
dotplot(c10.kegg)
cnetplot(c10.kegg)
browseKEGG(c10.kegg, pathID = c10.kegg@result$ID[3])
# biocLite("pathview")
library("pathview")
pathview(gene.data = my.clusters.eg_only[[10]], pathway.id = c10.kegg@result$ID[1:2], species = "cel")

##### KEGG clustering with a formula interface #####
# need to translate IDs to KEGG
head(my.df)
head(my.df.merge)
my.df.bitr.kegg <-  bitr_kegg(my.df.merge$ENTREZID, fromType = "ncbi-geneid", toType = "kegg", organism = "cel")
# Warning message:
#   In bitr_kegg(my.df.merge$ENTREZID, fromType = "ncbi-geneid", toType = "kegg",  :
#                  0.04% of input gene IDs are fail to map...
head(my.df.bitr.kegg)
my.df.kegg.merge <- merge(my.df.merge, my.df.bitr.kegg, by.x="ENTREZID", by.y="ncbi-geneid")
head(my.df.kegg.merge)

my.df.enrichKEGG <- compareCluster(kegg~DE+HrsPHS, data = my.df.kegg.merge, fun="enrichKEGG", organism="cel")
as.data.frame(my.df.enrichKEGG)
dotplot(my.df.enrichKEGG)#, showCategory=NULL)

write.csv(my.df.enrichKEGG, file = "enrichKEGG_timecourse_20180507.csv")

pathview(gene.data = my.df.eg$ENTREZID, pathway.id = my.df.enrichKEGG@compareClusterResult$ID[146], species = "cel")



