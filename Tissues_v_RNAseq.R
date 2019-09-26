##### Tissue specific (enriched) gene expression vs RNAseq differential gene expression #####

# Goal: data table/heatmap of significant enrichment between tissue type gene lists and DEG lists

# 1) need DEG lists
head(RNAseq.DEGs)
head(Up.3)
length(Up.3) #vecotr of 4231 genes

# 2) need lists of tissues
# I'm using the same lists of selectively enriched gene sets as the gene expression plots previously made, from Spencer et al. 2011

selectively_enriched <- read.csv("Selectively_enriched_genes.csv", stringsAsFactors = F, na.strings = "")
View(selectively_enriched)
names(selectively_enriched)[1] <- "EE.BAG.neuron"
names(selectively_enriched)

tissue.list <- as.list(selectively_enriched)
tissue.list <- lapply(tissue.list, function(x) x[!is.na(x)])

tissue.counts <- data.frame(names(selectively_enriched))
for (i in 1:length(tissue.list)) {
  tissue.counts$num.genes[i] <- length(tissue.list[[i]])
}

tissue.counts
# make a separate data frame for each up/downreg gene set and numbers of overlapping genes
Up.3.tissues <- tissue.counts
Up.3.tissues

for (i in 1:length(tissue.list)) {
  tmp <- tissue.list[[i]]
  overlap <- tmp[tmp %in% Up.3]
  Up.3.tissues$num.overlap[i] <- length(overlap)
  Up.3.tissues$num.tissue.notDEG[i] <- Up.3.tissues$num.genes[i] - length(overlap)
  Up.3.tissues$num.NOTtissue.DEG[i] <- length(Up.3) - length(overlap)
  Up.3.tissues$num.NOTtissue.notDEG[i] <- length(RNAseq.DEGs$WBGene) - length(Up.3) - Up.3.tissues$num.genes[i] + length(overlap)
  Up.3.tissues$pc.tissue[i] <- length(overlap)/length(tmp)
  Up.3.tissues$pc.DEG[i] <- length(overlap)/length(Up.3)
  ct <- matrix(c(Up.3.tissues$num.overlap[i],Up.3.tissues$num.NOTtissue.DEG[i],
                 Up.3.tissues$num.tissue.notDEG[i],Up.3.tissues$num.NOTtissue.notDEG[i]), nrow = 2)
  ct.ft <- fisher.test(ct)
  Up.3.tissues$p.value[i] <- ct.ft$p.value
}

Up.3.tissues


##### compiling above code into easily executable function #####

gene.lists.overlap.func <- function(gene.list, DEGs){
  return.df <- tissue.counts
  for (i in 1:length(gene.list)) {
    tmp <- gene.list[[i]]
    overlap <- tmp[tmp %in% DEGs]
    return.df$num.overlap[i] <- length(overlap)
    return.df$num.tissue.notDEG[i] <- return.df$num.genes[i] - length(overlap)
    return.df$num.NOTtissue.DEG[i] <- length(DEGs) - length(overlap)
    return.df$num.NOTtissue.notDEG[i] <- length(RNAseq.DEGs$WBGene) - length(DEGs) - return.df$num.genes[i] + length(overlap)
    return.df$pc.tissue[i] <- length(overlap)/length(tmp)
    return.df$pc.DEG[i] <- length(overlap)/length(DEGs)
    ct <- matrix(c(return.df$num.overlap[i],return.df$num.NOTtissue.DEG[i],
                   return.df$num.tissue.notDEG[i], return.df$num.NOTtissue.notDEG[i]), nrow = 2)
    ct.ft <- fisher.test(ct)
    return.df$p.value[i] <- ct.ft$p.value
    }
  return(return.df)
}

#############################################################################
Up.3.tissues.func.result <- gene.lists.overlap.func(gene.list = tissue.list, DEGs = Up.3)
Up.3.tissues.func.result
Up.3.tissues == Up.3.tissues.func.result # ALL TRUE "beautiful" ;)

Up.6.tissues <- gene.lists.overlap.func(tissue.list, Up.6)
Up.12.tissues <- gene.lists.overlap.func(tissue.list, Up.12)
Up.30.tissues <- gene.lists.overlap.func(tissue.list, Up.30)

Down.3.tissues <- gene.lists.overlap.func(tissue.list, Down.3)
Down.6.tissues <- gene.lists.overlap.func(tissue.list, Down.6)
Down.12.tissues <- gene.lists.overlap.func(tissue.list, Down.12)
Down.30.tissues <- gene.lists.overlap.func(tissue.list, Down.30)

#############################################################################################

####################### P - V A L U E  H E A T  M A P #######################################

#############################################################################################

## Need to create a matrix of data for heatmap. Can also use `as.matrix()` on a data frame.
library(gplots)
# help(heatmap.2)

name.heatmap.cols <- c("Up 3","Up 6","Up 12","Up 30","Down 3","Down 6","Down 12","Down 30")


pvals.matrix.tissue <- data.frame(Up.3.tissues$p.value, Up.6.tissues$p.value, Up.12.tissues$p.value, Up.30.tissues$p.value,
                                  Down.3.tissues$p.value, Down.6.tissues$p.value, Down.12.tissues$p.value, Down.30.tissues$p.value)
rownames(pvals.matrix.tissue) <- tissue.counts$names.selectively_enriched.
colnames(pvals.matrix.tissue) <- name.heatmap.cols
pvals.matrix.tissue <- as.matrix(pvals.matrix.tissue)
is.matrix(pvals.matrix.tissue)
pvals.matrix.tissue
heatmap.2(pvals.matrix.tissue)

write.csv(pvals.matrix.tissue, "Tissues.DEGs.pvalues.csv")

pvals.matrix.tissue.log <- -log(pvals.matrix.tissue)
is.infinite(pvals.matrix.tissue.log) # all FALSE
heatmap.2(pvals.matrix.tissue.log)

heatmap.2(pvals.matrix.tissue.log, trace = "none", dendrogram = "row", Colv = F, col = colorpanel(25,"black","green","white"),
          scale = "none", key = T, keysize = 1.1, density.info = "none", key.title = NA, key.xlab = "-log(p-value)",
          colsep = 4, margins = c(6,10))


