test <- data.frame(cbind(elt7geneTPM_5tp[2:3],elt7geneTPM_5tp[12:15],elt7geneTPM_5tp[4:6]), row.names = elt7geneTPM_5tp$WBGene)

test <- data.frame(cbind(elt7geneTPM_5tp[2:4],elt7geneTPM_5tp[14:17],elt7geneTPM_5tp[5:13]), row.names = elt7geneTPM_5tp$WBGene)

class(test)


x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
# compute the list mean for each list element
lapply(x, mean)
# median and quartiles for each list element
lapply(x, quantile, probs = 1:3/4)
sapply(x, quantile)
i39 <- sapply(3:9, seq) # list of vectors
sapply(i39, fivenum)
vapply(i39, fivenum,
       c(Min. = 0, "1st Qu." = 0, Median = 0, "3rd Qu." = 0, Max. = 0))

###############################################################################################################################################

# We can now plot the log2 fold changes in a heatmap (figure below).
topGenes <- head(order(resTC$padj),30)
mat <- betas[topGenes, -c(1,2)]
thr <- 5
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, gaps_col = 4)

write.csv(rownames(resTC[topGenes,]), file = "top30genes_log2FC heatmap.csv") #correct gene names, but not in the cluster order in the fig

# map.data <- pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
#          cluster_col=FALSE, gaps_col = 4, silent = T)

#took WBGene names from .csv > wormbase.org simple mine to get common names > excel copypasta transpose text join "" to paste list below...
mat.rownames <- c("rmd-2","fbxc-29","tmc-1","F36F2.1","R12C12.7","jph-1","hrpf-1","Y71A12B.12","sorb-1","dsbn-1","ddi-1","K02B12.2",
                  "C27A12.6","F59B10.3","C08B6.3","rde-1","unc-37","C27H6.8","B0252.3","nhr-17","unc-49","W04A8.5","clec-49",
                  "Y39A1A.16","nas-4","F09B9.4","vps-24","C35E7.5","pes-10","clec-196")

rownames(mat) <- mat.rownames
thr <- 5
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, gaps_col = 4)


# We can now plot the log2 fold changes in a heatmap (figure below).
topGenes <- head(order(resTC$padj),300)
mat <- betas[topGenes, -c(1,2)]
thr <- 5
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, gaps_col = 4, show_rownames = F)

# We can now plot the log2 fold changes in a heatmap (figure below).
# topGenes <- head(order(resTC$padj),300)
# topGenes <- resTC[resTC$lfcSE > 0,]
mat <- betas[, -c(1,2)]
mat <- mat[complete.cases(mat),] #need to remove any rows with NA
summary(mat)

thr <- 3 #threshold level for the range of FC values to plot
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
summary(mat)
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, gaps_col = 4, show_rownames = F)

# with kmeans_k clustering
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, gaps_col = 4, show_rownames = F, kmeans_k = 100)


################################################################################################
# ggplot2 cheat sheet
f <- ggplot(mpg, aes(class, hwy))
f + geom_violin(scale = "area")
e <- ggplot(mpg, aes(cty, hwy))
e + geom_quantile()
e + geom_point(aes(color = class))
e + geom_rug()
e + geom_text(aes(label = cty))
e + geom_smooth(method = lm)
e + stat_quantile

ggplot(mtcars) +
  stat_qq(geom = "line", aes(sample = mpg))
