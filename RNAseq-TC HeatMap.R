###Heat Map of temporal gene expression###
#
#
#install.packages("RColorBrewer")
library("RColorBrewer")

## create a sequential palette for usage and show colors
mypalette<-brewer.pal(7,"Greens")
image(1:7,1,as.matrix(1:7),col=mypalette,xlab="Greens (sequential)",
      ylab="",xaxt="n",yaxt="n",bty="n")
## display a divergent palette
display.brewer.pal(11,"Spectral")

display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE,
                   colorblindFriendly=FALSE)

##test heatmap for linear scale data - too skewed for a good heatmap
heatmap.2(NormELT7matrix[80:100,], Colv = F, dendrogram = "row", colsep = c(3,6,9), trace = "none", scale = "row")

#labels for the twelve columns
heatnames <- c("0hr_1","0hr_2","0hr_3","6hr_1","6hr_2","6hr_3","12hr_1","12hr_2","12hr_3","24hr_1","24hr_2","24hr_3")

#log2 gene expression heat map
logNormELT7 <- log2(NormELT7matrix+1)
head(logNormELT7)
heatmap.2(logNormELT7[2000:3000,], Colv = F, Rowv = T, dendrogram = "none", colsep = c(3,6,9),
          trace = "none", labRow = NA, labCol = heatnames, adjCol = 0.8, col = rev(brewer.pal(9, "BuGn")),
          key.title = NA, key.xlab = "log2 Gene Expression", density.info = "density")

###full log2 expression heat map - LONG RUN TIME
heatmap.2(logNormELT7, Colv = F, Rowv = T, dendrogram = "none", colsep = c(3,6,9),
          trace = "none", labRow = NA, labCol = heatnames, adjCol = 0.8, col = rev(brewer.pal(9, "PuBu")),
          key.title = NA, key.xlab = "log2 Gene Expression", density.info = "density")
###

###Duplicate for control data (different color heatmap)
logNormGFP <- log2(NormGFPmatrix+1)
head(logNormGFP)

heatmap.2(logNormGFP, Colv = F, Rowv = T, dendrogram = "none", colsep = c(3,6,9),
          trace = "none", labRow = NA, labCol = heatnames, adjCol = 0.8, col = rev(brewer.pal(9, "BuGn")),
          key.title = NA, key.xlab = "log2 Gene Expression", density.info = "density")
###