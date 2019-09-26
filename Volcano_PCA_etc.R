head(resTC)

df.resTC <- data.frame(resTC)
head(df.resTC)


##### Volcano Plot #####
ggplot(df.resTC[df.resTC$padj<0.01,], aes(log2FoldChange, -log10(padj))) +
  geom_point(alpha = 0.1) + 
  scale_y_continuous(limits = c(0,100))

## 3 hrs
df.res3 <- data.frame(res3hr)

ggplot(df.res3[df.res3$padj<0.01,], aes(log2FoldChange, -log10(padj))) +
  geom_point(alpha = 0.1) + 
  scale_y_continuous(limits = c(0,100)) +
  ggtitle("3 hours PHS")

##6 hrs
df.res6 <- data.frame(res6hr)

ggplot(df.res6[df.res6$padj<0.01,], aes(log2FoldChange, -log10(padj))) +
  geom_point(alpha = 0.1) + 
  scale_y_continuous(limits = c(0,100)) +
  ggtitle("6 hours PHS")

##12 hrs
df.res12 <- data.frame(res12hr)

ggplot(df.res12[df.res12$padj<0.01,], aes(log2FoldChange, -log10(padj))) +
  geom_point(alpha = 0.1) + 
  scale_y_continuous(limits = c(0,100)) +
  ggtitle("12 hours PHS")

##30 hrs
df.res30 <- data.frame(res30hr)

ggplot(df.res30[df.res30$padj<0.01,], aes(log2FoldChange, -log10(padj))) +
  geom_point(alpha = 0.1) + 
  scale_y_continuous(limits = c(0,100)) +
  ggtitle("30 hours PHS")


##################################################################################################################################
# https://github.com/mikelove/DESeq2/blob/master/R/plots.R
vignette("DESeq2")


# 5.4 PCA plot
plotPCA(rld, intgroup = c("Strain", "HrsPHS"))

### Code for PCA plot in DESeq2 ###
plotPCA.DESeqTransform = function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
}



###############################################################################################
## Custom PCA analysis - modifying above code #####

# calculate the variance for each gene
rv <- rowVars(assay(rld))
head(rv)
summary(rv)

# select the ntop genes by variance
ntop = 500
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
head(select)

assay(rld)[10694,]

plot(rv)
plot(rv[select])

#### which genes have the greatest variance?
var.genes <- rownames(assay(rld)[select,])
write.csv(var.genes, file = "MostVariableGenes.csv")

# perform a PCA on the data in assay(x) for the selected genes
?prcomp
pca <- prcomp(t(assay(rld)[select,]))
pca$x
plot(pca) #"scree" plot ?
screeplot(pca)
summary(pca)
biplot(pca)

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
plot(percentVar)

pca$rotation[,1]

## some basic PCA plots looking at different components

plot(pca$x)
plot(pca$x[,2:3])
plot(pca$x[,3:4])

# install.packages("factoextra")
library("factoextra")

get_eigenvalue(pca)
sum(get_eigenvalue(pca)$eigenvalue)

fviz_eig(pca, addlabels = TRUE, ylim = c(0, 70))
fviz_screeplot(pca, addlabels = T, ggtheme = theme_classic())

var.pca <- get_pca_var(pca)
var.pca
# Coordinates
head(var.pca$coord)
# Cos2: quality on the factore map
head(var.pca$cos2)
# Contributions to the principal components
head(var.pca$contrib)

# To plot PCA variables
fviz_pca_var(pca)
select.pca.var <- rownames(pca$rotation[1:50,])
fviz_pca_var(pca, select.var = list(name = select.pca.var), repel = T)

# Color by cos2 values: quality on the factor map
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             # repel = TRUE # Avoid text overlapping
)

# Change the transparency by cos2 values
fviz_pca_var(pca, alpha.var = "cos2")


# You can visualize the cos2 of variables on all the dimensions using the corrplot package:
# install.packages("corrplot")  
library("corrplot")
# corrplot(var.pca$cos2, is.corr=FALSE) # default is way overplotted
corrplot(var.pca$cos2[1:30,], is.corr = F)
corrplot(var.pca$cos2[31:60,], is.corr = F)
corrplot(var.pca$cos2[41:60,1:5], is.corr = F)

# It’s also possible to create a bar plot of variables cos2 using the function fviz_cos2()[in factoextra]:
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(pca, choice = "var", axes = 1)
fviz_cos2(pca, choice = "var", axes = 1, top = 100, sort.val = "desc")



# The contribution of variables can be extracted as follow :
head(var.pca$contrib, 4)
# The contributions of variables in accounting for the variability in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 (i.e., Dim.1) and PC2 (i.e., Dim.2) are the most important
#   in explaining the variability in the data set.
# Variables that do not correlated with any PC or correlated with the last dimensions are variables
#   with low contribution and might be removed to simplify the overall analysis.

# It’s possible to use the function corrplot() [corrplot package] to highlight the most contributing variables for each dimension:
corrplot(var.pca$contrib[1:30,], is.corr=FALSE) 
corrplot(var.pca$contrib[31:60,], is.corr=FALSE) 
corrplot(var.pca$contrib[61:90,], is.corr=FALSE) 
corrplot(var.pca$contrib[91:120,], is.corr=FALSE) 
# Contributions of variables to PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 500)
pca.contrib <- fviz_contrib(pca, choice = "var", axes = 1, top = 250)
pca.contrib.data <- pca.contrib$data
order(pca.contrib.data$contrib, decreasing = TRUE)
sort(pca.contrib.data$contrib, decreasing = TRUE)
pca.contrib.data <- pca.contrib.data[order(pca.contrib.data$contrib, decreasing = TRUE),]
head(pca.contrib.data)
plot(pca.contrib.data)
write.csv(pca.contrib.data, file = "PCA dim1 gene contributions.csv")


# The most important (or, contributing) variables can be highlighted on the correlation plot as follow:
fviz_pca_var(pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))



### **** Below code does not work with the packages / functions used above ****
#
# the function dimdesc() [in FactoMineR], for dimension description, can be used to identify the most
# significantly associated variables with a given principal component.
# install.packages("FactoMineR")  
library("FactoMineR")
# pca.res.desc <- dimdesc(pca, axes = c(1,2), proba = 0.05)
pca.res.desc <- dimdesc(pca)
# Description of dimension 1
pca.res.desc$Dim.1
