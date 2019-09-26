# Tissue Gene Expression Changes
library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library("doBy")
# Import lists of Tissue Expressed Genes from Spreadsheet
library("readxl")
Tissues <- read_excel("~/UCSB/Rothman Lab/Data/Published Data/Wormbase Tissue Ontology/Wormbase Tissue Ontology.xlsx", sheet = "Tissues")
View(Tissues)

# Tissue.list <- vector("list", length = length(colnames(Tissues)))
Tissue.list <- as.list(Tissues)
Tissue.list <- lapply(Tissue.list, function(x) x[!is.na(x)])
str(Tissue.list)

Tissue.list$`coelomocyte (WBbt:0005751) (894 direct)`
Tissue.list[[1]]
names(Tissue.list)[1]

##### Tissue gene list plotting #####
Tissue.list[[1]] %in% names(ddsTC)

goi <- Tissue.list[[1]][Tissue.list[[1]] %in% names(ddsTC)]
stopifnot(all((goi) %in% names(ddsTC)))

tcounts <- (scale(t(counts(ddsTC[goi], normalized=TRUE, replaced=FALSE)))) %>%
  merge(colData(ddsTC), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

summary(tcounts)
head(tcounts)
tcounts.sum <- summaryBy(expression ~ Strain + HrsPHS, data = tcounts, FUN = c(length,mean,sd))
tcounts.sum

ggplot(tcounts, aes(x = HrsPHS, y = expression, group = interaction(gene, Strain), color = Strain)) + 
  geom_line(stat = "smooth", method = "loess", alpha = .5) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_grid(Strain ~ .) +
  labs(x="Hours PHS", 
       y="Expression (scaled normalized counts)", 
       fill="Strain", 
       title=names(Tissue.list)[1]) +
  geom_line(data = tcounts.sum, aes(x = HrsPHS, y = expression.mean, group = Strain), color = "black")

###
names(Tissue.list[1])

##### Custom function to plot tissue specific gene expression changes! #####
TissuePlot <- function(x) {
  goi <- x[x %in% names(ddsTC)]
  
  tcounts <- scale(t(counts(ddsTC[goi], normalized = T, replaced = F))) %>%
    merge(colData(ddsTC), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
  
  tcounts.sum <- summaryBy(expression ~ Strain + HrsPHS, data = tcounts, FUN = c(length,mean,sd))
  
  p <- ggplot(tcounts, aes(x = HrsPHS, y = expression, group = interaction(gene, Strain), color = Strain)) + 
    geom_line(stat = "smooth", method = "loess", alpha = .4) +
    theme_minimal() +
    theme(legend.position = "none") +
    facet_grid(Strain ~ .) +
    labs(x="Hours PHS", 
         y="Expression (scaled normalized TPM)", 
         fill="Strain", 
         title=names(x)) +
    geom_line(data = tcounts.sum, aes(x = HrsPHS, y = expression.mean, group = Strain), color = "black")
  
  return(p)
  
}

names(Tissue.list)

hypodermal <- TissuePlot(Tissue.list$`hypodermal cell (WBbt:0007846) (38 direct)`)
hypodermal + ggtitle(as.character(names(Tissue.list[3])))

tissue2 <- TissuePlot(Tissue.list[[2]])
tissue2 + ggtitle(as.character(names(Tissue.list[2])))


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
##########################################################################################################


### Loop to create ggplot data for each tissue
# Container for TissuePlot data
AllTissuePlots <- vector("list",length = length(Tissue.list))
# Looping TissuePlot through list of tissue expressed genes
for (i in 1:length(AllTissuePlots)) {
  tp <- TissuePlot(Tissue.list[[i]])
  tp2 <- tp + ggtitle(as.character(names(Tissue.list[i])))
  AllTissuePlots[[i]] <- tp2
  
}

# Test plotting tissue 1
AllTissuePlots[[1]]
# ggsave(paste0("tissue",1,".tiff"), plot = AllTissuePlots[[1]], width = 5, height = 5)

### Exporting image files for all tissues ###
for (i in 1:16) {
  ggsave(paste0("tissue",i,".tiff"), plot = AllTissuePlots[[i]], width = 4, height = 4)
}









######################################################################################################################################

######################################################################################################################################

######################## E N R I C H E D . T I S S U E S . . S P E N C E R . 2 0 1 1 #################################################

######################################################################################################################################

######################################################################################################################################

##### Lists of Selectively Enriched Genes sets for various tissues - from Spencer et al. 2011 #####

Enrich.Tissues <- read_excel("~/UCSB/Rothman Lab/Data/Published Data/Spencer et al 2011/Selectively_enriched_genes/Selectively_enriched_genes.xlsx")
View(Enrich.Tissues)

Enrich.Tissues.list <- as.list(Enrich.Tissues)
Enrich.Tissues.list <- lapply(Enrich.Tissues.list, function(x) x[!is.na(x)])
str(Enrich.Tissues.list)

### Loop to create ggplot data for each tissue
# Container for TissuePlot data
AllEnriched <- vector("list",length = length(Enrich.Tissues.list))
# Looping TissuePlot through list of tissue expressed genes
for (i in 1:length(AllEnriched)) {
  tp <- TissuePlot(Enrich.Tissues.list[[i]])
  tp2 <- tp + ggtitle(paste0(as.character(names(Enrich.Tissues.list[i]))," (",length(Enrich.Tissues.list[[i]])," genes)"))
  AllEnriched[[i]] <- tp2
  
}
# Test plot first tissue
AllEnriched[[6]]
### Exporting image files for all tissues ###
for (i in 1:length(AllEnriched)) {
  ggsave(paste0("EnrichedTissue",i,".tiff"), plot = AllEnriched[[i]], width = 4, height = 4)
}









##########################################################################################################################################

######################################################################################################################################

################################### F O L D  . C H A N G E .  P L O T S ##############################################################

######################################################################################################################################

######################################################################################################################################




##### Plotting enriched genes sets by Fold Change #####

# need to create master table of WBGenes by FC for each timepoint....

FC.df <- as.data.frame(rownames(res3hr))
names(FC.df) <- "Gene"
FC.df$FC.3 <- res3hr$log2FoldChange
FC.df$FC.6 <- res6hr$log2FoldChange
FC.df$FC.12 <- res12hr$log2FoldChange
FC.df$FC.30 <- res30hr$log2FoldChange
head(FC.df)
# now need to convert from wide data to long data
# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)
# data_long <- gather(olddata_wide, condition, measurement, control:cond2, factor_key=TRUE)
FC.long <- gather(data = FC.df, key = HrsPHS, value = FC, FC.3:FC.30, factor_key = TRUE)
FC.long <- orderBy(~ Gene, FC.long)
FC.long[1:10,]
# Get the average FC for a subset of genes
FC.avg <- summaryBy(FC ~ HrsPHS, data = FC.long[1:4000,], FUN = c(length,mean,sd))
FC.avg
# Plot of FC for genes of intrest, avg FC with labels
ggplot(FC.long[1:4000,], aes(x = HrsPHS, y = FC, group = Gene)) +
  theme_classic() +
  geom_line(alpha = 0.4, color = "grey", size = 1.2) +
  geom_line(data = FC.avg, inherit.aes = F, aes(x = HrsPHS, y = FC.mean, group = NA), color = "black", size = 1.4) +
  geom_text(aes(x = HrsPHS, y = FC.mean, label = round(FC.mean, digits = 3)), inherit.aes = F, data = FC.avg, nudge_y = 0.8) +
  ggtitle("Genes of Intrest")

### Above pipeline works for plotting FC for groups of genes

Enrich.Tissues.list[17]

FC.enrich.17 <- FC.long[FC.long$Gene %in% Enrich.Tissues.list[[17]],]
head(FC.enrich.17)
FC.avg <- summaryBy(FC ~ HrsPHS, data = FC.enrich.17, FUN = c(length,mean,sd))
FC.avg

ggplot(FC.enrich.17, aes(x = HrsPHS, y = FC, group = Gene)) +
  theme_classic() +
  coord_cartesian(ylim = c(-2,2)) +
  # geom_line(alpha = 0.7, color = "grey", size = 1.2) +
  geom_line(data = FC.avg, inherit.aes = F, aes(x = HrsPHS, y = FC.mean, group = NA), color = "black", size = 1.4) +
  geom_text(aes(x = HrsPHS, y = FC.mean, label = round(FC.mean, digits = 3)), inherit.aes = F, data = FC.avg, nudge_y = 0.5) +
  labs(x="Hours PHS", 
     y="log2 Fold Change", 
     title=names(Enrich.Tissues.list[17])) +
  scale_x_discrete(labels = c("3","6","12","30"))

## plot just the mean ##
ggplot(FC.avg, aes(x = HrsPHS, y = FC.mean, group = NA)) +
  theme_classic() +
  geom_line(size = 1.4) +
  geom_errorbar(aes(ymin = FC.mean-FC.sd, ymax = FC.mean+FC.sd), width = 0.2) +
  labs(x="Hours PHS", 
       y="log2 Fold Change", 
       title=names(Enrich.Tissues.list[17])) +
  scale_x_discrete(labels = c("3","6","12","30"))



#### custom function: input gene list, output ggplot graph data ####
FoldChangeTissuePlot <- function(x, title){
  FCdata <- FC.long[FC.long$Gene %in% x,]
  
  FCdata.avg <- summaryBy(FC ~ HrsPHS, data = FCdata, FUN = mean)
  
  p <- ggplot(FCdata, aes(x = HrsPHS, y = FC, group = Gene)) +
    theme_classic() +
    geom_line(alpha = 0.4, color = "grey", size = 1.2) +
    geom_line(data = FCdata.avg, inherit.aes = F, aes(x = HrsPHS, y = FC.mean, group = NA), color = "black", size = 1.4) +
    geom_text(aes(x = HrsPHS, y = FC.mean, label = round(FC.mean, digits = 3)), inherit.aes = F, data = FCdata.avg, nudge_y = 1) +
    labs(x="Hours PHS", 
         y="log2 Fold Change", 
         title=title) +
    scale_x_discrete(labels = c("3","6","12","30"))
  
  return(p)
}
# # # # # # # # end function # # # # # # # #

tissue1 <- FoldChangeTissuePlot(Enrich.Tissues.list[[1]], title = names(Enrich.Tissues.list[1]))
tissue1
# it works!


################################################################################
################################################################################
########## 2019.09.09 new code - all genes FC plot for pub #####################

allgenesFC <- FoldChangeTissuePlot(names(ddsTC), title = "All Genes (20094)")

allgenesFC # doesn't plot the mean for some reason...
            # 13 genes have NA fold changes - same as had NaN expression data

FCdata.avg <- summaryBy(FC ~ HrsPHS, data = FC.long[!is.na(FC.long$FC),], FUN = mean) # subsetting FC.long to remove genes with NA
FCdata.avg


p <- ggplot(FC.long, aes(x = HrsPHS, y = FC, group = Gene)) +
  theme_classic() +
  geom_line(alpha = 0.2, color = "grey", size = 1.2) +
  geom_line(data = FCdata.avg, inherit.aes = F, aes(x = HrsPHS, y = FC.mean, group = NA), color = "black", size = 1.4) +
  geom_text(aes(x = HrsPHS, y = FC.mean, label = round(FC.mean, digits = 3)), inherit.aes = F, data = FCdata.avg, nudge_y = 1) +
  labs(x="Hours PHS", 
       y="log2 Fold Change", 
       title="All Genes (20094)") +
  scale_x_discrete(labels = c("3","6","12","30"))

p

ggsave(filename = "All Genes FC plot.tiff", width = 4, height = 4)

#####################################################################################
#####################################################################################
#####################################################################################



### loop for all tissue lists

# Container for TissuePlot data
AllEnrichedFC <- vector("list",length = length(Enrich.Tissues.list))
# Looping FoldChangeTissuePlot through list of tissue expressed genes
for (i in 1:length(AllEnrichedFC)) {
 
    AllEnrichedFC[[i]] <- FoldChangeTissuePlot(Enrich.Tissues.list[[i]], title = names(Enrich.Tissues.list[i]))
  
}
# Test plot first tissue
AllEnrichedFC[[6]]
### Exporting image files for all tissues ###
for (i in 1:length(AllEnrichedFC)) {
  ggsave(paste0("EnrichedTissue",i,".FC.tiff"), plot = AllEnrichedFC[[i]], width = 4, height = 4)
}


##########################################################################################################
##########################################################################################################

## separate plots for all genes in a tissue, and just the means

GeneFCSummary <- function(x){
  FCdata <- FC.long[FC.long$Gene %in% x,]
  FCdata.avg <- summaryBy(FC ~ HrsPHS, data = FCdata, FUN = c(length,mean,sd))
  return(FCdata.avg)
}

# test of new function: GeneFCSummary
tissue1.summary <- GeneFCSummary(Enrich.Tissues.list[[1]])
tissue1.summary

Tissue.Summaries <- vector("list", length = length(Enrich.Tissues.list))
for (i in 1:length(Tissue.Summaries)) {
  Tissue.Summaries[[i]] <- GeneFCSummary(Enrich.Tissues.list[[i]])
}
names(Tissue.Summaries) <- names(Enrich.Tissues.list)
Tissue.Summaries[[2]]$FC.length[1]

# tissues of intrest to plot
toi <- c(2,3,4,7,12,13,17,18,19,20,21,24)


## plot just the mean ##
FCmean.plots <- vector("list", length = length(toi))
for (i in toi) {
p <- ggplot(Tissue.Summaries[[i]], aes(x = HrsPHS, y = FC.mean, group = NA)) +
  # theme_classic() +
  # coord_cartesian(ylim = c(-1.5,1.5)) +
  geom_line(size = 2.5) +
  geom_errorbar(aes(ymin = FC.mean-FC.sd, ymax = FC.mean+FC.sd), width = 0.2) +
  labs(x="Hours PHS", 
       y="log2 Fold Change", 
       title=paste0(names(Tissue.Summaries[i]), " (", Tissue.Summaries[[i]]$FC.length[1], " genes)")) +
  scale_x_discrete(labels = c("3","6","12","30"))

FCmean.plots[[i]] <- p
}
# test ggplot loop output
FCmean.plots[[12]]

# save plots
for (i in 1:length(FCmean.plots)) {
  ggsave(paste0("EnrichedTissue",i,".FC.means.SD.tiff"), plot = FCmean.plots[[i]], width = 3.5, height = 3.5)
  
}



######################################################################################################

# multi gene plots without means


FoldChangeTissuePlot2 <- function(x, title){
  FCdata <- FC.long[FC.long$Gene %in% x,]
  
  # FCdata.avg <- summaryBy(FC ~ HrsPHS, data = FCdata, FUN = mean)
  
  p <- ggplot(FCdata, aes(x = HrsPHS, y = FC, group = Gene)) +
    theme_classic() +
    geom_line(alpha = 0.7, color = "grey", size = 1.2) +
    # geom_line(data = FCdata.avg, inherit.aes = F, aes(x = HrsPHS, y = FC.mean, group = NA), color = "black", size = 1.4) +
    # geom_text(aes(x = HrsPHS, y = FC.mean, label = round(FC.mean, digits = 3)), inherit.aes = F, data = FCdata.avg, nudge_y = 0.5) +
    labs(x="Hours PHS", 
         y="log2 Fold Change", 
         title=title) +
    scale_x_discrete(labels = c("3","6","12","30"))
  
  return(p)
}
# # # # # # # # end function # # # # # # # #

tissue1 <- FoldChangeTissuePlot2(Enrich.Tissues.list[[1]], title = paste0(names(Enrich.Tissues.list[1]), " (", length(Enrich.Tissues.list[[1]]), " genes)"))
tissue1

### loop for all tissue lists

# Container for TissuePlot data
AllEnrichedFC2 <- vector("list",length = length(toi))
# Looping FoldChangeTissuePlot through list of tissue expressed genes
for (i in toi) {
  
  AllEnrichedFC2[[i]] <- FoldChangeTissuePlot2(Enrich.Tissues.list[[i]],
                                               title = paste0(names(Enrich.Tissues.list[i]), " (", length(Enrich.Tissues.list[[i]]), " genes)"))
  
}
# Test plot first tissue
AllEnrichedFC2[[7]]
### Exporting image files for all tissues ###
for (i in toi) {
  ggsave(paste0("EnrichedTissue",i,".FC.noMean.tiff"), plot = AllEnrichedFC2[[i]], width = 3.5, height = 3.5)
}



