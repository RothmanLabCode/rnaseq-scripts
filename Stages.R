##### Stage specific expressed genes, differential expression. Gene lists from Spencer et al. 2011 #####

library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(doBy)
library(readxl)
# Import lists of Tissue Expressed Genes from Spreadsheet
Stages <- read_excel("~/UCSB/Rothman Lab/Data/Published Data/Spencer et al 2011/Whole_Animal_Expressed_genes/Whole_Animal_Expressed_Genes.xlsx")
View(Stages)

Stages <- as.list(Stages)
Stages <- lapply(Stages, function(x) x[!is.na(x)])
str(Stages)

### Loop to create ggplot data for each Stage
# Container for TissuePlot data
AllStages <- vector("list",length = length(Stages))
# Looping TissuePlot through list of Stage expressed genes
for (i in 1:length(AllStages)) {
  tp <- TissuePlot(Stages[[i]])
  tp2 <- tp + ggtitle(as.character(names(Stages[i])))
  AllStages[[i]] <- tp2
  
}
# Test plot
AllStages[[1]]
### Exporting image files for all tissues ###
for (i in 1:length(AllStages)) {
  ggsave(paste0("Stage",i,".tiff"), plot = AllStages[[i]], width = 4, height = 4)
}

### Plotting fold change for stage expressed genes
# Container for FoldChangeTissuePlot data
AllStagesFC <- vector("list",length = length(Stages))
# Looping FoldChangeTissuePlot through list of tissue expressed genes
for (i in 1:length(AllStagesFC)) {
  AllStagesFC[[i]] <- FoldChangeTissuePlot(Stages[[i]], title = names(Stages[i]))
}

# Test plot first tissue
AllStagesFC[[6]]
### Exporting image files for all tissues ###
for (i in 1:length(AllStagesFC)) {
  ggsave(paste0("Stage",i,".FC.tiff"), plot = AllStagesFC[[i]], width = 3.5, height = 3.5)
}




#################################################################################################################

#################################################################################################################

################################## N. P a r i s i i ( B a k o w s k i _ 2 0 1 4 ) ###############################

#################################################################################################################



# N. Parisii infection DEGs from Bakowski 2014
# Import lists of microsporidia DEGs from Spreadsheet
microsporidia <- read_excel("~/UCSB/Rothman Lab/Data/Published Data/Bakowski et al 2014 Microsporidia/N.Parisii_DEGs.xlsx", sheet = "WBGene")
View(microsporidia)

microsporidia <- as.list(microsporidia)
microsporidia <- lapply(microsporidia, function(x) x[!is.na(x)])
str(microsporidia)

### Loop to create ggplot data for each Stage
# Container for TissuePlot data
Allmicrosporidia <- vector("list",length = length(microsporidia))
# Looping TissuePlot through list of microsporidia DEGs
for (i in 1:length(Allmicrosporidia)) {
  tp <- TissuePlot(microsporidia[[i]])
  tp2 <- tp + ggtitle(as.character(names(microsporidia[i])))
  Allmicrosporidia[[i]] <- tp2
}
# Test plot
Allmicrosporidia[[7]]

### Exporting image files for all tissues ###
for (i in 1:length(Allmicrosporidia)) {
  ggsave(paste0("Microsporidia",i,".tiff"), plot = Allmicrosporidia[[i]], width = 3.5, height = 3.5)
}

### Plotting fold change for stage expressed genes
# Container for FoldChangeTissuePlot data
AllmicrosporidiaFC <- vector("list",length = length(microsporidia))
# Looping FoldChangeTissuePlot through list of tissue expressed genes
for (i in 1:length(AllmicrosporidiaFC)) {
  AllmicrosporidiaFC[[i]] <- FoldChangeTissuePlot(microsporidia[[i]], title = names(microsporidia[i]))
}

# Test plot first tissue
AllmicrosporidiaFC[[7]]
### Exporting image files for all tissues ###
for (i in 1:length(AllmicrosporidiaFC)) {
  ggsave(paste0("Microsporidia",i,".FC.tiff"), plot = AllmicrosporidiaFC[[i]], width = 3.5, height = 3.5)
}




############################################################################################################################

#############################################################################################################################

#############################################################################################################################

#############################################################################################################################

# 20181205 DEG overview plots

# Import lists of DEGs up and downreg and not DEG overlaps from Spreadsheet
DEG.overview <- read_excel("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Rstudio/20180105 RNAseq-TC/DEG.overview.xlsx")
View(DEG.overview)

DEG.overview <- as.list(DEG.overview)
DEG.overview <- lapply(DEG.overview, function(x) x[!is.na(x)])
str(DEG.overview)

### Loop to create ggplot data for each Stage
# Container for TissuePlot data
DEG.tissueplot.data <- vector("list",length = length(DEG.overview))
# Looping TissuePlot through list of microsporidia DEGs
for (i in 1:length(DEG.tissueplot.data)) {
  tp <- TissuePlot(DEG.overview[[i]])
  tp2 <- tp + ggtitle(as.character(names(DEG.overview[i])))
  DEG.tissueplot.data[[i]] <- tp2
}
# Test plot
DEG.tissueplot.data[[1]]

###############

### loop for all tissue lists

# Container for TissuePlot data
DEG.FCtissueplot.data <- vector("list",length = length(DEG.overview))
# Looping FoldChangeTissuePlot through list of tissue expressed genes
for (i in 1:length(DEG.FCtissueplot.data)) {
  
  DEG.FCtissueplot.data[[i]] <- FoldChangeTissuePlot2(DEG.overview[[i]],
                                               title = as.character(names(DEG.overview[i])))#paste0(names(Enrich.Tissues.list[i]), " (", length(Enrich.Tissues.list[[i]]), " genes)"))
  
}
# Test plot first tissue
DEG.FCtissueplot.data[[1]]
### Exporting image files for all tissues ###
for (i in 1:length(DEG.FCtissueplot.data)) {
  ggsave(paste0(as.character(names(DEG.overview[i])),".FC.noMean.tiff"), plot = DEG.FCtissueplot.data[[i]], width = 3.5, height = 3.5)
}

##########################################################################################################

## separate plots for all genes in a list, and just the means

GeneFCSummary <- function(x){
  FCdata <- FC.long[FC.long$Gene %in% x,]
  FCdata.avg <- summaryBy(FC ~ HrsPHS, data = FCdata, FUN = c(length,mean,sd))
  return(FCdata.avg)
}

# test of new function: GeneFCSummary
tissue1.summary <- GeneFCSummary(Enrich.Tissues.list[[1]])
tissue1.summary

DEG.overview.Summaries <- vector("list", length = length(DEG.overview))
for (i in 1:length(DEG.overview.Summaries)) {
  DEG.overview.Summaries[[i]] <- GeneFCSummary(DEG.overview[[i]])
}
names(DEG.overview.Summaries) <- names(DEG.overview)
DEG.overview.Summaries[[2]]$FC.length[1]
DEG.overview.Summaries


## plot just the mean ##
DEG.FCmean.plots <- vector("list", length = length(DEG.overview.Summaries))
for (i in 1:length(DEG.overview.Summaries)) {
  p <- ggplot(DEG.overview.Summaries[[i]], aes(x = HrsPHS, y = FC.mean, group = NA)) +
    # theme_classic() +
    # coord_cartesian(ylim = c(-1.5,1.5)) +
    geom_line(size = 2.5) +
#    geom_errorbar(aes(ymin = FC.mean-FC.sd, ymax = FC.mean+FC.sd), width = 0.2) +
    labs(x="Hours PHS", 
         y="log2 Fold Change", 
         title=as.character(names(DEG.overview.Summaries[i]))) +
    scale_x_discrete(labels = c("3","6","12","30"))
  
  DEG.FCmean.plots[[i]] <- p
}
# test ggplot loop output
DEG.FCmean.plots[[4]]

# save plots
for (i in 1:length(DEG.FCmean.plots)) {
  ggsave(paste0(as.character(names(DEG.overview.Summaries[i])),".FC.means.SD.tiff"), plot = DEG.FCmean.plots[[i]], width = 3.5, height = 3.5)
  
}

# save plots (no error bars)
for (i in 1:length(DEG.FCmean.plots)) {
  ggsave(paste0(as.character(names(DEG.overview.Summaries[i])),".FC.means.tiff"), plot = DEG.FCmean.plots[[i]], width = 3.5, height = 3.5)
  
}




