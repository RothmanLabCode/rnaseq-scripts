##### Project to cross reference "complete" list of C. elegans ChIPseq peaks with nearby TSS of DEGs from 2018 RNAseq-TC #####

## Need to import data:
### ChIPseq peak file (large); RNAseq DEGs; gene TSSs


chipseq <- read.delim(file = "SupplementalTable6Bwormsites.txt", header = F)
chipcolnames <- c("Chr","PkStart","PkEnd","ExptName","V5","V6","SignalVal","V8","qval","PkSummittoStart","PkWidth","ce","ExptName2","Strain","TF","Stage","WA","ChrSummit")
colnames(chipseq) <- chipcolnames
head(chipseq)
# most critical data from "chipseq" = V1("Chr"), V15("TF"), V18("ChrSummit"), V4("ExptName"), V16("Stage")
chipseqChromPeaks <- chipseq[,c("Chr","TF","ChrSummit","ExptName","Stage")]
chipseqChromPeaks
head(chipseqChromPeaks)[,3]
# need to remove the "stn " from "ChrSummit" column

library(dplyr)
library(tidyr)
# split "ChrSummit" column into two parts using 'separate()', and take away the "STN" part
chipseqChromPeaks <- chipseqChromPeaks %>% separate(ChrSummit, into = c("STN", "Summit"), sep = " ", remove = T)
# head(chipseqChromPeaks)[,-3]
#str(chipseqChromPeaks)
chipseqChromPeaks <- chipseqChromPeaks[,-3]

# what are all of the chromosome names and transcription factor names in the chipseq data table?
ChromosomeNames <- levels(chipseqChromPeaks[,1])
TranscriptionFactorNames <- levels(chipseqChromPeaks[,2])
length(TranscriptionFactorNames)
# there are 217 TFs with ChIPseq data

#remove "chr..." from Chomosome designations
chipseqChromPeaks$Chr <- substr(chipseqChromPeaks$Chr, 4, 10)
chipseqChromPeaks$Chr <- as.factor(chipseqChromPeaks$Chr)
levels(chipseqChromPeaks$Chr)
head(chipseqChromPeaks)
str(chipseqChromPeaks)
summary(chipseqChromPeaks)
#split up whole genome file into each chromosome
chr1.chip <- chipseqChromPeaks[chipseqChromPeaks$Chr=="I",]
chr2.chip <- chipseqChromPeaks[chipseqChromPeaks$Chr=="II",]
chr3.chip <- chipseqChromPeaks[chipseqChromPeaks$Chr=="III",]
chr4.chip <- chipseqChromPeaks[chipseqChromPeaks$Chr=="IV",]
chr5.chip <- chipseqChromPeaks[chipseqChromPeaks$Chr=="V",]
chrX.chip <- chipseqChromPeaks[chipseqChromPeaks$Chr=="X",]
chrMt.chip <- chipseqChromPeaks[chipseqChromPeaks$Chr=="MtDNA",]


### gene TSSs: can I extract gene locations from the canonical_geneset.gtf.gz file? ###
geneset <- read.delim(file = "c_elegans.PRJNA13758.WS263.canonical_geneset.gtf.gz", header = F)
# only want "gene" (or CDS?) from feature col V3
genesetcolnames <- c("seqname","source","feature","start","end","score","strand","frame","group")
colnames(geneset) <- genesetcolnames
# remove first row comment from geneset table
geneset <- geneset[-1,]
# subset only desired rows
# "gene" feature includes non-coding RNAs and pseudogenes, but "CDS" has a new line for every exon of every gene...?
levels(geneset[,3])
head(geneset[geneset$feature=="gene",],30)
head(geneset[geneset$feature=="CDS",],30)
levels(geneset$seqname)
geneset.feat.gene <- geneset[geneset$feature=="gene",]
# geneset.feat.CDS <- geneset[geneset$feature=="CDS",]

## create table of gene_biotype = protein_coding and paried with gene_id and "start", etc.
geneset.feat.gene <- separate(data = geneset.feat.gene, col = group, into = c("id","gene_id","source2", "gene_source","biotype", "gene_biotype"), sep = " ")

geneset.protein_coding <- geneset.feat.gene[geneset.feat.gene$gene_biotype=="protein_coding;",]
geneset.protein_coding <- geneset.protein_coding[,1:10]

geneset.protein_coding <- separate(geneset.protein_coding, gene_id, into = c("gene_id","semi"), sep = ";")
geneset.protein_coding <- geneset.protein_coding[,c(1:8,10)]
summary(geneset.protein_coding)

#split up whole genome gene set into each chomosome
chr1.pcg <- geneset.protein_coding[geneset.protein_coding$seqname=="I",]
chr2.pcg <- geneset.protein_coding[geneset.protein_coding$seqname=="II",]
chr3.pcg <- geneset.protein_coding[geneset.protein_coding$seqname=="III",]
chr4.pcg <- geneset.protein_coding[geneset.protein_coding$seqname=="IV",]
chr5.pcg <- geneset.protein_coding[geneset.protein_coding$seqname=="V",]
chrX.pcg <- geneset.protein_coding[geneset.protein_coding$seqname=="X",]
chrMt.pcg <- geneset.protein_coding[geneset.protein_coding$seqname=="MtDNA",]



##### Cross referencing protein coding genes with called ChIPseq peaks #####
head(chr1.chip)
head(chr1.pcg)
# First: which of the rows of "Summit" from '.chip' are within a promoter region for each gene in '.pcg'
#chr1.pcg$strand <- as.character(chr1.pcg$strand)
# chr1.pcg$strand[1]=="+"


# if(chr1.pcg$strand[2]=="+") print("plus")

# set the distances from the TSS to look for ChIPseq summits
promoter.upstream <- 1000
promoter.downstream <- 250
# create columns for each end of the defined promoter interval
# !critical! where the TSS is for a gene depends on the strand!
chr1.pcg$promoter.up <- if_else(chr1.pcg$strand=="+", true = chr1.pcg$start-promoter.upstream, false = chr1.pcg$end+promoter.upstream)
chr1.pcg$promoter.down <- if_else(chr1.pcg$strand=="+", true = chr1.pcg$start+promoter.downstream, false = chr1.pcg$end-promoter.downstream)
# comparing a promoter interval to ChIPseq summits would be easier if ordered left-to-right
chr1.pcg$prom.left <- if_else(chr1.pcg$promoter.up<chr1.pcg$promoter.down, 
        true = chr1.pcg$promoter.up, 
        false = chr1.pcg$promoter.down)
chr1.pcg$prom.right <- if_else(chr1.pcg$promoter.up<chr1.pcg$promoter.down,
                               true = chr1.pcg$promoter.down,
                               false = chr1.pcg$promoter.up)
stopifnot(chr1.pcg$prom.left<chr1.pcg$prom.right)

# 10000 > chr1.pcg$promoter.down[1]

# between(as.numeric(chr1.chip$Summit), chr1.pcg$prom.left[1], chr1.pcg$prom.right[1])

# row1test <- chr1.chip[between(as.numeric(chr1.chip$Summit), chr1.pcg$prom.left[10], chr1.pcg$prom.right[10]),]
# row1test$Gene <- chr1.pcg$gene_id[10]
# row1test
# row2test <- chr1.chip[between(as.numeric(chr1.chip$Summit), chr1.pcg$prom.left[200], chr1.pcg$prom.right[200]),]
# row2test$Gene <- chr1.pcg$gene_id[200]
# row2test
# row3test <- rbind(row1test,row2test)
# row3test
# typeof(chr1.chip$Summit)
# typeof(chr1.pcg$promoter.up[1])
# chr1.pcg$promoter.up[1]
# length(chr1.pcg$seqname)


chr1.Summit2Gene <- data.frame()

for (i in 1:length(chr1.pcg$seqname)) {
  new.df <- chr1.chip[between(as.numeric(chr1.chip$Summit), chr1.pcg$prom.left[i], chr1.pcg$prom.right[i]),]
  if(length(new.df[,1])>1){
  new.df$Gene <- chr1.pcg$gene_id[i]
  chr1.Summit2Gene <- rbind(chr1.Summit2Gene,new.df)
  }
}


# new.df <- chr1.chip[between(as.numeric(chr1.chip$Summit), chr1.pcg$prom.left[6], chr1.pcg$prom.right[6]),]
# if(length(new.df[,1])>1){
# new.df$Gene <- chr1.pcg$gene_id[6]
# chr1.Summit2Gene <- rbind(chr1.Summit2Gene,new.df)}


#### defining a funciton Summit2Gene to replicate mappings from pcg and chip data frames ####

Summit2Gene <- function(pcg, chip, upstream, downstream){
  pcg$promoter.up <- if_else(pcg$strand=="+", true = pcg$start-upstream, false = pcg$end+upstream)
  pcg$promoter.down <- if_else(pcg$strand=="+", true = pcg$start+downstream, false = pcg$end-downstream)
  # comparing a promoter interval to ChIPseq summits would be easier if ordered left-to-right
  pcg$prom.left <- if_else(pcg$promoter.up<pcg$promoter.down, 
                                true = pcg$promoter.up, 
                                false = pcg$promoter.down)
  pcg$prom.right <- if_else(pcg$promoter.up<pcg$promoter.down,
                                 true = pcg$promoter.down,
                                 false = pcg$promoter.up)
  stopifnot(pcg$prom.left<pcg$prom.right)
  #create container for queried data
  S2G <- data.frame()
  # loop for all protein coding gene promoters
  for (i in 1:length(pcg$seqname)) {
    new.df <- chip[between(as.numeric(chip$Summit), pcg$prom.left[i], pcg$prom.right[i]),]
    if(length(new.df[,1])>1){
      new.df$Gene <- pcg$gene_id[i]
      S2G <- rbind(S2G,new.df)
    }
  }
  return(S2G)
}

### Summit2Gene on MtDNA ###
chrMt.Summit2Gene <- Summit2Gene(pcg = chrMt.pcg, chip = chrMt.chip, upstream = 1000, downstream = 250)
head(chrMt.Summit2Gene)
chrMt.Summit2Gene

### remap chr1 using Summit2Gene function ###
chr1.Summit2Gene.fun <- Summit2Gene(pcg = chr1.pcg,
                                    chip = chr1.chip,
                                    upstream = 1000,
                                    downstream = 250)
summary(chr1.Summit2Gene.fun)

### chr2
chr2.Summit2Gene <- Summit2Gene(pcg = chr2.pcg,
                                chip = chr2.chip,
                                upstream = 1000,
                                downstream = 250)
### chr3
chr3.Summit2Gene <- Summit2Gene(pcg = chr3.pcg,
                                chip = chr3.chip,
                                upstream = 1000,
                                downstream = 250)
head(chr3.Summit2Gene, 50)
### chr4
chr4.Summit2Gene <- Summit2Gene(pcg = chr4.pcg,
                                chip = chr4.chip,
                                upstream = 1000,
                                downstream = 250)
head(chr4.Summit2Gene, 50)
### chr5
chr5.Summit2Gene <- Summit2Gene(pcg = chr5.pcg,
                                chip = chr5.chip,
                                upstream = 1000,
                                downstream = 250)
head(chr5.Summit2Gene, 100)
### chrX
chrX.Summit2Gene <- Summit2Gene(pcg = chrX.pcg,
                                chip = chrX.chip,
                                upstream = 1000,
                                downstream = 250)
head(chrX.Summit2Gene, 100)


##### combining and mining all TF::gene data #####
genome.Summit2Gene <- rbind(chr1.Summit2Gene,chr2.Summit2Gene,chr3.Summit2Gene,chr4.Summit2Gene,
                            chr5.Summit2Gene,chrX.Summit2Gene,chrMt.Summit2Gene)
head(genome.Summit2Gene, 100)
str(genome.Summit2Gene)
summary(genome.Summit2Gene)
summary(genome.Summit2Gene$Stage)

## subselect for Stages of intrest
# L3, L4, YA, yAd
stages <- c("L3","L4","YA","yAd")
# stages <- "L4"
# genome.Summit2Gene.stages <- genome.Summit2Gene[genome.Summit2Gene$Stage==stages,] #doesn't work correctly??
# genome.Summit2Gene.L3 <- genome.Summit2Gene[genome.Summit2Gene$Stage=="L3",]


genome.Summit2Gene.stages <- data.frame()
for (i in stages) {
  new.df <- genome.Summit2Gene[genome.Summit2Gene$Stage==i,]
  genome.Summit2Gene.stages <- rbind(genome.Summit2Gene.stages,new.df)
}

# 49438+70931+8353+44988
# = 173710

# Export .csv for TF-gene mappings
write.csv(genome.Summit2Gene.stages, file = "genome.summit2Gene.L3-YA.csv")

############################################################################################
#### cross reference RNAseq DEGs ####
## importing DEGs ##

RNAseq.DEGs <- read.csv(file = "20190124_FoldChange.AllHrsPHS.csv", header = TRUE)
head(RNAseq.DEGs)
dimnames(RNAseq.DEGs)[[2]][1] <- "WBGene"
colnames(RNAseq.DEGs[1])
# colnames(RNAseq.DEGs[1]) <- "WBGene"
nrow(RNAseq.DEGs)
str(RNAseq.DEGs)
# order(RNAseq.DEGs$ï..WBGene)
RNAseq.DEGs <- RNAseq.DEGs[order(RNAseq.DEGs$WBGene),]


# subsetting DEGs
nrow(RNAseq.DEGs[RNAseq.DEGs$padj.30hr<0.01 & RNAseq.DEGs$padj.12hr<0.01 & RNAseq.DEGs$log2FC.12hr>0 & RNAseq.DEGs$log2FC.30hr>0,])
RNAseq.DEGs[RNAseq.DEGs$padj.30hr<0.01 & RNAseq.DEGs$padj.12hr<0.01 & RNAseq.DEGs$log2FC.12hr>0 & RNAseq.DEGs$log2FC.30hr>0,]

plot(RNAseq.DEGs$log2FC.3hr[RNAseq.DEGs$baseMean>90000],RNAseq.DEGs$baseMean[RNAseq.DEGs$baseMean>90000])

RNAseq.DEGs$ï..WBGene[RNAseq.DEGs$baseMean>90000 & RNAseq.DEGs$log2FC.6hr< 0]

### What is the "background" frequency of binding sites for each TF in the ChIP-seq data subset? ###
summary(genome.Summit2Gene.stages)
str(genome.Summit2Gene.stages)

library("ggplot2")
#polar label angles: 90-360/113
polarlabelangels <- 90-360/113*0:112
polarlabelangels[113/4]
genome.Summit2Gene.stages$TF
length(unique(genome.Summit2Gene.stages$TF))
# genome.Summit2Gene.stages_TF freq_polar
ggplot(data = genome.Summit2Gene.stages, aes(x = TF)) + geom_bar() + coord_polar() + 
  theme_light() +
  theme(axis.text.x = element_text(angle = polarlabelangels))
# coordcart
ggplot(data = genome.Summit2Gene.stages, aes(x = TF)) + 
  geom_bar() + 
  theme_bw() +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,8500))
  # scale_x_discrete(limits = rev(levels(genome.Summit2Gene.stages$TF))) #reverses labels, but also includes TFs with 0 count...

levels(genome.Summit2Gene.stages$TF)[5]
length(levels(genome.Summit2Gene.stages$TF)) # 217
head(genome.Summit2Gene.stages)
TF.background.counts <- count(genome.Summit2Gene.stages, vars = genome.Summit2Gene.stages$TF)
TF.background.counts[1:20,]

write.csv(x = TF.background.counts, file = "TF.background.counts.L3-YA.csv")

# TFtest <- as.character(TF.background.counts$vars[12])
# TFtest
# "CEH-9"
# genome.Summit2Gene.stages$Gene
# genome.Summit2Gene.stages$Gene == TFtest
### all the genes which have a CEH-9 ChIP-seq peak (at L3-Ad stage) ###
# genome.Summit2Gene.stages$Gene[genome.Summit2Gene.stages$TF == TFtest]

TF.vector <- as.vector(TF.background.counts$vars)
TF.vector


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


### For genes upregulated at 6.hrsPHS (padj < 0.01), how many TF ChIP peaks are present? ###
Up.6 <- RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.6hr > 0 & RNAseq.DEGs$padj.6hr < 0.01]
Up.6 <- as.vector(Up.6)# a few thousand gene names
Up.6[1]
# copy chip-seq mapping and rename to something easier
ChIP.map <- genome.Summit2Gene.stages
head(ChIP.map)
# ChIP.map[ChIP.map$TF==TF.vector[1] & ChIP.map$Gene==Up.6,]
ChIP.map[ChIP.map$Gene=="WBGene00000001",]
ChIP.map$Gene=="WBGene00022277"
ChIP.map[ChIP.map$Gene=="WBGene00022277",]
ChIP.map[ChIP.map$Gene==Up.6[c(2,5)],]
TF.vector[2]
ChIP.map[ChIP.map$Gene==Up.6[5] & ChIP.map$TF==TF.vector[2],]

### cross comparisons for a single gene set (Up.6) and single TF (ATF-7) ###
TF.vector
ATF7.map <- ChIP.map[ChIP.map$TF=="ATF-7",]
ATF7.map
Up.6.ATF7 <- ATF7.map[ATF7.map$Gene %in% Up.6,]
Up.6.ATF7
length(Up.6.ATF7$Gene) # 259 ATF-7 binding sites
length(unique(Up.6.ATF7$Gene)) # 248 genes
length(Up.6.ATF7$Gene)/length(ATF7.map$Gene) # 0.259 (26% of background sites in Up.6)
length(unique(Up.6.ATF7$Gene))/length(Up.6) # 0.0602 (6% of genes upregulated 6hrsPHS have ATF-7 ChIP peaks...)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# ATF7.map
# TF.vector
TF.list <- as.list(TF.vector)
TF.list

# Break up the genome wide list of all TF:Gene mappings into separate data frames for each TF, all contained within one List object #
for (tf in 1:length(TF.list)) {
  TF.list[[tf]] <- ChIP.map[ChIP.map$TF==TF.vector[tf],]
}

# head(TF.list)
## rename lists by TF name
names(TF.list) <- TF.vector
str(TF.list)
TF.list$`FKH-3`


### Need the background number of genes which are bound by each TF ###

# TF.background.genes

length(TF.list[[10]]$Gene) # 2523 = number of TF binding sites (from L3-YA)
length(unique(TF.list[[10]]$Gene)) # 1849 = number of GENES with TF binding sites (from L3-YA)

# starting data.frame container for summarized TF::gene binding data
TF.background.genes <- data.frame(TF.vector)
# looping through the list of TF binding sites and counting the number of unique genes for each
for (i in 1:length(TF.background.genes$TF.vector)) {
  TF.background.genes$num.genes[i] <- length(unique(TF.list[[i]]$Gene))
}

head(TF.background.genes)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Create a data frame for a given DEG set (Up.6), where each TF has a row, and columns for:
## (1) number of TF binding sites "num.bs"
## (2) number of genes w/ >= 1 ChIP peak "num.genes"
## (3) percent of the background sites present in the DEG set "pc.bkgnd"
## (4) percent of the DEG set genes which have >= 1 ChIP peak "pc.genes"

Up.6.df <- data.frame(TF = TF.vector)

for (i in 1:length(Up.6.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Up.6,]
  Up.6.df$num.bs[i] <- length(overlap$Gene)
  Up.6.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.6.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.6.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.6)
}

Up.6.df
### PLOT: % of DEGs with TF binding peaks
ggplot(Up.6.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Upregulated 6 hrsPHS (", length(Up.6), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.3))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Up.6.df[Up.6.df$pc.genes>0.15,], mapping = aes(x = TF, y = pc.genes*100,
                                                      size = Up.6.df$num.bs[Up.6.df$pc.genes>0.15],
                                                      colour = Up.6.df$pc.bs[Up.6.df$pc.genes>0.15]*100)) + 
  geom_point() +
  ggtitle(paste0("Upregulated 6 hrsPHS (", length(Up.6), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(15,30)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")
## This figure potentially illustrates the `correlation` between TFBS and DEGs, the `magnitude` of binding sites for DEGs
## (with larger dots = more ChIP peaks), and the `specificity` of a transcription factor's binding profile for the given set of DEGs
## (brighter blue means a greater proportion of the overall, genome-wide binding profile for a given TF is present within the DEG set).

######################################################################################################################################
######################################################################################################################################

##### Continuing analysis for major DEG sets, Up3,6,12,30 and Down3,6,12,30 #####
colnames(RNAseq.DEGs)
Up.3 <- as.vector(RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.3hr > 0 & RNAseq.DEGs$padj.3hr < 0.01])
Up.3[1:20] # a few thousand genes
Up.12 <- as.vector(RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.12hr > 0 & RNAseq.DEGs$padj.12hr < 0.01])
Up.30 <- as.vector(RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.30hr > 0 & RNAseq.DEGs$padj.30hr < 0.01])

Down.3 <- as.vector(RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.3hr < 0 & RNAseq.DEGs$padj.3hr < 0.01])
Down.6 <- as.vector(RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.6hr < 0 & RNAseq.DEGs$padj.6hr < 0.01])
Down.12 <- as.vector(RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.12hr < 0 & RNAseq.DEGs$padj.12hr < 0.01])
Down.30 <- as.vector(RNAseq.DEGs$WBGene[RNAseq.DEGs$log2FC.30hr < 0 & RNAseq.DEGs$padj.30hr < 0.01])

## Create a data frame for a given DEG set (Up.6), where each TF has a row, and columns for:
## (1) number of TF binding sites "num.bs"
## (2) number of genes w/ >= 1 ChIP peak "num.genes"
## (3) percent of the background sites present in the DEG set "pc.bkgnd"
## (4) percent of the DEG set genes which have >= 1 ChIP peak "pc.genes"

Up.3.df <- data.frame(TF = TF.vector)

for (i in 1:length(Up.3.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Up.3,]
  Up.3.df$num.bs[i] <- length(overlap$Gene)
  Up.3.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.3.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.3.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.3)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Up.12.df <- data.frame(TF = TF.vector)

for (i in 1:length(Up.12.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Up.12,]
  Up.12.df$num.bs[i] <- length(overlap$Gene)
  Up.12.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.12.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.12.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.12)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Up.30.df <- data.frame(TF = TF.vector)

for (i in 1:length(Up.30.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Up.30,]
  Up.30.df$num.bs[i] <- length(overlap$Gene)
  Up.30.df$num.genes[i] <- length(unique(overlap$Gene))
  Up.30.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Up.30.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Up.30)
}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.3.df <- data.frame(TF = TF.vector)

for (i in 1:length(Down.3.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Down.3,]
  Down.3.df$num.bs[i] <- length(overlap$Gene)
  Down.3.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.3.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.3.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.3)
}
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.6.df <- data.frame(TF = TF.vector)

for (i in 1:length(Down.6.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Down.6,]
  Down.6.df$num.bs[i] <- length(overlap$Gene)
  Down.6.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.6.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.6.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.6)
}
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.12.df <- data.frame(TF = TF.vector)

for (i in 1:length(Down.12.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Down.12,]
  Down.12.df$num.bs[i] <- length(overlap$Gene)
  Down.12.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.12.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.12.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.12)
}
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
Down.30.df <- data.frame(TF = TF.vector)

for (i in 1:length(Down.30.df$TF)) {
  tmp <- TF.list[[i]]
  overlap <- tmp[tmp$Gene %in% Down.30,]
  Down.30.df$num.bs[i] <- length(overlap$Gene)
  Down.30.df$num.genes[i] <- length(unique(overlap$Gene))
  Down.30.df$pc.bkgnd[i] <- length(overlap$Gene)/length(tmp$Gene)
  Down.30.df$pc.genes[i] <- length(unique(overlap$Gene))/length(Down.30)
}
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

##############################################################################################################################

############################## D O T  P L O T S ##############################################################################

##############################################################################################################################

###### Bar graphs & Dotplots for additional DEG sets ######

#### Up.3

### PLOT: % of DEGs with TF binding peaks
ggplot(Up.3.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Upregulated 3 hrsPHS (", length(Up.3), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.3))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Up.3.df[Up.3.df$pc.genes>0.15,], mapping = aes(x = TF, y = pc.genes*100,
                                                      size = Up.3.df$num.bs[Up.3.df$pc.genes>0.15],
                                                      colour = Up.3.df$pc.bkgnd[Up.3.df$pc.genes>0.15]*100)) + 
  geom_point() +
  ggtitle(paste0("Upregulated 3 hrsPHS (", length(Up.3), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(15,30)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Up.3.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#### Up.12

### PLOT: % of DEGs with TF binding peaks
ggplot(Up.12.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Upregulated 12 hrsPHS (", length(Up.12), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.3))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Up.12.df[Up.12.df$pc.genes>0.15,], mapping = aes(x = TF, y = pc.genes*100,
                                                      size = Up.12.df$num.bs[Up.12.df$pc.genes>0.15],
                                                      colour = Up.12.df$pc.bkgnd[Up.12.df$pc.genes>0.15]*100)) + 
  geom_point() +
  ggtitle(paste0("Upregulated 12 hrsPHS (", length(Up.12), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(15,30)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Up.12.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#### Up.30

### PLOT: % of DEGs with TF binding peaks
ggplot(Up.30.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Upregulated 30 hrsPHS (", length(Up.30), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.2))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Up.30.df[Up.30.df$pc.genes>0.05,], mapping = aes(x = TF, y = pc.genes*100,
                                                        size = Up.30.df$num.bs[Up.30.df$pc.genes>0.05],
                                                        colour = Up.30.df$pc.bkgnd[Up.30.df$pc.genes>0.05]*100)) + 
  geom_point() +
  ggtitle(paste0("Upregulated 30 hrsPHS (", length(Up.30), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(5,20)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Up.30.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#### Down.3

### PLOT: % of DEGs with TF binding peaks
ggplot(Down.3.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Downregulated 3 hrsPHS (", length(Down.3), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.3))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Down.3.df[Down.3.df$pc.genes>0.15,], mapping = aes(x = TF, y = pc.genes*100,
                                                        size = Down.3.df$num.bs[Down.3.df$pc.genes>0.15],
                                                        colour = Down.3.df$pc.bkgnd[Down.3.df$pc.genes>0.15]*100)) + 
  geom_point() +
  ggtitle(paste0("Downregulated 3 hrsPHS (", length(Down.3), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(15,30)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Down.3.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#### Down.6

### PLOT: % of DEGs with TF binding peaks
ggplot(Down.6.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Downregulated 6 hrsPHS (", length(Down.6), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.3))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Down.6.df[Down.6.df$pc.genes>0.15,], mapping = aes(x = TF, y = pc.genes*100,
                                                          size = Down.6.df$num.bs[Down.6.df$pc.genes>0.15],
                                                          colour = Down.6.df$pc.bkgnd[Down.6.df$pc.genes>0.15]*100)) + 
  geom_point() +
  ggtitle(paste0("Downregulated 6 hrsPHS (", length(Down.6), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(15,30)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Down.6.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#### Down.12

### PLOT: % of DEGs with TF binding peaks
ggplot(Down.12.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Downregulated 12 hrsPHS (", length(Down.12), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.4))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Down.12.df[Down.12.df$pc.genes>0.2,], mapping = aes(x = TF, y = pc.genes*100,
                                                          size = Down.12.df$num.bs[Down.12.df$pc.genes>0.2],
                                                          colour = Down.12.df$pc.bkgnd[Down.12.df$pc.genes>0.2]*100)) + 
  geom_point() +
  ggtitle(paste0("Downregulated 12 hrsPHS (", length(Down.12), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(20,40)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Down.12.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#### Down.30

### PLOT: % of DEGs with TF binding peaks
ggplot(Down.30.df, mapping = aes(x = TF, y = pc.genes)) + 
  geom_col() +
  ggtitle(paste0("Downregulated 30 hrsPHS (", length(Down.30), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.45))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Down.30.df[Down.30.df$pc.genes>0.2,], mapping = aes(x = TF, y = pc.genes*100,
                                                           size = Down.30.df$num.bs[Down.30.df$pc.genes>0.2],
                                                           colour = Down.30.df$pc.bkgnd[Down.30.df$pc.genes>0.2]*100)) + 
  geom_point() +
  ggtitle(paste0("Downregulated 30 hrsPHS (", length(Down.30), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(20,45)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "ChIP peaks", colour = "% of background")

ggsave(filename = "Down.30.TF dotplot.tiff", device = "tiff", width = 14, height = 10, units = "cm")


#######################################################################################################################################

########################################################### H E A T  M A P ############################################################

#######################################################################################################################################

##### Compiling results and plotting all DEG lists and TFs in a Heatmap #####

## Need to create a matrix of data for heatmap. Can also use `as.matrix()` on a data frame.
library(gplots)
help(heatmap.2)

name.heatmap.cols <- c("Up 3","Up 6","Up 12","Up 30","Down 3","Down 6","Down 12","Down 30")
pc.genes.df <- data.frame(TF.vector, Up.3.df$pc.genes, Up.6.df$pc.genes, Up.12.df$pc.genes, Up.30.df$pc.genes,
                          Down.3.df$pc.genes, Down.6.df$pc.genes, Down.12.df$pc.genes, Down.30.df$pc.genes)
View(pc.genes.df)

pc.genes.matrix <- data.frame(Up.3.df$pc.genes, Up.6.df$pc.genes, Up.12.df$pc.genes, Up.30.df$pc.genes,
                              Down.3.df$pc.genes, Down.6.df$pc.genes, Down.12.df$pc.genes, Down.30.df$pc.genes)
rownames(pc.genes.matrix) <- TF.vector
colnames(pc.genes.matrix) <- name.heatmap.cols
pc.genes.matrix <- as.matrix(pc.genes.matrix)
is.matrix(pc.genes.matrix)
pc.genes.matrix[1:10,1:6]

heatmap.2(pc.genes.matrix)

### Full heatmap for all 113 TFs
heatmap.2(pc.genes.matrix, trace = "none", dendrogram = "row", Colv = F, col = bluered(20), 
          scale = "none", key = T, keysize = 1.1, density.info = "none")

## trimming TFs with the least % genes overall (down to 58)
# rowSums(pc.genes.matrix)
# plot(rowSums(pc.genes.matrix))
z <- rowSums(pc.genes.matrix) > 0.3
pc.genes.matrix.trim <- pc.genes.matrix[z,]

## trimmed heatmap
heatmap.2(pc.genes.matrix.trim, trace = "none", dendrogram = "row", Colv = F, col = bluered(20), 
          scale = "none", key = T, keysize = 1.1, density.info = "none")

### matrix for percent background

pc.bkgnd.matrix <- data.frame(Up.3.df$pc.bkgnd, Up.6.df$pc.bkgnd, Up.12.df$pc.bkgnd, Up.30.df$pc.bkgnd,
                              Down.3.df$pc.bkgnd, Down.6.df$pc.bkgnd, Down.12.df$pc.bkgnd, Down.30.df$pc.bkgnd)
rownames(pc.bkgnd.matrix) <- TF.vector
colnames(pc.bkgnd.matrix) <- name.heatmap.cols
pc.bkgnd.matrix[1:15,1:6]
pc.bkgnd.matrix <- as.matrix(pc.bkgnd.matrix)
is.matrix(pc.bkgnd.matrix)

## heatmap for Percent of Background TF binding sites [NOT genes] found within the DEG sets ##
heatmap.2(pc.bkgnd.matrix, trace = "none", dendrogram = "row", Colv = F, col = bluered(20), density.info = "none")

### overlay comments on cells with % background values
heatmap.2(pc.genes.matrix, trace = "none", dendrogram = "row", Colv = F, col = bluered(20), 
          scale = "none", key = T, keysize = 1.1, density.info = "none",
          cellnote = round(pc.bkgnd.matrix, 2), notecol = "black")

pc.bkgnd.matrix.trim <- pc.bkgnd.matrix[z,]

### trimmed heatmap overlay comments on cells with % background values
heatmap.2(pc.genes.matrix.trim, trace = "none", dendrogram = "row", Colv = F, col = bluered(20), 
          scale = "none", key = T, keysize = 1.1, density.info = "none",
          cellnote = round(pc.bkgnd.matrix.trim, 2), notecol = "black")

#########################################################################################################################################

################################### S T A T S ###########################################################################################

#########################################################################################################################################

##### Statistical analysis of correlation significance using Fisher's Exact tests #####

### test using data for Up.3 and ELT-2 ChIP

Up.3.df[Up.3.df$TF=="ELT-2",]
Up.3.df$num.genes[Up.3.df$TF=="ELT-2"]
# (a) num.genes = 1084
length(Up.3) # (v) 4231
length(RNAseq.DEGs$WBGene) # (z) 20094 - all genes tested

## !!!!!!!!!!!! computational error !!!!! background counts is the number of TF binding sites, NOT the number of Genes !!!!!
TF.background.counts$n[TF.background.counts$vars=="ELT-2"] # (x) 5063
## !!!!!!!!!!!!

TF.background.genes$num.genes[TF.background.genes$TF.vector=="ELT-2"] # (x) 4544

4544-1084 # (b) 3460
4231-1084 # (c) 3147
20094-1084-3460-3147 # (d) 12403
# 2x2 contingency table
ct <- matrix(c(1084,3147,3460,12403),nrow = 2)
ct
ct.ft <- fisher.test(ct)
attributes(ct.ft)
ct.ft
barplot(ct, beside = F, legend.text = T)
dimnames(ct) <- list(TF.bound = c("Yes","No"),DEG = c("Upreg","not DE"))

## test2: Up.12 and NHR-80
a <- Up.12.df$num.genes[Up.12.df$TF=="NHR-80"]
v <- length(Up.12)
z <- length(RNAseq.DEGs$WBGene)
# x <- TF.background.counts$n[TF.background.counts$vars=="NHR-80"] # WRONG BACKGROUND!
x <- TF.background.genes$num.genes[TF.background.genes$TF.vector=="NHR-80"]
b <- x-a
c <- v-a
d <- z-a-b-c
ct <- matrix(c(a,c,b,d),nrow=2)
ct.ft <- fisher.test(ct)
ct.ft


### for loop to compute p-value for all TFs in a DEG set ###

# Up.3.df$pval <- 1
# Up.3.df
# as.vector(Up.3.df$TF[4]) == TF.background.counts$vars[4]

# # # # # # # # # # # # # # # #
# # # Up.3 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Up.3.df$num.genes[TF]
  v <- length(Up.3)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.3.df$pval[TF] <- ct.ft$p.value
}

Up.3.df

# # # # # # # # # # # # # # # #
# # # Up.6 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Up.6.df$num.genes[TF]
  v <- length(Up.6)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.6.df$pval[TF] <- ct.ft$p.value
}

Up.6.df

# # # # # # # # # # # # # # # #
# # # Up.12 fisher tests  # # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Up.12.df$num.genes[TF]
  v <- length(Up.12)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.12.df$pval[TF] <- ct.ft$p.value
}

Up.12.df

# # # # # # # # # # # # # # # #
# # # Up.30 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Up.30.df$num.genes[TF]
  v <- length(Up.30)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Up.30.df$pval[TF] <- ct.ft$p.value
}

Up.30.df

# # # # # # # # # # # # # # # #
# # # Down.3 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Down.3.df$num.genes[TF]
  v <- length(Down.3)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.3.df$pval[TF] <- ct.ft$p.value
}

Down.3.df

# # # # # # # # # # # # # # # #
# # # Down.6 fisher tests # # # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Down.6.df$num.genes[TF]
  v <- length(Down.6)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.6.df$pval[TF] <- ct.ft$p.value
}

Down.6.df

# # # # # # # # # # # # # # # #
# # # Down.12 fisher tests  # # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Down.12.df$num.genes[TF]
  v <- length(Down.12)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.12.df$pval[TF] <- ct.ft$p.value
}

Down.12.df

# # # # # # # # # # # # # # # #
# # # Down.30 fisher tests  # #
# # # # # # # # # # # # # # # #

for (TF in 1:113) {
  a <- Down.30.df$num.genes[TF]
  v <- length(Down.30)
  z <- length(RNAseq.DEGs$WBGene)
  x <- TF.background.genes$num.genes[TF]
  b <- x-a
  c <- v-a
  d <- z-a-b-c
  ct <- matrix(c(a,c,b,d),nrow=2)
  ct.ft <- fisher.test(ct)
  Down.30.df$pval[TF] <- ct.ft$p.value
}

Down.30.df

#### Plotting p-values ####

### PLOT: -log of p-val of significant association between DEGs with TF binding peaks
ggplot(Up.3.df, mapping = aes(x = TF, y = -log(pval))) + 
  geom_col() +
  ggtitle(paste0("Upregulated 3 hrsPHS (", length(Up.3), " genes)"), subtitle = "-log(p-value) of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,30))

### Bubblechart - Plot with dot size = number of ChIP peaks, color = % background
## Subset the dataframe to plot only TFs with peaks in at least 15% of the DEGs
ggplot(Up.3.df[Up.3.df$pc.genes>0.1,], mapping = aes(x = TF, y = pc.genes*100,
                                                     size = Up.3.df$pc.bkgnd[Up.3.df$pc.genes>0.1]*100,
                                                      # size = Up.3.df$num.bs[Up.3.df$pc.genes>0.1],
                                                      colour = -log(Up.3.df$pval[Up.3.df$pc.genes>0.1]))) + 
  geom_point() +
  ggtitle(paste0("Upregulated 3 hrsPHS (", length(Up.3), " genes)"), subtitle = "Percent of DEGs with TF ChIP-seq peaks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(10,30)) +
  scale_size_continuous(range = c(2.5,10)) +
  labs(x = "Transcription Factor", y = "% of genes", size = "% of binding sites", colour = "-log10(p-value)")

ggsave(filename = "Up.3.TF dotplot_pvals.tiff", device = "tiff", width = 14, height = 10, units = "cm")


#############################################################################################

####################### P - V A L U E  H E A T  M A P #######################################

#############################################################################################


## Need to create a matrix of data for heatmap. Can also use `as.matrix()` on a data frame.
library(gplots)
# help(heatmap.2)

name.heatmap.cols <- c("Up 3","Up 6","Up 12","Up 30","Down 3","Down 6","Down 12","Down 30")


pvals.matrix <- data.frame(Up.3.df$pval, Up.6.df$pval, Up.12.df$pval, Up.30.df$pval,
                              Down.3.df$pval, Down.6.df$pval, Down.12.df$pval, Down.30.df$pval)
rownames(pvals.matrix) <- TF.vector
colnames(pvals.matrix) <- name.heatmap.cols
pvals.matrix <- as.matrix(pvals.matrix)
is.matrix(pvals.matrix)
pvals.matrix[1:10,1:6]
heatmap.2(pvals.matrix)

write.csv(pvals.matrix, "TF.DEGs.pvalues.csv")

pvals.matrix.log <- -log(pvals.matrix)

pvals.matrix.log[is.infinite(pvals.matrix.log)] = 666
is.infinite(pvals.matrix.log)
heatmap.2(pvals.matrix.log)


### Full heatmap of p-values for all 113 TFs
heatmap.2(pvals.matrix.log, trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 4)

### Full heatmap of p-values for all 113 TFs | Version 2
heatmap.2(pvals.matrix.log, trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "col", key = T, keysize = 1.1, density.info = "none", colsep = 4)



p <- pvals.matrix.log > 200
pvals.matrix.log.bin <- pvals.matrix.log
pvals.matrix.log.bin[p] = 200

### Full heatmap of p-values for all 113 TFs | Version 3
heatmap.2(pvals.matrix.log.bin, trace = "none", dendrogram = "row", Colv = F, col = bluered(90), 
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 4)


## trimming TFs with larger pvalues (down to 60)
rowSums(pvals.matrix.log)
plot(rowSums(pvals.matrix.log))
z <- rowSums(pvals.matrix.log) > 200
pvals.matrix.log.trim <- pvals.matrix.log[z,]

## trimmed heatmap
heatmap.2(pvals.matrix.log.trim, trace = "none", dendrogram = "row", Colv = F, col = bluered(20), 
          scale = "none", key = T, keysize = 1.1, density.info = "none", colsep = 4)
