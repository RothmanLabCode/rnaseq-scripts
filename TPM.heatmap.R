###### R Script Objective: Publication Figure - Heatmap for TPMs for all samples and all genes #####
## Import table of TPMs for all samples
TPM.import <- read.csv("C:/Users/spickard/Documents/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Galaxy Output/Genes TPM Table_Samples 1-32.csv")

head(TPM.import)

TPM.import[1:10,1:10]

TPM.dataframe <- data.frame(TPM.import[,-1], row.names = TPM.import[,1])

TPM.dataframe[1:10,1:10]

TPM.matrix <- as.matrix(TPM.dataframe)

TPM.matrix

library("gplots")

heatmap.2(TPM.matrix) #ugly, takes a long time to map

heatmap.2(TPM.matrix[1:1000,],
          trace = "none",
          labRow = NA,
          dendrogram = "column",
          col = redgreen(31),
#          col = rev(blues9),
          scale = "row")

### Need to reorder columns to match with samples and experimental design. Replicates cluster well

my.sample.order <- c(1,2,3,29,30,31,32,7,8,9,13,14,15,19,20,21,4,5,6,25,26,27,28,10,11,12,16,17,18,22,23,24)
TPM.matrix[10001:10005,my.sample.order]

# scaled and red-green colored
heatmap.2(TPM.matrix[1:1000,my.sample.order],
          trace = "none",
          labRow = NA,
          Colv = F,
          dendrogram = "none",
          col = redgreen(31),
          colsep = 16,
          #          col = rev(blues9),
          scale = "row")


# log transformed and black-blue colored
heatmap.2(log2(TPM.matrix[,my.sample.order]+1),
          trace = "none",
          labRow = NA,
          Colv = F,
          dendrogram = "none",
          colsep = 16,
          col = colorpanel(21, "black", "blue","yellow"))
          #col = redgreen(31))
#          col = rev(blues9))
#          scale = "row")


log2(TPM.matrix[1:10,1:10]+0.01)
TPM.matrix[1:10,1:10]
