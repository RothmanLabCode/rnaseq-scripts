# #EBseqHMM
# #20180108
# ###
# #Gene and transcript level analysis of expression changes over time points.
# #Algorithm does not appear to take control data into account - I will likely need to split my data into hs-elt-7 and hs-gfp samples,
# #then perform parallel analyses, and compare expression paths...?
# 
# #First install required packages from bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("EBSeqHMM")
# browseVignettes("EBSeqHMM")
# 
# #load in required libraries
library(EBSeq)
library(EBSeqHMM)
# 
# Required inputs: Data in a 'Gene' by 'Sample' matrix *obtained via Salmon implemented in Galaxy*
# Need to split 32 samples into two sets of 16
library(readr)
GenesTPM_5TP <- read_csv("~/UCSB/Sequencing Analysis/2018 mRNAseq-TC/Galaxy Output/Genes TPM Table_Samples 1-32.csv", col_names = T)
head(GenesTPM_5TP)
#c("c",rep("n",32))
GenesTPM_5TP$`1-S1_T1_R1`
str(GenesTPM_5TP)

e7samples <- c(1,2,3,4,8,9,10,14,15,16,20,21,22,30,31,32,33)

gfpSamples <- c(1,5,6,7,11,12,13,17,18,19,23,24,25,26,27,28,29)

head(GenesTPM_5TP[,e7samples],n = 20)
elt7geneTPM_5tp <- GenesTPM_5TP[,e7samples]
gfpGeneTPM_5tp <- GenesTPM_5TP[,gfpSamples]
head(elt7geneTPM_5tp)
# head(gfpGeneTPM)
# #nice
# 
# #Required intput: Conditions - factor of time points. I now have 5 time points, 3,4,3,3 replicates each = 16 samples
# CondVecE7 <- rep(paste("t",1:4,sep = ""), each = 3)
# CondVecE7
conditionvector <- c(rep("0hr",3),rep("3hr",4),rep("6hr",3),rep("12hr",3),rep("24hr",3))
conditionvector
# conditions <- factor(CondVecE7, levels = c("t1","t2","t3","t4"))
conditions <- factor(conditionvector, levels = c("0hr","3hr","6hr","12hr","24hr"))
str(conditions)

##########library size factor - adjust for differences in sequencing depth between samples
# #Sizes <- MedianNorm(elt7geneTPM)
# #Need to make the first column into row names
##also need to resort columns arranged in chronological order
row.names(elt7geneTPM_5tp)
colnames(elt7geneTPM_5tp)

elt7df_5tp <- data.frame(cbind(elt7geneTPM_5tp[2:4],elt7geneTPM_5tp[14:17],elt7geneTPM_5tp[5:13]), row.names = elt7geneTPM_5tp$WBGene)

gfpdf_5tp <- data.frame(cbind(gfpGeneTPM_5tp[2:4],gfpGeneTPM_5tp[14:17],gfpGeneTPM_5tp[5:13]), row.names = gfpGeneTPM_5tp$WBGene)


Sizes_elt7 <- MedianNorm(elt7df_5tp)
Sizes_elt7
plot(Sizes_elt7)
 
Sizes_GFP <- MedianNorm(gfpdf_5tp)
Sizes_GFP
plot(Sizes_GFP)

### Generate a normalized gene expression matrix
# NormELT7matrix <- GetNormalizedMat(Data = elt7geneTPM12, Sizes = Sizes)
NormE7_5tp <- GetNormalizedMat(Data = elt7df_5tp, Sizes = Sizes_elt7)
# NormGFPmatrix <- GetNormalizedMat(gfpGeneTPM12, SizesGFP)
NormGFP_5tp <- GetNormalizedMat(gfpdf_5tp, Sizes_GFP)


##############################
#Visualize expression path of a particular gene of interest
# #elt-7 = WBGene00015981 
par(mfrow = c(3,2))
PlotExp(NormE7_5tp, conditions, Name = "WBGene00015981")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00015981")


##########################################
##Expression of potential gut factors relevant to Ethan###
#odd-1: WBGene00003845
PlotExp(NormE7_5tp, conditions, Name = "WBGene00003845")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00003845")

#odd-2: WBGene00003846
PlotExp(NormE7_5tp, conditions, Name = "WBGene00003846")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00003846")
#########################################################
# WBGene00006843 (unc-119)
par(mfrow = c(2,1))
PlotExp(NormE7_5tp, conditions, Name = "WBGene00006843")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00006843")


PlotExp(NormE7_5tp, conditions, Name = "WBGene00010819")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00010819")

# END-3 WBGene00001311
PlotExp(NormE7_5tp, conditions, Name = "WBGene00001311")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00001311")


# ###Running EBSeqHMM on gene expression estimates
# ##Can also be repeated with various number of iterations, other variable parameters, etc...
# ?EBSeqHMMTest
# 
# EBSeqHMMgeneoutELT7 <- EBSeqHMMTest(Data=as.matrix(elt7geneTPM12), Conditions = conditions, sizeFactors = Sizes)
EBSeqOutELT7_5tp <- EBSeqHMMTest(Data = as.matrix(elt7df_5tp), Conditions = conditions, sizeFactors = Sizes_elt7)
str(EBSeqOutELT7_5tp)
# 
# EBSeqHMMgeneoutGFP <- EBSeqHMMTest(Data = as.matrix(gfpGeneTPM12), Conditions = conditions, sizeFactors = SizesGFP)
EBSeqOutGFP_5tp <- EBSeqHMMTest(Data = as.matrix(gfpdf_5tp), Conditions = conditions, sizeFactors = Sizes_GFP)
 
#################################################################
###Detection of DE genes and inference of gene’s most likely path
###using more stringent FDR of 0.01
elt7geneDEcalls_5tp <- GetDECalls(EBSeqOutELT7_5tp, FDR = 0.01)
head(elt7geneDEcalls_5tp)

PlotExp(NormE7_5tp, conditions, Name = "WBGene00006956")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00006956")

gfpGeneDEcalls_5tp <- GetDECalls(EBSeqOutGFP_5tp, FDR = 0.01)
head(gfpGeneDEcalls_5tp)


######################################################
##### Clustering DE genes into expression paths ######
# elt7GeneConfCalls <- GetConfidentCalls(EBSeqHMMOut = EBSeqHMMgeneoutELT7, FDR = .05, cutoff = .5, OnlyDynamic = T)
elt7GeneConfCalls_5tp <- GetConfidentCalls(EBSeqHMMOut = EBSeqOutELT7_5tp, FDR = .01, cutoff = .75, OnlyDynamic = F)
elt7GeneConfCalls_5tp$NumEach
sum(elt7GeneConfCalls_5tp$NumEach)
#1,037 genes!
elt7GeneConfCalls_5tp$EachPathNames$`Up-Up-Up-Up`

PlotExp(NormE7_5tp, conditions, Name = "WBGene00011063")
PlotExp(NormGFP_5tp, conditions, Name = "WBGene00011063")

### GFP Control
gfpGeneConfCalls_5tp <- GetConfidentCalls(EBSeqHMMOut = EBSeqOutGFP_5tp, FDR = .01, cutoff = .75, OnlyDynamic = F)
gfpGeneConfCalls_5tp$EachPath
gfpGeneConfCalls_5tp$NumEach
sum(gfpGeneConfCalls_5tp$NumEach)
#528 genes




# ####Diagnostic Plots
par(mfrow = c(2,2))
QQP(EBOut = EBSeqOutELT7_5tp, GeneLevel = TRUE)
QQP(EBSeqOutGFP_5tp, GeneLevel = T)
#looks pretty good for both
DenNHist(EBOut = EBSeqOutELT7_5tp, GeneLevel = T)  #Does not look great
DenNHist(EBSeqOutGFP_5tp, GeneLevel = T) #Somewhat better, still very jagged data
# ####???What do these diagnostics show??? Is my data/analysis okay?! Generates different plots each time...??



############################
### Data Visualization #####
#Top 12 genes
top12elt7genenames <- rownames(elt7GeneConfCalls_5tp$Overall[1:12,])
print(top12elt7genenames)
par(mfrow = c(3,4))
for (i in 1:12) {
  PlotExp(NormE7_5tp, conditions, Name = top12elt7genenames[i])
}


#Top 20 elt-7 genes
top20elt7genenames <- rownames(elt7GeneConfCalls_5tp$Overall[1:20,])
print(top20elt7genenames)
par(mfrow = c(4,5))
for (i in 1:20) {
  PlotExp(NormE7_5tp, conditions, Name = top20elt7genenames[i])
}

#Top 20 GFP genes
top20gfpGenenames <- rownames(gfpGeneConfCalls_5tp$Overall[1:20,])
print(top20gfpGenenames)
par(mfrow = c(4,5))
for (i in 1:20) {
  PlotExp(NormGFP_5tp, conditions, Name = top20gfpGenenames[i])
}




sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DESeq2_1.18.1              SummarizedExperiment_1.8.1 DelayedArray_0.4.1         matrixStats_0.52.2         Biobase_2.38.0            
# [6] GenomicRanges_1.30.0       GenomeInfoDb_1.14.0        IRanges_2.12.0             S4Vectors_0.16.0           BiocGenerics_0.24.0       
# [11] readr_1.1.1                tximport_1.6.0            
# 
# loaded via a namespace (and not attached):
#   [1] locfit_1.5-9.1         Rcpp_0.12.14           lattice_0.20-35        digest_0.6.14          R6_2.2.2               plyr_1.8.4            
# [7] backports_1.1.2        acepack_1.4.1          RSQLite_2.0            ggplot2_2.2.1          pillar_1.1.0           zlibbioc_1.24.0       
# [13] rlang_0.1.6            lazyeval_0.2.1         rstudioapi_0.7         data.table_1.10.4-3    annotate_1.56.1        blob_1.1.0            
# [19] rpart_4.1-11           Matrix_1.2-12          checkmate_1.8.5        splines_3.4.3          BiocParallel_1.12.0    geneplotter_1.56.0    
# [25] stringr_1.2.0          foreign_0.8-69         htmlwidgets_0.9        RCurl_1.95-4.10        bit_1.1-12             munsell_0.4.3         
# [31] compiler_3.4.3         pkgconfig_2.0.1        base64enc_0.1-3        htmltools_0.3.6        nnet_7.3-12            tibble_1.4.1          
# [37] gridExtra_2.3          htmlTable_1.11.1       GenomeInfoDbData_1.0.0 Hmisc_4.1-1            XML_3.98-1.9           bitops_1.0-6          
# [43] grid_3.4.3             xtable_1.8-2           gtable_0.2.0           DBI_0.7                magrittr_1.5           scales_0.5.0          
# [49] stringi_1.1.6          XVector_0.18.0         genefilter_1.60.0      latticeExtra_0.6-28    Formula_1.2-2          RColorBrewer_1.1-2    
# [55] tools_3.4.3            bit64_0.9-7            hms_0.4.0              survival_2.41-3        AnnotationDbi_1.40.0   colorspace_1.3-2      
# [61] cluster_2.0.6          memoise_1.1.0          knitr_1.18            
#