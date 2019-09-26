#EBseqHMM
#20180108
###
#Gene and transcript level analysis of expression changes over time points.
#Algorithm does not appear to take control data into account - I will likely need to split my data into hs-elt-7 and hs-gfp samples,
#then perform parallel analyses, and compare expression paths...?

#First install required packages from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("EBSeqHMM")
browseVignettes("EBSeqHMM")

#load in required libraries
library(EBSeq)
library(EBSeqHMM)

#Required inputs: Data in a 'Gene' by 'Sample' matrix *obtained via Salmon implemented in Galaxy*
#Need to split 24 samples into two sets of 12

GeneTPM$X1.S1_T1_R1
str(GeneTPM)  

e7samples <- c(1,2,3,4,8,9,10,14,15,16,20,21,22)

gfpSamples <- c(1,5,6,7,11,12,13,17,18,19,23,24,25)

GeneTPM[,e7samples]
elt7geneTPM <- GeneTPM[,e7samples]
gfpGeneTPM <- GeneTPM[,gfpSamples]
head(elt7geneTPM)
head(gfpGeneTPM)
#nice

#Required intput: Conditions - factor of time points. I have 4 time points, 3 replicates each = 12 samples
CondVecE7 <- rep(paste("t",1:4,sep = ""), each = 3)
CondVecE7

conditions <- factor(CondVecE7, levels = c("t1","t2","t3","t4"))
str(conditions)

#library size factor - adjust for differences in sequencing depth between samples
#Sizes <- MedianNorm(elt7geneTPM)
#Need to make the first column into row names
row.names(elt7geneTPM)
elt7geneTPM12 <- data.frame(elt7geneTPM[,2:13], row.names = elt7geneTPM$WBGene)
Sizes <- MedianNorm(elt7geneTPM12)
Sizes
plot(Sizes)

gfpGeneTPM12 <- data.frame(gfpGeneTPM[,2:13], row.names = gfpGeneTPM$WBGene)
SizesGFP <- MedianNorm(gfpGeneTPM12)
SizesGFP
plot(SizesGFP)

#Generate a normalized gene expression matrix
NormELT7matrix <- GetNormalizedMat(Data = elt7geneTPM12, Sizes = Sizes)

NormGFPmatrix <- GetNormalizedMat(gfpGeneTPM12, SizesGFP)

#Visualize expression path of a particular gene of interest
#elt-7 = WBGene00015981

PlotExp(NormELT7matrix, conditions, Name = "WBGene00015981")
PlotExp(NormGFPmatrix, conditions, Name = "WBGene00015981")

###Running EBSeqHMM on gene expression estimates
##Can also be repeated with various number of iterations, other variable parameters, etc...
?EBSeqHMMTest

EBSeqHMMgeneoutELT7 <- EBSeqHMMTest(Data=as.matrix(elt7geneTPM12), Conditions = conditions, sizeFactors = Sizes)
str(EBSeqHMMgeneoutELT7)

EBSeqHMMgeneoutGFP <- EBSeqHMMTest(Data = as.matrix(gfpGeneTPM12), Conditions = conditions, sizeFactors = SizesGFP)

###Detection of DE genes and inference of gene’s most likely path
elt7geneDEcalls <- GetDECalls(EBSeqHMMgeneoutELT7, FDR = 0.05)
head(elt7geneDEcalls)
#WBGene00003057 "Down-Down-Up"   "1"
PlotExp(NormELT7matrix, conditions, Name = "WBGene00003057")
PlotExp(NormGFPmatrix, conditions, Name = "WBGene00003057")

gfpGeneDEcalls <- GetDECalls(EBSeqHMMgeneoutGFP, FDR = 0.05)
head(gfpGeneDEcalls)
#WBGene00014125 top DE call "Down-Up-Down"
PlotExp(NormGFPmatrix, conditions, Name = "WBGene00014125")
PlotExp(NormELT7matrix, conditions, Name = "WBGene00014125")


####Clustering DE genes into expression paths
elt7GeneConfCalls <- GetConfidentCalls(EBSeqHMMOut = EBSeqHMMgeneoutELT7, FDR = .05, cutoff = .5, OnlyDynamic = T)
elt7GeneConfCalls$EachPath[1:4]
elt7GeneConfCalls$NumEach
#elt7GeneConfCalls$NumEach
#Up-Up-Up     Down-Up-Up     Up-Down-Up   Down-Down-Up     Up-Up-Down   Down-Up-Down   Up-Down-Down Down-Down-Down 
#54            615            652            348             93            869           1032            162 


gfpGeneConfCalls <- GetConfidentCalls(EBSeqHMMOut = EBSeqHMMgeneoutGFP, FDR = .05, cutoff = .5, OnlyDynamic = T)
gfpGeneConfCalls$EachPath
gfpGeneConfCalls$NumEach
# gfpGeneConfCalls$NumEach
#Up-Up-Up     Down-Up-Up     Up-Down-Up   Down-Down-Up     Up-Up-Down   Down-Up-Down   Up-Down-Down Down-Down-Down 
#801            384            188            101            113            323           1201            275 
sum(gfpGeneConfCalls$NumEach) #3386


elt7GeneConfCalls$Overall
elt7GeneConfCalls$EachPath
PlotExp(NormELT7matrix, conditions, Name = "WBGene00000215")
elt7GeneConfCalls$NumEach
sum(elt7GeneConfCalls$NumEach)  #3825
elt7GeneConfCalls$EachPathNames

####Diagnostic Plots
par(mfrow = c(2,2))
QQP(EBOut = EBSeqHMMgeneoutELT7, GeneLevel = TRUE)
QQP(EBSeqHMMgeneoutGFP, GeneLevel = T)
#looks pretty good for both
DenNHist(EBOut = EBSeqHMMgeneoutELT7, GeneLevel = T)  #Does not look great
DenNHist(EBSeqHMMgeneoutGFP, GeneLevel = T) #Somewhat better, still very jagged data
####???What do these diagnostics show??? Is my data/analysis okay?!

###Data Visualization
#Top 12 genes
top12elt7genenames <- rownames(elt7GeneConfCalls$Overall[1:12,])
print(top12elt7genenames)
par(mfrow = c(3,4))
for (i in 1:12) {
  PlotExp(NormELT7matrix, conditions, Name = top12elt7genenames[i])
}

#############################################################################
########### DE calling with higher cutoffs and all paths ####################
#############################################################################






