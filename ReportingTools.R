source("https://bioconductor.org/biocLite.R")
biocLite("ReportingTools")
library(ReportingTools)


# Genome wide annotation for Worm, primarily based on mapping using Entrez Gene identifiers.
# Citation (from within R, enter citation("org.Ce.eg.db")):
# Carlson M (2017). org.Ce.eg.db: Genome wide annotation for Worm. R package version 3.5.0.
biocLite("org.Ce.eg.db")
library("org.Ce.eg.db")

org.Ce.egCHRLENGTHS
org.Ce.eg.db
org.Ce.egGENENAME
rownames(ddsTC)[1:100]
columns(org.Ce.eg.db)
keys(org.Ce.eg.db)
select(org.Ce.eg.db, rownames(ddsTC)[1:20], columns = "GENENAME")


# Now the results can be written to a report using the DESeqDataSet object.
des2Report <- HTMLReport(shortName = '20180209_RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./reports")


publish(ddsTC, des2Report, pvalueCutoff=0.05,
          factor = colData(ddsTC),
          reportDir="./reports")
# annotation.db="org.Ce.eg.db", 

finish(des2Report)
