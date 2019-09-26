genome.Summit2Gene
genome.Summit2Gene[genome.Summit2Gene$Gene=="WBGene00003624",]

##### What are the genes which are both bound by ELT-2 and Downregulated after hs-ELT-7? #####

### First starting with L3-YA stages of ChIP-seq data ###
TF.list$`ELT-2`$Stage
summary(TF.list$`ELT-2`$Stage) # all L3 (5063 peaks)
TF.list$`ELT-2`$Gene
length(TF.list$`ELT-2`$Gene)
length(unique(TF.list$`ELT-2`$Gene)) # 4544 genes


# vector of genes with >= 1 ELT-2 binding peak in promoter
genes.boundby.ELT2.L3 <- unique(TF.list$`ELT-2`$Gene)
write.csv(x = genes.boundby.ELT2.L3, file = "genes.boundby.ELT2.L3.csv")


# DEGs downreg 6 hrsPHS with ELT-2 TFBS
genes.boundby.ELT2.L3[genes.boundby.ELT2.L3 %in% Down.6]
ELT2.Down.6 <- genes.boundby.ELT2.L3[genes.boundby.ELT2.L3 %in% Down.6]
write.csv(ELT2.Down.6, "ELT2.down.6.csv")
