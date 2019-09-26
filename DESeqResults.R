# Wald tests for the log2 fold changes at individual time points can be investigated using the test argument to results:

resultsNames(ddsTC)

# [1] "Intercept"               "Strain_JR3642_vs_CL2070" "HrsPHS_3_vs_0"           "HrsPHS_6_vs_0"           "HrsPHS_12_vs_0"         
# [6] "HrsPHS_30_vs_0"          "StrainJR3642.HrsPHS3"    "StrainJR3642.HrsPHS6"    "StrainJR3642.HrsPHS12"   "StrainJR3642.HrsPHS30"  


res3hr <- results(ddsTC, name="StrainJR3642.HrsPHS3", test="Wald")
res3hr[which.min(resTC$padj),]
head(res3hrPHS)

res6hr <- results(ddsTC, name = "StrainJR3642.HrsPHS6", test = "Wald")
head(res6hr)

res12hr <- results(ddsTC, name = "StrainJR3642.HrsPHS12", test = "Wald")
head(res12hr)

res30hr <- results(ddsTC, name = "StrainJR3642.HrsPHS30", test = "Wald")
head(res30hr)

# Exporting results to CSV files
# write.csv(as.data.frame(resOrdered), 
#           file="condition_treated_results.csv")

write.csv(as.data.frame(res3hr), file = "20180209_StrainJR3642.HrsPHS3.csv")
write.csv(as.data.frame(res6hr), file = "20180209_StrainJR3642.HrsPHS6.csv")
write.csv(as.data.frame(res12hr), file = "20180209_StrainJR3642.HrsPHS12.csv")
write.csv(as.data.frame(res30hr), file = "20180209_StrainJR3642.HrsPHS30.csv")

############################################
### My own plotting of genes of interest ###
############################################

# fbxc-29
WBGene00019923 <- plotCounts(ddsTC, "WBGene00019923", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00019923, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00019923 fbxc-29")

# elt-7	WBGene00015981
WBGene00015981 <- plotCounts(ddsTC, "WBGene00015981", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00015981, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00015981 elt-7") #+ scale_y_log10()


# end-3	WBGene00001311
WBGene00001311 <- plotCounts(ddsTC, "WBGene00001311", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00001311, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00001311 end-3")


# elt-2	WBGene00001250
WBGene00001250 <- plotCounts(ddsTC, "WBGene00001250", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00001250, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00001250 elt-2")


# ifb-2	WBGene00002054
WBGene00002054 <- plotCounts(ddsTC, "WBGene00002054", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00002054, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00002054 ifb-2")


# pha-4	WBGene00004013
WBGene00004013 <- plotCounts(ddsTC, "WBGene00004013", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00004013, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00004013 pha-4")


# smo-1	WBGene00004888
WBGene00004888 <- plotCounts(ddsTC, "WBGene00004888", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00004888, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00004888 smo-1")


# epc-1	WBGene00007030
WBGene00007030 <- plotCounts(ddsTC, "WBGene00007030", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00007030, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00007030 epc-1")


# pyp-1	WBGene00008149
WBGene00008149 <- plotCounts(ddsTC, "WBGene00008149", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00008149, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00008149 pyp-1")


# cep-1	WBGene00000467
WBGene00000467 <- plotCounts(ddsTC, "WBGene00000467", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00000467, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00000467 cep-1")

# egl-1	WBGene00001170
WBGene00001170 <- plotCounts(ddsTC, "WBGene00001170", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00001170, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00001170 egl-1")


# ced-3	WBGene00000417
WBGene00000417 <- plotCounts(ddsTC, "WBGene00000417", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00000417, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00000417 ced-3")


# ced-4	WBGene00000418
WBGene00000418 <- plotCounts(ddsTC, "WBGene00000418", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00000418, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00000418 ced-4")


# ced-9	WBGene00000423
WBGene00000423 <- plotCounts(ddsTC, "WBGene00000423", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00000423, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00000423 ced-9")

# WBGene00001641 which.max(resTC$padj)
WBGene00001641 <- plotCounts(ddsTC, "WBGene00001641", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00001641, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00001641")

# htz-1 WBGene00019947
WBGene00019947 <- plotCounts(ddsTC, "WBGene00019947", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00019947, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00019947 htz-1")

# sod-1 WBGene00004930
WBGene00004930 <- plotCounts(ddsTC, "WBGene00004930", intgroup = c("Strain","HrsPHS"), returnData = T)

ggplot(WBGene00004930, aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00004930 sod-1")

# fln-1 WBGene00022048
plotCounts(ddsTC, "WBGene00022048", intgroup = c("Strain","HrsPHS"), returnData = T) %>%
ggplot(aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00022048 fln-1")

resTC["WBGene00022048",]

# ddn-1 WBGene00015227 - downstream of daf nineteen
plotCounts(ddsTC, "WBGene00015227", intgroup = c("Strain","HrsPHS"), returnData = T) %>%
  ggplot(aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00015227 ddn-1")

# WBGene00015224
plotCounts(ddsTC, "WBGene00015224", intgroup = c("Strain","HrsPHS"), returnData = T) %>%
  ggplot(aes(x = as.numeric(HrsPHS), y = count, color = Strain, group = Strain)) + 
  geom_point() + geom_smooth(se = F, method = "auto") + ggtitle("WBGene00015224")


