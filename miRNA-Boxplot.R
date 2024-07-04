library(ggpubr)
load("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/miRNA/Survival_CHOL_miRNA_Median_Sd.RData")
tmp <- as.data.frame(t(miRNAExp))
rownames(tmp) <- substr(rownames(tmp), 1, 16)
tmp.select <- tmp %>% select("hsa-mir-551b", "hsa-mir-22")
colnames(tmp.select) <- c("Mir551b", "Mir22")
tmp.select$Mir551b <- log2(tmp.select$Mir551b+1)
tmp.select$Mir22 <- log(tmp.select$Mir22+1)/2
tmp.select$Sample <- ifelse(substr(rownames(tmp.select), 14,15)=="01", "Tumor", "Normal")



p <- ggboxplot(tmp.select, x = "Sample", y = "Mir551b", color = "Sample",
               palette = c( "red4","blue4"), add = "jitter", ylab = "log2 (mir-551b)")
p + stat_compare_means(method = "t.test", label.y = max(tmp.select$Mir551b)+1)
dev.print(pdf, 'miR551b Boxplot.pdf ', width = 6, height = 6)



p <- ggboxplot(tmp.select, x = "Sample", y = "Mir22", color = "Sample",
               palette = c( "red4","blue4"), add = "jitter", ylab = "log2 (mir-22)")
p + stat_compare_means(method = "t.test", label.y = max(tmp.select$Mir551b)+1)
dev.print(pdf, 'miR22 Boxplot.pdf ', width = 6, height = 6)

