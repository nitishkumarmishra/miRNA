library(pacman)
p_load(DESeq2, ggrepel, VennDiagram, clusterProfiler, pathview, org.Hs.eg.db, edgeR, limma)

adj.method = "BH"; adj.pval = 0.01; raw.pval = 0.01; logFC = 1; hmTopUpN = 15; hmTopDownN = 15
x.cut=2; y.cut=0.01

miRNAExp <- read.csv("CHOL.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
miRNAExp <- miRNAExp[,seq(1,135,3)]
cancerID <- grep("-01A", colnames(miRNAExp))
normalID <- grep("-11A", colnames(miRNAExp))
cnts <- cbind(miRNAExp[,cancerID], miRNAExp[,normalID])
cnts <- data.matrix(cnts)
cnts <- cnts[apply(cnts,1,function(x) sum(x==0))<ncol(cnts)*0.8,]
cnts <- cnts[apply(cnts,2,function(x) sum(x==0))<nrow(cnts)*0.8,]

keep <- rowSums(cpm(cnts)>1) >= ncol(cnts)*0.1 #### Select only genes which have have CPM > 1 for >=50% samples
cnts <- cnts[keep,]
############DESeq2 ###################
cond <- factor(c(rep(2, length(cancerID)), rep(1, length(normalID))))
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
dds <- DESeq(dds)
resmiRNA <- results(dds, lfcThreshold = 1, pAdjustMethod = "BH")
resmiRNAOrdered <- resmiRNA[order(resmiRNA$padj),]
resmiRNASig <- subset(resmiRNAOrdered, padj <= 0.05)
tmp <- as.data.frame(resmiRNASig)
tmp$foldchange <- 2^tmp$log2FoldChange
write.csv(as.data.frame(tmp),file="miRNA_results.csv")

################# edgeR ###############
factors <- factor(c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(cnts, group=factors)
tmm <- calcNormFactors(cnts.dgelist, method = "TMM")
tmmScaleFactors <- cnts.dgelist$samples$lib.size * tmm$samples$norm.factors # equation from the edgeR documentation for estimating normalized absolute expression from their scaling factors
tmmExp <- round(t(t(tmm$counts)/tmmScaleFactors) * mean(tmmScaleFactors)) # equation from the edgeR documentation for estimating normalized absolute expression from their scaling factors
#tmm <- estimateCommonDisp(tmm)
#tmm <- estimateTagwiseDisp(tmm)
tmm <- estimateDisp(tmm) ## estimateDisp command do both CommonDisp and TagwiseDis
tmm.DE <- exactTest(tmm)
#tmm.DE.top <- topTags(tmm.DE, n=nrow(tmm.DE), sort.by="logFC", p.value = 0.01)
tmm.DE.top <- topTags(tmm.DE, n=nrow(tmm.DE), sort.by="logFC")
topTags.DEG <- tmm.DE.top[tmm.DE.top$table$FDR < 0.01 & abs(tmm.DE.top$table$logFC) >=1,]
diffExp <- topTags.DEG$table
rownames(diffExp) <- rownames(diffExp)
write.csv(diffExp, file = "DiffExp-miRNA_EdgeR_TMM.csv")

############ limma+voom ##############
design <- model.matrix(~0 + factor(c(rep(2, length(cancerID)), rep(1, length(normalID)))))
colnames(design) <- c("Normal", "Tumor")
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
cnf <- calcNormFactors(cnts, method = "TMM") ## TMM normalization 
v <- voom(cnts, design, plot = TRUE, lib.size=colSums(cnts) * cnf) ##Transform count data to log2-counts per million (logCPM)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aradeger <- topTable(fit2, adjust.method = "BH", number = length(fit2), sort.by='logFC')
aradeger <- data.frame(aradeger[aradeger$adj.P.Val < 0.01 & aradeger$P.Value < 0.01 & abs(aradeger$logFC) >=1,])
diffExp.voom <- aradeger
rownames(diffExp.voom) <- rownames(diffExp.voom)
write.csv(diffExp.voom, file = "DiffExp-miRNA_VoomLimma_TMM.csv")
save.image(file = "miRNA.RData")

### Volcano plot ####
volcanoplot(fit2, names = fit2$genes$ID, xlab = "Log Fold Change", ylab = "Log Odds", pch = 16, cex = 0.35)
### Heatmap plot ####
if (nrow(aradeger) > 2) {
  aradeger <- aradeger[order(aradeger[, 1], decreasing = TRUE), ]
  if (nrow(aradeger) >= (hmTopDownN + hmTopUpN)) {
    if (hmTopUpN > 0) {
      topgenes <- rownames(aradeger)[1:hmTopUpN]
    } else {
      topgenes <- NULL
    }
    if (hmTopDownN > 0) {
      bottomgenes <- rownames(aradeger)[(nrow(aradeger) - (hmTopDownN - 1)):nrow(aradeger)]
    } else {
      bottomgenes <- NULL
    }
    bluered <- colorRampPalette(c("blue", "white", "red"))(256)
    v <- v[c(topgenes, bottomgenes), ]
    v <- apply(v, 2, as.numeric)
    rownames(v) <- c(topgenes, bottomgenes)
    try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", Colv = NA), silent = FALSE)
  } else {
    bluered <- colorRampPalette(c("blue", "white", "red"))(256)
    v <- v[rownames(aradeger), ]
    v <- apply(v, 2, as.numeric)
    rownames(v) <- rownames(aradeger)
    try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", Colv = NA), silent = FALSE)
  }
}

#### Volcano plot with name and color ########
genes <- tmm.DE.top$table
genes$Gene <- rownames(genes)
Significant <- ifelse(genes$logFC >= x.cut & genes$FDR < y.cut, "Upregulated", ifelse(genes$logFC <= -x.cut & genes$FDR < y.cut, "Downregulated", "Insignificant"))
genes$Significant <- factor(Significant, levels = c("Downregulated", "Upregulated", "Insignificant"), labels =c("Downregulated", "Upregulated", "Insignificant")) ## Reoder, Down, Up and then Insignificant
x.cut=2;y.cut=0.01
p <- ggplot(genes, aes(x = logFC, y = -log10(FDR)),size = 4) +
  geom_point(aes(color = Significant), size = 1) +
  scale_color_manual(values = c("red", "green4", "gray20")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  theme(legend.position="right")+
  ggtitle("Volcano plot Tumor Vs. Normal")+
  theme(plot.title = element_text(size = 14, face = "bold",hjust=0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ 
  theme(legend.title=element_text(colour="black", size = 12, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) +
  geom_text_repel(
    ### used p value 0.00001 log2FC >=abs(6) for name in plot
    data = subset(genes, (FDR < 1e-2 & abs(logFC) >=2)),
    aes(label = Gene),size = 3, 
    show.legend = FALSE, fontface = "bold", color = "black", 
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  )
ggsave(p, filename = "VolcanoPlotmiRNA.pdf", width = 12, height = 8, dpi = 800)