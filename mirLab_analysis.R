library(miRLAB)
library(edgeR)
library(miRComb)

geneExp <- read.csv("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = 1, check.names = FALSE)
geneExp <- geneExp[-c(1),which(geneExp[1,]=="raw_count")]
rownames(geneExp) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(geneExp))
geneExp <- round(data.matrix(geneExp),0)
geneExp <- geneExp[-grep("\\?", rownames(geneExp)),] ## Remove genes which HGNC name is known i.e. ?
rownames(geneExp) <- sapply(strsplit(rownames(geneExp),"\\|"),'[[',1)
colnames(geneExp) <- substr(colnames(geneExp), 1, 16)
geneExp <- geneExp[rowSums(geneExp==0)< ncol(geneExp)*0.25,]
cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))
geneExp <- cbind(geneExp[,cancerID], geneExp[,normalID])
#save(geneExp, file = "CHOL_mRNA_Firehose.Rda")
miRNAExp <- read.csv("CHOL.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
miRNAExp <- miRNAExp[,seq(1,135,3)]
colnames(miRNAExp) <- substr(colnames(miRNAExp), 1, 16)
cancerID <- grep("-01A", colnames(miRNAExp))
normalID <- grep("-11A", colnames(miRNAExp))
cnts <- cbind(miRNAExp[,cancerID], miRNAExp[,normalID])
cnts <- data.matrix(cnts)
cnts <- cnts[apply(cnts,1,function(x) sum(x==0))<ncol(cnts)*0.8,]
cnts <- cnts[apply(cnts,2,function(x) sum(x==0))<nrow(cnts)*0.8,]
miRNAExp <- cnts
#####################################
factors <- factor(c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(geneExp, group=factors)
tmm <- calcNormFactors(cnts.dgelist, method = "TMM")
tmm <- estimateDisp(tmm) ## estimateDisp command do both CommonDisp and TagwiseDis
tmm.DE <- exactTest(tmm)
tmm.DE.top <- topTags(tmm.DE, n=nrow(tmm.DE), sort.by="logFC", p.value = 0.01)
topTags.DEG <- tmm.DE.top[tmm.DE.top$table$FDR < 0.01 & abs(tmm.DE.top$table$logFC) >=1.5,]
diffExp.edgeR <- topTags.DEG$table

factors <- factor(c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(miRNAExp, group=factors)
tmm <- calcNormFactors(cnts.dgelist, method = "TMM")
tmm <- estimateDisp(tmm) ## estimateDisp command do both CommonDisp and TagwiseDis
tmm.DE <- exactTest(tmm)
tmm.DE.top <- topTags(tmm.DE, n=nrow(tmm.DE), sort.by="logFC", p.value = 0.01)
topTags.DEG <- tmm.DE.top[tmm.DE.top$table$FDR < 0.01 & abs(tmm.DE.top$table$logFC) >=1,]
diffExpmiRNA.edgeR <- topTags.DEG$table
#########################################

bindData <- rbind(miRNAExp[rownames(diffExpmiRNA.edgeR),], geneExp[rownames(diffExp.edgeR),])
bindData <- t(bindData)
bindData.log <- log(bindData+1)
write.csv(bindData.log, file = "BindData_log.csv", row.names = FALSE)
#cause=1:99
#effect=100:4368
PearsonCHOL.1 <- Pearson("BindData_log.csv", cause=1:99, effect=100:4368, targetbinding = "TargetScan.csv")

####################################################
tmp <- read.csv("BindData_log.csv", header = TRUE, check.names = FALSE)
x <- t(tmp[,1:99])
y <- t(tmp[,100:4368])
#### Make correlation matrix and pvalue matrix with NA
correlation.matrix <- matrix(NA, nrow = nrow(x), ncol = nrow(y))
colnames(correlation.matrix) <- rownames(y)
rownames(correlation.matrix) <- rownames(x)
pval.matrix <- matrix(NA, nrow = nrow(x), ncol = nrow(y))
colnames(pval.matrix) <- rownames(y)
rownames(pval.matrix) <- rownames(x)
##### Fill correlation and pvalue matrix ##############
for (i in 1:nrow(x)) {
  for (j in 1:nrow(y)) {
    x1 <- x[i, ]
    y1 <- y[j, ]
    cor.t <- cor.test(x1, y1, method = "pearson", alternative = "less")
    correlation.matrix[i, j] <- cor.t$estimate
    pval.matrix[i, j] <- cor.t$p.value
  }
}
corr <- t(correlation.matrix)
pvals <- t(pval.matrix)
ObjNet <- data.frame(miRNA = rep(colnames(corr), each = nrow(corr)), 
                     mRNA = rep(rownames(corr), ncol(corr)), cor = as.vector(corr), 
                     pval = as.vector(pvals))
ObjNet$BH <- p.adjust(ObjNet$pval, method = "BH") 
##### Multiple testing ##################################
ObjNet.0.01 <- ObjNet[which(ObjNet$BH <= 0.01&ObjNet$cor < -0.75),]
#ObjNet.0.01$adj.pval <- p.adjust(ObjNet.0.01$pval, method = "BH")
#ObjNet.0.01 <- ObjNet.0.01[which(ObjNet.0.01$adj.pval <= 0.01),]
rownames(ObjNet.0.01) <- paste(ObjNet.0.01$miRNA, ObjNet.0.01$mRNA, sep = ":")
data("targetScan_v6.2_18")
data("microCosm_v5_18")
intersect(rownames(ObjNet.0.01), rownames(targetScan_v6.2_18))
#########################################################
temp <- read.csv("../Correlation/Merged/results.dms.txt")
miRNA_DNAmeth <- unique(intersect(ObjNet.0.01$mRNA, temp$gene))
eg <- bitr(miRNA_DNAmeth, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
eg <- eg[order(eg$ENTREZID, decreasing = TRUE),]
kk <- enrichKEGG(gene = eg$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
ggo <- groupGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", level=3, readable = TRUE)
barplot(kk, showCategory = 8)## only 8 pathway enriched
###############
eg1 <- bitr(unique(ObjNet.0.01$mRNA), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
kk2 <- enrichKEGG(gene = eg1$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
barplot(kk2, showCategory = 10)


save.image("mirLab_analysis.RData")
