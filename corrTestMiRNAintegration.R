library(miRLAB)
library(TCGAbiolinks)
#tmp <- read.csv("CHOL-miRLAB-Logit.csv", header = TRUE)
## we use rawcount data for scale, logit are already transferred
tmp <- read.csv("CHOL-miRLAB.csv", header = TRUE)
tmp <- scale(tmp)
x <- t(tmp[,1:122])
y <- t(tmp[,123:4743])
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
ObjNet$miRNA <- gsub("\\.","-", ObjNet$miRNA)
##### Multiple testing ##################################
ObjNet.0.01 <- ObjNet[which(ObjNet$pval <= 0.01&ObjNet$cor < 0),]
ObjNet.0.01$adj.pval <- p.adjust(ObjNet.0.01$pval, method = "BH")
ObjNet.0.01 <- ObjNet.0.01[which(ObjNet.0.01$adj.pval <= 0.01),]
#########################################################
Validation(ObjNet.0.01, "TargetScan.csv")
Validation(ObjNet.0.01, "Groundtruth-all.csv")
###### Pathway analysis of genes controlled by miRNA ####
GeneList <- as.character(unique(ObjNet.0.01$mRNA))
system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Tumor Vs Normal",GeneList1))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = GeneList1,
                        nBar = 10)

