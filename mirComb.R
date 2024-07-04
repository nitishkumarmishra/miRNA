#### Pipeline for miRNA-mRNA integration by using mirComb tool
## In this pipeline, we used log2 of reads and limma for DE analysis
library(miRComb)
library(edgeR)
geneExp <- read.csv("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = 1, check.names = FALSE)
geneExp <- geneExp[-c(1),which(geneExp[1,]=="raw_count")]
rownames(geneExp) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(geneExp))
geneExp <- round(data.matrix(geneExp),0)
geneExp <- geneExp[-grep("\\?", rownames(geneExp)),] ## Remove genes which HGNC name is known i.e. ?
rownames(geneExp) <- sapply(strsplit(rownames(geneExp),"\\|"),'[[',1)
colnames(geneExp) <- substr(colnames(geneExp), 1, 16)

miRNAExp <- read.csv("CHOL.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
miRNAExp <- miRNAExp[,seq(1,135,3)]
colnames(miRNAExp) <- substr(colnames(miRNAExp), 1, 16)
miRNAExp <- data.matrix(miRNAExp)


geneExp <- geneExp[rowSums(geneExp==0)< ncol(geneExp)*0.25,]
miRNAExp <- miRNAExp[rowSums(miRNAExp==0)< ncol(miRNAExp)*0.25,]
#miRNAExp <- miRNAExp[rowMedians(miRNAExp) >= 10, ]
#geneExp <- geneExp[rowMedians(geneExp) > 10, ]
cancerID <- grep("01A", colnames(miRNAExp))
normalID <- grep("11A", colnames(miRNAExp))
miRNAExp <- cbind(miRNAExp[,normalID], miRNAExp[,cancerID])
cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))
geneExp <- cbind(geneExp[,normalID], geneExp[,cancerID])
####################################################################
### log2 transformation of data
#mRNA <- log2(geneExp+1)
#miRNA <- log2(miRNAExp+1)
###############################################################
############# Normalized with voom ##################################
design <- model.matrix(~0 + factor(c(rep(1, length(normalID)), rep(2, length(cancerID)))))
colnames(design) <- c("Normal", "Tumor")
factors <- factor(c(rep("Normal", length(normalID)), rep("Tumor", length(cancerID))))
cnts.dgelist <- DGEList(miRNAExp, group=factors)
cnf <- calcNormFactors(cnts.dgelist, method = "TMM")
v <- voom(cnf, design, plot = FALSE) ##Transform count data to log2-counts per million (logCPM)
miRNA <- v$E

cnts.dgelist <- DGEList(geneExp, group=factors)
cnf <- calcNormFactors(cnts.dgelist, method = "TMM")
v <- voom(cnf, design, plot = FALSE)
mRNA <- v$E


pheno.mRNA <- data.frame(group = c(rep("H", length(normalID)), rep("D", length(cancerID))),
                         DvH = as.numeric(c(rep("0", length(normalID)), rep("1", length(cancerID)))))
pheno.miRNA <- data.frame(group = c(rep("H", length(normalID)), rep("D", length(cancerID))),
                          DvH = as.numeric(c(rep("0", length(normalID)), rep("1", length(cancerID)))))
rownames(pheno.mRNA) <- colnames(mRNA)
rownames(pheno.miRNA) <- colnames(miRNA)


#mRNA <- as.matrix(mRNA)
#miRNA <- as.matrix(miRNA)
#####################################################################
data.obj<-new("corObject",dat.miRNA=as.matrix(miRNA),dat.mRNA=as.matrix(mRNA), pheno.miRNA=pheno.miRNA,pheno.mRNA=pheno.mRNA)
### Plot to check data
plotCordist(data.obj,subset="mRNA",type="dist")
plotCordist(data.obj,subset="miRNA",type="dist")
plotCordist(data.obj,subset="mRNA",type="dist", method.dist = "manhattan")
plotCordist(data.obj,subset="miRNA",type="dist", method.dist = "manhattan")
plotCordist(data.obj,subset="mRNA",type="cor", method.cor="spearman")
plotCordist(data.obj,subset="miRNA",type="cor", method.cor="spearman")
boxplotSamples(data.obj,subset="mRNA")
boxplotSamples(data.obj,subset="miRNA")

#####################################################################
data.obj<-addDiffexp(data.obj,"miRNA",classes="DvH",method.dif="limma")
data.obj<-addDiffexp(data.obj,"mRNA",classes="DvH",method.dif="limma")
data.obj<-addSig(data.obj,"mRNA",adj.pval=0.05,FC=1)
data.obj<-addSig(data.obj,"miRNA",adj.pval=0.05, FC=1.5)
data.obj<-addCorrelation(data.obj,method = "pearson",alternative="less")

data.obj<-addNet(data.obj)
data(microCosm_v5_18)
data(targetScan_v6.2_18)
data.obj<-addDatabase(data.obj,database=c("microCosm_v5_18","targetScan_v6.2_18"))

data.obj<-correctPval(data.obj, pval="pval",method.adj="BH")
data.obj <- addFoldchanges(data.obj)
data.obj<-addScore(data.obj)

plotHeatmap(data.obj,"mRNA", col.color  = 1, colors = c("green4", "red4"))
plotHeatmap(data.obj,"miRNA")
plotPca(data.obj,subset="mRNA", colors = c("green4", "red4"))
plotPca(data.obj,subset="miRNA", colors =   c("green4", "red4"))
plot3d(data.obj,"miRNA", angle = 45, colors = c("green4", "red4"))
plot3d(data.obj,"mRNA", angle = 45, colors = c("green4", "red4"))
##########################################
## I make function plot3Dpca to change color
plot3Dpca(obj = data.obj, "miRNA", angle = 70)
plot3Dpca(obj = data.obj, "mRNA", angle = 45)
##############################################
temp <- data.obj@net
net.sig <- temp[temp$pval < 0.01 & abs(temp$cor) >= 0.8,]
##############################################
library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
df1 <- data.frame(seqnames=names(hm450),
                   starts=start(hm450),
                   ends=end(hm450),
                   names=c(rep(".", length(hm450))),
                   strand=strand(hm450))
probenames <- df1$seqnames
probes <- hm450[probenames]
TSS.probes <- getNearestTSS(probes)
TSS.probes$queryHits <- NULL; TSS.probes$subjectHits <- NULL
TSS.probes.1.5Kb <- TSS.probes[which(TSS.probes$distance <= 1500),]
##################################################
diffMeth <- read.csv("../Correlation/Merged/results.dms.txt", row.names = 1)
net.sig.diffMeth <- net.sig[net.sig$mRNA%in%diffMeth$gene,]
diffMeth_TSS <- diffMeth[rownames(diffMeth)%in% rownames(TSS.probes.1.5Kb),]
net.sig.TSS <- net.sig[net.sig$mRNA%in%diffMeth_TSS$gene,]
Spearman.TSS.sig <- read.csv("Spearman.TSS.sig.txt", header = TRUE, row.names = 1)## TSS probes and gene correlation
Spearman.TSS.sig.neg <- Spearman.TSS.sig[Spearman.TSS.sig$cor.value < 0,]
net.sig.TSS.neg <- net.sig[net.sig$mRNA%in%Spearman.TSS.sig.neg$gene,]

plotCorrelation(data.obj,miRNA="hsa-mir-101-1",mRNA="SKA3",type="cor",col.color="group",sample.names=FALSE, colors = c("green4", "red4"))## Figure6
boxplotCorrelation(data.obj,miRNA="hsa-mir-101-1",mRNA="SKA3", colors = c("red4", "green4"))## Figure *
######################################################################
## net.sig and DiffMerhylated in TSS1.5 and negative spearman's correlation
miRNA_DNAmeth_Spearman <- unique(intersect(net.sig$mRNA, Spearman.TSS.sig.neg$gene))
eg <- bitr(miRNA_DNAmeth_Spearman, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
eg <- eg[order(eg$ENTREZID, decreasing = TRUE),]
kk <- enrichKEGG(gene = eg$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
ggo <- groupGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", level=3, readable = TRUE)
barplot(kk, showCategory = 3)## only 3 pathway enriched
#####################################################################

eg1 <- bitr(unique(net.sig$mRNA), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
kk2 <- enrichKEGG(gene = eg1$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
barplot(kk2, showCategory = 10)
####################################################################
## ne.sig and DiffMet in TSS 1.5 kb
miRNA_DNAmeth <- unique(intersect(net.sig$mRNA, diffMeth_TSS$gene))
eg2 <- bitr(miRNA_DNAmeth, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
kk3 <- enrichKEGG(gene = eg2$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
################################
rm(seqinfo.hg19, hm450, hm450.controls, df1, probenames, probes, TSS.probes, v, first_time)
save.image("mirComb.RData")
