library(miRLAB)
library(TCGAbiolinks)
## CHOL.miRLAB is file for miRLAB analysis
## This file have 122 DE miRNA and 4621 DE genes
cause=1:122
effect=123:4743
PearsonCHOL.1 <- Pearson("CHOL-miRLAB-Logit.csv", cause=1:122, effect=123:4743, targetbinding = "TargetScan.csv")
Top2000 <- Extopk(PearsonCHOL.1, 2000)
#IDACHOL <- IDA("CHOL-miRLAB.csv", cause=1:122, effect=123:4743, "stable", 0.01)
#Top1000.GO <- GOBPenrichment(unique(Top1000[,2]), 0.05)
GeneList <- unique(Top2000[,2])
system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Tumor Vs Normal",GeneList))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = GeneList, 
                        nBar = 10)