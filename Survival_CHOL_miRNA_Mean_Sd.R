library(pacman)
p_load('TCGAbiolinks', 'survival', 'survminer')

setwd("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/miRNA/")
load("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffExp/miRNA.RData")
clin <- GDCquery_clinic("TCGA-CHOL","clinical")
rna <- log2(cnts+1)
cancerID <- grep("-01A", colnames(rna))
normalID <- grep("-11A", colnames(rna))
rna <- cbind(rna[,cancerID], rna[,normalID])
rm(list= ls()[!(ls() %in% c("rna", "diffExp", "diffExp.voom","clin","miRNAExp"))])
#rownames(rna) <- sapply(strsplit(rownames(rna),"\\|"),'[[',1)
common_names <- rownames(diffExp)

# Create directory (if it doesn't exist)
if (!dir.exists("Survival_Plot_miRNA_Mean")) dir.create("Survival_Plot_miRNA_Mean")   ### change directory name here

### R function for the survival analysis
plot_surv <- function (clinical_patient, dataGE, Genelist, Survresult = FALSE, Median = TRUE, Mean = TRUE, Sdev=0, p.cut = 0.05) 
{
  Genelist <- intersect(rownames(dataGE), Genelist)
  group1 <- colnames(dataGE[,grep("-11A", colnames(dataGE))])
  group2 <- colnames(dataGE[,grep("-01A", colnames(dataGE))])
  dataNormal <- dataGE[Genelist, group1, drop = FALSE]
  dataCancer <- dataGE[Genelist, group2, drop = FALSE]
  colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)
  cfu <- clinical_patient[clinical_patient[, "bcr_patient_barcode"] %in% substr(colnames(dataCancer), 1, 12), ]
  if ("days_to_last_followup" %in% colnames(cfu)) 
    colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <- "days_to_last_follow_up"
  cfu <- as.data.frame(subset(cfu, select = c("bcr_patient_barcode", "days_to_death", "days_to_last_follow_up", "vital_status")))
  if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 0) 
    cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), "days_to_death"] <- "-Inf"
  if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 0) 
    cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), "days_to_last_follow_up"] <- "-Inf"
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  followUpLevel <- FALSE
  tabSurv_Matrix <- matrix(0, nrow(as.matrix(rownames(dataNormal))), 10)
  colnames(tabSurv_Matrix) <- c("mRNA", "pvalue", "Cancer Deaths", "Cancer Deaths with Top", "Cancer Deaths with Down", "Mean Tumor Top", "Mean Tumor Down", "Mean Normal", "Surv Mean Top Months", "Surv Mean Down Months")
  tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)
  cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
  cfu$days_to_last_follow_up <- as.numeric(as.character(cfu$days_to_last_follow_up))
  rownames(cfu) <- cfu[, "bcr_patient_barcode"]
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  cfu_complete <- cfu
  ngenes <- nrow(as.matrix(rownames(dataNormal)))
  for (i in 1:nrow(as.matrix(rownames(dataNormal)))) {
    cat(paste0((ngenes - i), "."))
    mRNAselected <- as.matrix(rownames(dataNormal))[i]
    mRNAselected_values <- dataCancer[rownames(dataCancer) == mRNAselected, ]
    mRNAselected_values_normal <- dataNormal[rownames(dataNormal) == mRNAselected, ]
    if (all(mRNAselected_values == 0)) 
      next
    tabSurv_Matrix[i, "mRNA"] <- mRNAselected
    mRNAselected_values_ordered <- sort(mRNAselected_values, decreasing = TRUE)
    if(Median==TRUE){
      sd <- sd(as.numeric(mRNAselected_values_ordered))
      mRNAselected_values_ordered_top <- median(as.numeric(mRNAselected_values_ordered))+ (sd*Sdev)
      mRNAselected_values_ordered_down <- median(as.numeric(mRNAselected_values_ordered))- (sd*Sdev)
    }
    if(Mean==TRUE){
      sd <- sd(as.numeric(mRNAselected_values_ordered))
      mRNAselected_values_ordered_top <- mean(as.numeric(mRNAselected_values_ordered))+ (sd*Sdev)
      mRNAselected_values_ordered_down <- mean(as.numeric(mRNAselected_values_ordered))- (sd*Sdev)
    }
    mRNAselected_values_newvector <- mRNAselected_values
    if (!is.na(mRNAselected_values_ordered_top)) {
      numberOfSamples <- length(mRNAselected_values_ordered)
      lastelementTOP <- max(which(mRNAselected_values_ordered > mRNAselected_values_ordered_top))
      firstelementDOWN <- min(which(mRNAselected_values_ordered <= mRNAselected_values_ordered_down))
      samples_top_mRNA_selected <- names(mRNAselected_values_ordered[1:lastelementTOP])
      samples_down_mRNA_selected <- names(mRNAselected_values_ordered[firstelementDOWN:numberOfSamples])
      samples_UNCHANGED_mRNA_selected <- names(mRNAselected_values_newvector[which((mRNAselected_values_newvector) > 
                                                                                     mRNAselected_values_ordered_down & mRNAselected_values_newvector < 
                                                                                     mRNAselected_values_ordered_top)])
      cfu_onlyTOP <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_top_mRNA_selected, ]
      cfu_onlyDOWN <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_down_mRNA_selected, ]
      cfu_onlyUNCHANGED <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected, ]
      cfu_ordered <- NULL
      cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
      cfu <- cfu_ordered
      ttime <- as.numeric(cfu[, "days_to_death"])
      sum(status <- ttime > 0)
      deads_complete <- sum(status <- ttime > 0)
      ttime_only_top <- cfu_onlyTOP[, "days_to_death"]
      deads_top <- sum(ttime_only_top > 0)
      if (dim(cfu_onlyDOWN)[1] >= 1) {
        ttime_only_down <- cfu_onlyDOWN[, "days_to_death"]
        deads_down <- sum(ttime_only_down > 0)
      }
      else {deads_down <- 0 }
      tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
      tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
      tabSurv_Matrix[i, "Cancer Deaths with Down"] <- deads_down
      tabSurv_Matrix[i, "Mean Normal"] <- mean(as.numeric(mRNAselected_values_normal))
      dataCancer_onlyTop_sample <- dataCancer[, samples_top_mRNA_selected, drop = FALSE]
      dataCancer_onlyTop_sample_mRNASelected <- dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == mRNAselected, ]
      dataCancer_onlyDown_sample <- dataCancer[, samples_down_mRNA_selected, drop = FALSE]
      dataCancer_onlyDown_sample_mRNASelected <- dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == mRNAselected, ]
      tabSurv_Matrix[i, "Mean Tumor Top"] <- mean(as.numeric(dataCancer_onlyTop_sample_mRNASelected))
      tabSurv_Matrix[i, "Mean Tumor Down"] <- mean(as.numeric(dataCancer_onlyDown_sample_mRNASelected))
      ttime[!status] <- as.numeric(cfu[!status, "days_to_last_follow_up"])
      ttime[which(ttime == -Inf)] <- 0
      ttime <- ttime*0.032854884083862 ## conver days in months
      tabSurv_Matrix[i, "Surv Mean Top Months"] <- mean(ttime[1:nrow(cfu_onlyTOP)])
      tabSurv_Matrix[i, "Surv Mean Down Months"] <- mean(tail(ttime, nrow(cfu_onlyDOWN)))
      ttime <- Surv(ttime, status)
      rownames(ttime) <- rownames(cfu)
      legendHigh <- paste(mRNAselected, "High")
      legendLow <- paste(mRNAselected, "Low")
      legendHigh <- paste0(legendHigh, " (",nrow(cfu_onlyTOP), ")")
      legendLow <- paste0(legendLow, " (",nrow(cfu_onlyDOWN), ")")
      print(paste0("Now running analysis for the Gene :: ",mRNAselected))
      tabSurv_pvalue <- tryCatch({
        tabSurv <- survdiff(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN))))
        tabSurv_chis <- unlist(tabSurv)$chisq
        tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
      }, error = function(e) {
        return(Inf)
      })
      tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue
      if (Survresult == TRUE) {
        
        titlePlot <- paste("Kaplan-Meier Survival analysis, pvalue=", round(tabSurv_pvalue, 3))
        plot(survfit(ttime ~ c(rep("high", nrow(cfu_onlyTOP)), rep("low", nrow(cfu_onlyDOWN)))), mark.time=TRUE, col = c("red4", "blue4"),  main = titlePlot, xlab = "Overall Survival Time", ylab = "Survival Probability", lwd = 1)
        legend("bottomleft", bty="n", lty=1, lwd=1.5, cex=0.8, legend = c(legendHigh, legendLow), col = c("red4","blue4"), text.col = c("red4","blue4"), pch = 15)
        #print(tabSurv)
        file_name = paste0("Survival_Plot_miRNA_Mean/KM_Plot_", mRNAselected, ".pdf")
        dev.print(pdf, file_name, width = 8, height = 8)
      }
    }
  }
  tabSurv_Matrix[tabSurv_Matrix == "-Inf"] <- 0
  tabSurvKM <- tabSurv_Matrix
  tabSurvKM <- tabSurvKM[tabSurvKM$mRNA != 0, ]
  tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < p.cut, ]
  tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA), ]
  rownames(tabSurvKM) <- tabSurvKM$mRNA
  tabSurvKM <- tabSurvKM[, -1]
  tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing = FALSE),]
  return(tabSurvKM)
}
############################################################################

Genelist <- rownames(rna)
results <- plot_surv(clin, dataGE = rna, Genelist = common_names, Survresult = TRUE, p.cut = 0.05, Mean = TRUE, Sdev = 0)
unlink("Survival_Plot_miRNA_Mean/KM_Plot_*") 
results1 <- plot_surv(clin, dataGE = rna, Genelist = rownames(results), Survresult = TRUE, p.cut = 0.05, Mean = TRUE, Sdev = 0)
#results1 <- results[results$`Cancer Deaths with Top` >=3 & results$`Cancer Deaths with Down` >=3,]
results1.DEG <- results1[rownames(results1)%in%common_names,]
#Genelist1 <- Genelist[1:500]
#results1 <- plot_surv(clin, dataGE = rna, Genelist = Genelist, Survresult = TRUE, p.cut = 0.01, Median = TRUE, Sdev = 0) # I already defines group1= normal and group2=cancer

save.image("Survival_CHOL_miRNA_Mean_Sd.RData")
