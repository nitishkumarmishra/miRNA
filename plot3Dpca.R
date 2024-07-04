## plot3Dpca.R# R program for 3d PCA plot in Red and green Color
library(miRComb)
plot3Dpca <- function (obj, subset, col.color = 1, angle = 45, colors = c("green4", 
                                                             "red4"), ...) 
{
  nf <- layout(mat = matrix(c(1), ncol = 1, nrow = 1))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  if (subset == "miRNA") {
    color <- as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[, 
                                                                      col.color])) + 1, labels = colors))
    groups <- levels(as.factor(obj@pheno.miRNA[, col.color]))
    sel <- which(apply(obj@dat.miRNA, 1, sd) == 0)
    if (length(sel) > 0) {
      pca <- prcomp(t(obj@dat.miRNA[-sel, ]), scale = T)
    }
    else {
      pca <- prcomp(t(obj@dat.miRNA), scale = T)
    }
    labels.dat <- colnames(obj@dat.miRNA)
  }
  if (subset == "mRNA") {
    color <- as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[, 
                                                                     col.color])) + 1, labels = colors))
    groups <- levels(as.factor(obj@pheno.mRNA[, col.color]))
    sel <- which(apply(obj@dat.mRNA, 1, sd) == 0)
    if (length(sel) > 0) {
      pca <- prcomp(t(obj@dat.mRNA[-sel, ]), scale = T)
    }
    else {
      pca <- prcomp(t(obj@dat.mRNA), scale = T)
    }
    labels.dat <- colnames(obj@dat.mRNA)
  }
  pca$pcolor <- color
  with(pca, {
    s3d <- scatterplot3d(pca$x[, 1], pca$x[, 2], pca$x[, 
                                                       3], color = pcolor, pch = 19, type = "h", lty.hplot = 0, 
                         cex.axis = 1.5, cex.lab = 1.5, cex.symbols=1.5, scale.y = 0.75, xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2) * 
                                                                          100, 1), "%", sep = ""), ylab = paste("Comp 2: ", 
                                                                                                                round(pca$sdev[2]^2/sum(pca$sdev^2) * 100, 1), 
                                                                                                                "%", sep = ""), zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2) * 
                                                                                                                                                                 100, 1), "%", sep = ""), angle = angle, ...)
    s3d.coords <- s3d$xyz.convert(pca$x[, 1], pca$x[, 2], 
                                  pca$x[, 3])
    legend("topleft", inset = 0.05, bty = "n", cex = 1.2, title = "Group", 
           groups, col = c(levels(as.factor(color))), pch = 19)
  })
}