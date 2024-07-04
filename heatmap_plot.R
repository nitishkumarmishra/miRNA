heatmap_plot <- function (obj, class, n = 50, col.color = 1, min.exp = NULL, 
          main = NULL, pval.cutoff = NULL) 
{
  triming <- function(X) {
    num <- X
    for (i in 1:dim(num)[1]) {
      x <- num[i, ]
      trim = 0.05
      lo = quantile(x, trim)
      hi = quantile(x, 1 - trim)
      x[x < lo] = lo
      x[x > hi] = hi
      num[i, ] <- x
    }
    return(num)
  }
  if (length(col.color) == 1) {
    nf <- layout(mat = matrix(c(1), ncol = 1, nrow = 1))
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    colors.cat <- rep(palette()[-1], length.out = nrow(obj@pheno.miRNA))
    col <- c("blue4", "red4")
    if (class == "miRNA") {
      if (length(levels(as.factor(obj@pheno.miRNA[, col.color]))) == 
          2) {
        cc.color <- as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[, 
                                                                             col.color])), labels = col))
      }
      else {
        cc.color <- as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[, 
                                                                             col.color])) + 1))
      }
      data <- obj@diffexp.miRNA
      if (!is.null(min.exp)) {
        data <- data[which(data$meanExp >= min.exp), 
                     ]
      }
      sel.sort <- rownames(data[with(data, order(pval)), 
                                ])
      xmat <- obj@dat.miRNA[sel.sort[1:n], ]
      rownames(xmat) <- gsub("hsa-", "", rownames(xmat))
      xmat <- xmat[which(apply(xmat, 1, sd) > 0), ]
      cc.row <- hclust(dist(xmat, method = "euclidean"))
      cc.col <- hclust(dist(t(xmat), method = "euclidean"))
      if (is.null(main)) {
        main <- paste("Top", n, class)### I can use change from main option
      }
      heatmap.2(triming(xmat), col = greenred(75), scale = "row", 
                ColSideColors = cc.color, key = TRUE, keysize = 1, symkey = FALSE, 
                density.info = "none", trace = "none", cexRow = 0.75, 
                main = main, labCol = NA, dendrogram = "both", 
                distfun = function(x) as.dist(1 - cor(t(x), method = "pearson")), 
                hclustfun = function(x) hclust(x, method = "average"))
      leg <- levels(as.factor(obj@pheno.miRNA[, col.color]))
      legend("topright", leg, col = levels(as.factor(cc.color)), 
             lwd = 7)
    }
    if (class == "mRNA") {
      if (length(levels(as.factor(obj@pheno.mRNA[, col.color]))) == 
          2) {
        cc.color <- as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[, 
                                                                            col.color])), labels = col))
      }
      else {
        cc.color <- as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[, 
                                                                            col.color])) + 1))
      }
      data <- obj@diffexp.mRNA
      if (!is.null(min.exp)) {
        data <- data[which(data$meanExp >= min.exp), 
                     ]
      }
      sel.sort <- rownames(data[with(data, order(pval)), 
                                ])
      xmat <- obj@dat.mRNA[sel.sort[1:n], ]
      xmat <- xmat[which(apply(xmat, 1, sd) > 0), ]
      cc.row <- hclust(dist(xmat, method = "euclidean"))$order
      cc.col <- hclust(dist(t(xmat), method = "euclidean"))$order
      if (is.null(main)) {
        main <- paste("Top", n, class)
      }
      heatmap.2(triming(xmat), col = greenred(75), scale = "row", 
                ColSideColors = cc.color, key = TRUE, keysize = 1, symkey = FALSE, 
                density.info = "none", trace = "none", cexRow = 0.75, 
                main = main, labCol = NA, dendrogram = "both", 
                distfun = function(x) as.dist(1 - cor(t(x), method = "pearson")), 
                hclustfun = function(x) hclust(x, method = "average"))
      leg <- levels(as.factor(obj@pheno.mRNA[, col.color]))
      legend("topright", leg, col = levels(as.factor(cc.color)), 
             lwd = 7)
    }
    if (class == "both") {
      cyto <- obj@net[with(obj@net, order(adj.pval)), ]
      names.miRNA <- unique(cyto[1:n, "miRNA"])
      names.mRNA <- unique(cyto[1:n, "mRNA"])
      xmat.miRNA <- obj@dat.miRNA[names.miRNA, ]
      rownames(xmat.miRNA) <- gsub("hsa-", "", rownames(xmat.miRNA))
      xmat.mRNA <- obj@dat.mRNA[names.mRNA, ]
      xmat <- rbind(xmat.miRNA, xmat.mRNA)
      xmat <- xmat[which(apply(xmat, 1, sd) > 0), ]
      if (length(levels(as.factor(obj@pheno.miRNA[, col.color]))) == 
          2) {
        cc.color <- as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[, 
                                                                             col.color])), labels = col))
      }
      else {
        cc.color <- as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[, 
                                                                             col.color])) + 1))
      }
      if (is.null(main)) {
        main <- paste("Top", n, class)
      }
      cc.row <- hclust(dist(xmat, method = "euclidean"))
      cc.col <- hclust(dist(t(xmat), method = "euclidean"))
      heatmap.2(triming(xmat), col = greenred(75), scale = "row", 
                ColSideColors = cc.color, key = TRUE, keysize = 1, symkey = FALSE, 
                density.info = "none", trace = "none", cexRow = 0.75, 
                main = main, labCol = NA, dendrogram = "both", 
                distfun = function(x) as.dist(1 - cor(t(x), method = "pearson")), 
                hclustfun = function(x) hclust(x, method = "average"))
      leg <- levels(as.factor(obj@pheno.mRNA[, col.color]))
      legend("topright", leg, col = levels(as.factor(cc.color)), 
             lwd = 7)
    }
  }
  if (length(col.color) > 1) {
    if (class == "miRNA") {
      data <- obj@diffexp.miRNA
      if (!is.null(min.exp)) {
        data <- data[which(data$meanExp >= min.exp), 
                     ]
      }
      sel.sort <- rownames(data[with(data, order(pval)), 
                                ])
      xmat <- obj@dat.miRNA[sel.sort[1:n], ]
      rownames(xmat) <- gsub("hsa-", "", rownames(xmat))
      xmat <- triming(xmat)
      xmat <- xmat[which(apply(xmat, 1, sd) > 0), ]
      if (is.null(main)) {
        main <- paste("Top", n, class)
      }
      pheatmap(xmat, annotation_col = data.frame(obj@pheno.miRNA[, 
                                                                 col.color]), color = greenred(75), main = main, 
               scale = "row")
    }
    if (class == "mRNA") {
      data <- obj@diffexp.mRNA
      if (!is.null(min.exp)) {
        data <- data[which(data$meanExp >= min.exp), 
                     ]
      }
      sel.sort <- rownames(data[with(data, order(pval)), 
                                ])
      xmat <- obj@dat.mRNA[sel.sort[1:n], ]
      xmat <- triming(xmat)
      xmat <- xmat[which(apply(xmat, 1, sd) > 0), ]
      if (is.null(main)) {
        main <- paste("Top", n, class)##### I can c
      }
      pheatmap(xmat, annotation_col = data.frame(obj@pheno.mRNA[, 
                                                                col.color]), color = greenred(75), main = main, 
               scale = "row")
    }
    if (class == "both") {
      cyto <- obj@net[with(obj@net, order(adj.pval)), ]
      names.miRNA <- unique(cyto[1:n, "miRNA"])
      names.mRNA <- unique(cyto[1:n, "mRNA"])
      xmat.miRNA <- obj@dat.miRNA[names.miRNA, ]
      rownames(xmat.miRNA) <- gsub("hsa-", "", rownames(xmat.miRNA))
      xmat.mRNA <- obj@dat.mRNA[names.mRNA, ]
      xmat <- rbind(xmat.miRNA, xmat.mRNA)
      xmat <- triming(xmat)
      xmat <- xmat[which(apply(xmat, 1, sd) > 0), ]
      if (is.null(main)) {
        main <- paste("Top", n, class)
      }
      pheatmap(xmat, annotation_col = data.frame(obj@pheno.miRNA[, 
                                                                 col.color]), color = greenred(75), main = main, 
               scale = "row")
    }
  }
}