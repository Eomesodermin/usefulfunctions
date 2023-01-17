#' @title Basic Heatmap function
#' @author Dillon corvino
#' @description Heatmap function
#'
#'
#' @param df dataframe
#' @param Col.cluster logical should columns be clustered
#' @param Row.cluster logical should rows be clustered
#' @param scale.var character either "row", "column", or "none"
#' @param title.var character for what tilte to give the plotted heatmap
#' @param colsepvar integer
#' @param row.size numeric
#' @param dist.method character for which distance metric to use available options are "euclidean", "maximum", "manhattan", "canberra", "binary"
#' @param hclust.method character for which cluster metric to use available options are "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
#' @return heatmap
#'
#' @examples
#' basic.heatmap(df,
#'              Col.cluster = FALSE,
#'              Row.cluster = FALSE,
#'              title.var = "",
#'              colsepvar = NULL,
#'              row.size = 0.4,
#'              dist.method = "euclidean",
#'              hclust.method = "ward.D")
#' @export


basic.heatmap <- function(df,
                          Col.cluster = FALSE,
                          Row.cluster = FALSE,
                          title.var = "",
                          colsepvar = NULL,
                          row.size = 0.4,
                          scale.var = "row",
                          dist.method = "euclidean",
                          hclust.method = "ward.D"){

  dissimfun <- function(x) {
    dist(x, method = paste0(dist.method))
  }

  clusterfun <- function(x) {
    hclust(x, method = paste0(hclust.method))
  }


  gplots::heatmap.2(as.matrix(log2(df + 1)),
                    col = grDevices::colorRampPalette(c("blue", "white", "red"))(100), #turquoise
                    scale = scale.var,
                    na.rm = TRUE,
                    trace = "none",
                    Rowv = Row.cluster,
                    Colv = Col.cluster,
                    distfun = dissimfun,
                    hclustfun = clusterfun,
                    dendrogram = "both",
                    margins = c(10, 7),
                    cexCol = 0.8,
                    cexRow = row.size,
                    main = title.var,
                    key = TRUE,
                    keysize = 1.3,
                    na.color = "yellow",
                    sepcolor = "black",
                    colsep = colsepvar,
                    srtCol = 90)

}

