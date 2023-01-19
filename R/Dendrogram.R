
#' @title Dendrogram function
#' @description Dendrogram function
#'
#' @param data.matrix A matrix of data
#' @param hc_metric.var A string of the distance metric to be used in the hierarchical clustering options include "euclidean", "maximum", "manhattan", "canberra", "binary"
#' @param hc_method.var A string of the clustering method to be used in the hierarchical clustering options include "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
#' @param type.var A string of the type of dendrogram to be plotted options include "phylogram”, “cladogram”, “fan”, “unrooted”, “radial”, if null, standard dendrogram plotted
#' @param output.dir A string of the output directory
#' @return saves a pdf file of a dendrogram

#' @examples
#' dendrogram(data.matrix = data.matrix,
#'                    hc_metric.var = "euclidean",
#'                    hc_method.var = "ward.D2",
#'                    type.var = "unrooted",
#'                    output.dir = "results/figures/")

#' @export

# Dendrogram function
dendrogram <- function(data.matrix,
                               hc_metric.var = "euclidean",
                               hc_method.var = "ward.D2",
                               type.var = NULL,
                               output.dir = "results/figures/"){

   # Create output.dir
  if(!dir.exists(output.dir)){
    dir.create(output.dir,
               recursive = TRUE)}



  # set plotting margins
  par(mar=c(6,4.1,4.1,2.1))


  #Sample Dendrogram
  res.hc <- factoextra::eclust(t(data.matrix),
                               stand = TRUE,
                               FUNcluster = "hclust",
                               hc_metric = hc_metric.var,
                               hc_method = hc_method.var,
                               k = 1,
                               verbose = TRUE,
                               graph = FALSE)



  x <- stats::as.dendrogram(res.hc)

  pdf(paste0(output.dir, "Dendrogram.pdf"))

  print(factoextra::fviz_dend(x,
                  show_labels = TRUE,
                  color_labels_by_k = FALSE,
                  type = "rectangle",
                  rect = TRUE,
                  rect_border = c("black"),
                  rect_lty = "solid",
                  rect_lwd = 1,
                  main = "Sample Dendrogram",
                  xlab = paste0("Dist = ", hc_metric.var, " & Clust = ",  hc_method.var),
                  cex = 0.6))
  dev.off()


  if(!is.null(type.var)){
  pdf(paste0(output.dir, "Dendrogram_", type.var, ".pdf"))

  print(plot(ape::as.phylo(res.hc),
             type = type.var,
             cex = 0.6,
             no.margin = TRUE))

  dev.off()
  }


}





