

# Define correlation functions - adapted from Hoelzel code
#' @title probe correlation
#' @description calculate correlation between goi and data
#' @author Adapted from Hoelzel code
#' @param em a matrix of expression data
#' @param gp string - gene of interest to prob correlation against
#' @param method a string
#' @return a correlation matrix
#' @examples
#' cor.res <- probeCor(filt.data,
#'                     gp = "NKG7",
#'                     method = "pearson")
#'
probeCor <- function(em, gp, method = "pearson") {
  prbcor <- apply(em, 1, stats::cor, em[gp,], method = method)
  return(prbcor)
}

#' @title Prob corelation test
#' @description
#' This function is used to calculate the statistical significance of correlation
#' @author Adapted from Hoelzel code
#' @param em a matrix of expression data
#' @param gp string - gene of interest to prob correlation against
#' @param method a string of method used to calculate the correlation
#' @return a list of correlation test result
#' @examples
#'cor.res.pval <- probeCor.test(filt.data,
#'                              gp = "NKG7",
#'                              method = "pearson")

probeCor.test <- function(em, gp, method = "pearson") {
  prbcort <- apply(em, 1, stats::cor.test, em[gp,], method = method)
  return(prbcort)
}



#' @title Correlation analysis
#' @description Calculate correlation between a gene of interest and all other genes in the dataset
#' @param data.slot dataframe from @data slot from a normalised and downsampled seurat object - already downsampled according to ident of interest to a reasonable cell number
#' @param goi gene to test correlation against
#' @param remove.quantile lowest X% expressed genes in dataset to remove from analysis
#' @param output.tables output directory for correlation results
#' @return correlation dataframe and saves a copy to designated output location
#' @export
#' @examples
#'
#'
#' # Set Idents
#' Idents(seurat.object) <- seurat.object@meta.data$cluster
#' # Downsample object as necessary
#' # extract normalised data table
#' normed.data <- as.matrix(seurat.object@assays$RNA@data)
#'
#' # calculate correlation
#' sc.correlation(data.slot = normed.data,
#'               goi = "NKG7",
#'               remove.quantile = 0.3,
#'               output.tables = "results/tables/Correlation_analysis/")

sc.correlation <- function(data.slot,
                           goi = "NKG7",
                           remove.quantile = 0.3,
                           output.tables = "results/tables/Correlation_analysis/",
                           debug.mode = FALSE){


  # Create output directory
  if(!dir.exists(paste0(output.tables))){
    dir.create(paste0(output.tables),
               recursive = T)
  }


  # Identify lowest X% expressed genes in dataset
  q.val <- quantile(rowMeans(data.slot), probs = remove.quantile)
  low.exprs <- rowMeans(data.slot) < q.val

  # Filter low expressed genes
  filt.data <- data.slot[!low.exprs, ]

  # Calculate pearson correlation for all genes vs. NKG7
  cor.res <- probeCor(filt.data,
                      gp = goi,
                      method = "pearson")

  # Sort correlation values
  cor.res <- as.data.frame(sort(cor.res,
                                decreasing = TRUE))

  colnames(cor.res) <- "Cor.val"


  # Calculate statistics
  cor.res.pval <- probeCor.test(filt.data,
                                gp = goi,
                                method = "pearson")

  cor.res.pval <- unlist(sapply(cor.res.pval, "[", "p.value"))

  names(cor.res.pval) <- gsub(".p.value", "", names(cor.res.pval))

  cor.res.pval <- as.data.frame(sort(cor.res.pval,
                                     decreasing = FALSE))

  colnames(cor.res.pval) <- "P.val"


  # Combine both pval and correlation statistic
  cor.df <- merge(cor.res, cor.res.pval, by = "row.names")

  colnames(cor.df)[1] <- c("GeneID")

  # Calculate FDR adjusted P value
  cor.df$FDR <- p.adjust(cor.df$P.val, method = "fdr")


  # Write data to file
  write.csv(cor.df, paste0(output.tables,  "Correlation_values_", goi, ".csv"))

  # return correlation dataframe
  return(cor.df)


}



















#' @title Correlation.heatmaps
#' @description Creates heatmaps of the top correlated genes to a gene of interest, output from sc.correlation function
#' @param seurat.object Seurat object, should be downsampled and idents and assay set to those of interest
#' @param correlation.df Dataframe of correlation values, output from sc.correlation function
#' @param up.n Number of up regulated genes to show
#' @param dn.n Number of downregulated genes to show
#' @param goi Gene of interest - used in calculating correlation
#' @param output.dir Output directory for saving figures
#' @export
#' @examples
#'
#'
#'DefaultAssay(downsampled.seurat) <- "SCT"
#'Idents(downsampled.seurat) <- downsampled.seurat@meta.data$celltype.l3
#'
#'small.seurat <- subset(downsampled.seurat,
#'                       downsample = 100)
#'
#'small.seurat <- ScaleData(small.seurat)
#'
#' correlation.heatmaps(seurat.object = small.seurat,
#'                     correlation.df = cor.df,
#'                     up.n = 20,
#'                     dn.n = 20,
#'                     goi = "NKG7",
#'                     output.dir = "results/figures/Correlation_analysis/")
#' @return print and export heatmaps with batlow.pal colour scheme

correlation.heatmaps <- function(seurat.object,
                                 correlation.df,
                                 up.n = 20,
                                 dn.n = 20,
                                 goi = "NKG7",
                                 output.dir = "results/figures/Correlation_analysis/"){


  # Create output directory
  if(!dir.exists(paste0(output.dir))){
    dir.create(paste0(output.dir),
               recursive = T)
  }


  # Top upregulated genes
  up.genes <- correlation.df %>%
    dplyr::filter(.data$FDR < 0.05) %>%
    dplyr::top_n(up.n, .data$Cor.val) %>%
    dplyr::pull(.data$GeneID)

  up.genes <- unique(c(goi, up.genes))


  print(Seurat::DoHeatmap(seurat.object,
                          features = up.genes,
                          angle = 90,
                          size = 3,
                          raster = FALSE) + scico::scale_fill_scico(palette = "batlow",
                                                                    direction = 1,
                                                                    na.value = "white"))

  dev.copy(pdf, paste0(output.dir, "Heatmap_top", up.n, "_corr_genes.pdf"))
  dev.off()


  # Top downregulated genes

  dn.genes <- correlation.df %>%
    dplyr::filter(.data$FDR < 0.05) %>%
    dplyr::top_n(-dn.n, .data$Cor.val) %>%
    dplyr::pull(.data$GeneID)

  dn.genes <- unique(c(goi, dn.genes))

  # Neg correlated genes
  print(Seurat::DoHeatmap(seurat.object,
                          features = dn.genes) + scico::scale_fill_scico(palette = "batlow",
                                                                         direction = 1,
                                                                         na.value = "white"))

  dev.copy(pdf, paste0(output.dir, "Heatmap_bottom", dn.n, "_corr_genes.pdf"))
  dev.off()


}



#' @title Plot correlation data
#'
#' @param cor.df dataframe of correlation values
#' @param output.dir output directory
#' @param goi gene of interest
#' @param top.n.val number of top values to plot
#'
#' @examples
#' plot.cor.data(cor.df = cor.df,
#'              output.dir = "output/figures/",
#'              goi = "NKG7",
#'              top.n.val = 30)
#' @export
#'
plot.cor.data <- function(cor.df,
                          output.dir = "output/figures/",
                          goi = "NKG7",
                          top.n.val = 30){

  # Create output directory
  if(!dir.exists(paste0(output.dir))){
    dir.create(paste0(output.dir),
               recursive = T)
  }



  # set up dataset and variables

  # order datset
  ordered.cor <- cor.df %>%
    dplyr::arrange(plyr::desc(.data$Cor.val))

  y.val <- nrow(ordered.cor)

  ordered.cor$Ypos <- c(y.val:1)

  # subset data for top X value
  top.corr <- ordered.cor[1:top.n.val, ]


  pdf(paste0(output.dir, "correlation.pdf"), width = 5, height= 9)

  # Plot overview of correlation analysis
  graphics::plot(ordered.cor$Cor.val,
             ordered.cor$Ypos,
             type = "l")

  # Highlight goi
  for(i in seq_along(goi)){
    goi.logic <- ordered.cor$GeneID == goi[i]
    goi.xpos <- ordered.cor$Cor.val[goi.logic]
    goi.ypos <- ordered.cor$Ypos[goi.logic]

    graphics::points(goi.xpos,
           goi.ypos)

    graphics::text(goi.xpos,
         goi.ypos,
         labels = goi[i],
         pos = 1)

  }


  # create line highlighting top n cor vals
  graphics::lines(top.corr$Cor.val,
        max(top.corr$Ypos):min(top.corr$Ypos),
        lwd = 2,
        col = "blue")

  dev.off()


  # Plot an inset of the top X

  # Plot top x
  pdf(paste0(output.dir, "correlation_top_", top.n.val, ".pdf"), width = 3, height= 6)

  graphics::plot(top.corr$Cor.val,
       max(top.corr$Ypos):min(top.corr$Ypos),
       type = "l",
       xlim = c(0.1,1),
       col = "blue",
       lwd = 2)

  graphics::points(top.corr$Cor.val,
         max(top.corr$Ypos):min(top.corr$Ypos),
         pch = 20)


  graphics::text(top.corr$Cor.val,
       max(top.corr$Ypos):min(top.corr$Ypos),
       top.corr$GeneID,
       pos = 4)

  dev.off()


}



