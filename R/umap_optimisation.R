#' @title UMAP optimisation function
#' @description This function will calculate UMAP for a range of min.dist and n.neighbor values.
#'
#'
#' @param input.seurat A Seurat object
#' @param output.dir A directory to save the output to. If NULL, will create a directory called "output/Optimising_UMAP"
#' @param init.min.dist The starting min.dist value. sensible values = 0.001 - 0.5
#' @param init.neigh The starting n.neighbor value. sensible values = 5 - 50
#' @param min.dist.step The step size for min.dist. Default = 2
#' @param neigh.step The step size for n.neighbor. Default = 20
#' @return
#' Numerous .pdf files are exported and saved to the desired output.dir

#' @examples
#' UMAP.optimise(input.seurat = pbmc_small,
#'              output.dir = NULL,
#'              init.min.dist = 0.001,
#'              init.neigh = 5,
#'              min.dist.step = 2,
#'              neigh.step = 20)
#' @export

UMAP.optimise <- function(input.seurat,
                          output.dir = "results/optimising_UMAP/",
                          init.min.dist = 0.001,
                          init.neigh = 5,
                          min.dist.step = 2,
                          neigh.step = 20){


  # Create directory if one isnt already created
  if(!dir.exists(output.dir)){
    dir.create(output.dir,
               recursive = TRUE)}


  # set starting min.dist
  min.dist.val <- init.min.dist

  while(min.dist.val < 1){

    # set starting neighbor value
    n.neigh.val <- init.neigh

    while(n.neigh.val <= 50){

      print(paste0("Calculating UMAP for min.dist = ", min.dist.val, " & n.neighbor = ", n.neigh.val))

      # Calculate UMAP
      temp <- Seurat::RunUMAP(object = input.seurat,
                              reduction = "pca",
                              dims = 1:20,
                              umap.method = "uwot",
                              n.neighbors = n.neigh.val,
                              min.dist = min.dist.val,
                              seed.use = 42)

      # Plot UMAP
      print(Seurat::UMAPPlot(object = temp,
                             label = TRUE,
                             label.size = 4) +
              ggtitle(paste0("UMAP Min dist = ", min.dist.val, " n.val = ", n.neigh.val)) +
              NoLegend())

      dev.copy(pdf, paste0(output.dir, "UMAP_Min_dist_", min.dist.val, "_neighval_", n.neigh.val, ".pdf"))
      dev.off()

      n.neigh.val <- n.neigh.val + neigh.step
    }
    min.dist.val <- min.dist.val * min.dist.step
  }

  rm(temp)

}








