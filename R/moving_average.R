
#' @title Moving average plot
#' @description This function plots the moving average of a gene of interest (GOI) against the expression of another GOI.
#'
#'
#' @param tcga.datasets A vector of TCGA datasets to use for the analysis.
#' @param output.dir The directory to save the output files.
#' @param primary.goi The primary GOI to use for ranking the samples.
#' @param goi A vector of GOI to plot against the primary GOI.
#' @param sample.width The sample window size for the moving average calculation.
#' @param laxis.col The colour of the dependent gene.
#' @return A plot of the moving average of a gene of interest (GOI) against the expression of another GOI.
#' @export
#' @examples
#'
#' moving.average()



moving.average <- function(tcga.datasets = c("laml_tcga", "luad_tcga"),
                           output.dir = "results/moving_average/",
                           primary.goi = "TLR2",
                           goi = c("TLR6", "CSF1R", "CD14", "CD8B", "NCAM1", "IFNG", "GZMB"),
                           sample.width = 20, # change sample window size for moving average calculation
                           laxis.col = "blue"  #change colour of dependent gene
){


  # load required packages
  #require("ggplot2")
  #require("gplots")
  #require("TCGAbiolinks")
  #require("SummarizedExperiment")
  #require("car")
  #require("psych")
  #require("GenomeInfoDb")
  #require("Biobase")
  #require("cgdsr")
 # require("GenomicDataCommons")

  #GenomicDataCommons::status()


  # ensure parent output dir is made
  # Create output directory
  if(!dir.exists(paste0(output.dir))){
    dir.create(paste0(output.dir), recursive = TRUE)}




  for(i in seq_along(tcga.datasets)){

    # which TCGA dataset to use
    tcga.cancer <- paste0(tcga.datasets[i])

    # Create output directory
    output.dir.temp <- paste0(output.dir, tcga.cancer)

    if(!dir.exists(paste0(output.dir.temp))){
      dir.create(paste0(output.dir.temp), recursive = TRUE)}


    # Get tcga data using CGDS R package
    # create CGDS object and connection to cBioportal
    mycgds <- cgdsr::CGDS("http://www.cbioportal.org/", verbose = T)

    # test Bioportal connection
    #test(mycgds)

    # Get gene expression data using getProfileData
    genetic_profile_id <-paste0(tcga.cancer,"_rna_seq_v2_mrna")
    case_list_id <- paste0(tcga.cancer, "_rna_seq_v2_mrna")


    # enter genes of interest to retrieve via cBioportal
    # Vector of genes of interest
    genes <- unique(c(primary.goi, goi))

    # get gene expression data frame
    gedf <- cgdsr::getProfileData(mycgds, genes, genetic_profile_id, case_list_id)

    # check normalized gene expression data (read counts) and log2 transform (avoid negative values)
    dim(gedf)
    ncol(gedf) == length(genes) # TRUE
    head(gedf)
    gedf[gedf < 1] <- 1
    gedf <- log2(gedf)
    head(gedf)
    n.val <- nrow(gedf)

    # define gene used for ranking of samples
    xSortGene <- gedf[[primary.goi]]
    xSortGeneID <- primary.goi

    xsort <- sort(xSortGene, decreasing = FALSE, index.return = TRUE)



    # Plot moving average for all GOI except first (CD226)
    for(gene.i in seq_along(genes)){

      if(primary.goi == genes[gene.i]){

        print("no need to plot primary.goi vs itself")
      }else{

        # define dependent gene of interest
        xDepGene <- gedf[ , genes[gene.i]]
        xDepGeneID <- genes[gene.i]

        xLeftGene <- xDepGene[xsort$ix] - median(xDepGene[xsort$ix])


        par(mar = c(5.1, 4.1, 4.1, 4.2))


        x <- xLeftGene

        ylim.plot <- c(summary(xLeftGene)[1], summary(xLeftGene)[6]) # change scale of dependent gene y-axis

        # Plot bars of dependent gene expression as a background trace
        plot(x,
             col = "lightgrey",
             type = "h",
             lwd = 0.5,
             ylim = ylim.plot,
             xlab = "",
             ylab ="",
             axes = F)
        box(which = "plot")
        axis(side = 1, outer = F)

        # Plot moving average line for the dependent gene
        axis(side = 2,
             col = laxis.col,
             col.axis = laxis.col,
             lwd = 2)
        abline(h = 0,
               col = "lightgrey",
               lwd = 0.5)
        lines(usefulfunctions::calc.moving.average(x,
                                  n = sample.width,
                                  centered =T),
              col = laxis.col,
              type = "l",
              lwd = 3)
        par(new=T)

        # Plot line for the ranking gene
        plot(xsort$x - median(xsort$x),
             col = "black",
             type = "l",
             lwd = 3,
             xlab = "",
             ylab ="",
             axes = F)
        axis(side = 4,
             col = "black",
             col.axis = "black",
             lwd = 2)

        # Add axes labels
        title(main = paste0(xDepGeneID, " in ", tcga.cancer),
              col.main = laxis.col,
              xlab = paste0("TCGA samples ranked by", primary.goi, "expression"),
              ylab = paste0("Moving average ", xDepGeneID),
              col.lab = laxis.col)

        mtext(paste0(xSortGeneID, " expression (log2)"), side = 4, line = 2.3)


        # Correlation values
        cor.val <- cor(xLeftGene, xsort$x, method = "spearman")
        cor.val <- signif(cor.val, 4)
        mtext(paste0("r = ", cor.val), adj = 0)


        cor.test.val <- cor.test(xLeftGene, xsort$x, method = "spearman")
        corr.pval <- signif(cor.test.val$p.value, 4)
        mtext(paste0("P = ", corr.pval), adj = 1)

        mtext(paste0("n = ", n.val), adj = 0.5)

        # Save plot to file
        dev.copy(pdf, paste0(output.dir.temp, "/", primary.goi, "_vs_", xDepGeneID, ".pdf"))
        dev.off()


      }




    }



  }


}





