
#' @title Gene ontology
#' @author Dillon Corvino wrote the wrapper for clusterProfiler function
#' @description
#' Function is designed to take only up or downregulated genes
#' Function will find GO terms enriched in the list of up or downregulated genes
#' input requires a column "gene" and a column "FDR"
#' gene should be as SYMBOL
#'
#'
#' @param markers a dataframe with a column "gene" and a column "FDR"
#' @param topn number of top genes to be used for GO enrichment
#' @param org character for which organism the genes come from either "human" or "mouse"
#' @param title.var character to use as title of the plot
#' @param plot logical, if TRUE, plot the result, if FALSE, return the result
#' @param ... other parameters to be passed to enrichGO function
#'
#' @return if plot is TRUE, return a plot, if plot is FALSE, return a dataframe
#'
#' @export



GO.function <- function(markers,
                        topn = 1000,
                        org = "human",
                        title.var = "",
                        plot = TRUE,
                        ...){



  gene_list <- markers %>%
    dplyr::arrange(.data$FDR) %>%
    head(topn) %>%
    dplyr::pull(.data$gene)

  if(org == "human"){
    db <- org.Hs.eg.db::org.Hs.eg.db
  }else{
    db <- org.Mm.eg.db::org.Mm.eg.db}

  res <- clusterProfiler::enrichGO(gene = gene_list,
                                   OrgDb = db,
                                   keyType = "SYMBOL", ...)

  df <- tibble::as_tibble(res@result) %>%
    dplyr::arrange(.data$p.adjust) %>%
    head(10) %>%
    dplyr::mutate(Description = as.factor(.data$Description)) %>%
    dplyr::mutate(Description = forcats::fct_reorder(.data$Description,
                                                     dplyr::desc(.data$p.adjust)))


  if(plot){
    ggplot(df,
           mapping = aes(x = .data$Description,
                         y = -log10(.data$p.adjust))) +
      geom_bar(aes(fill = .data$Count),
               stat = "identity") +
      scale_fill_gradient2("Gene Count",
                           low = "lightgrey",
                           mid = "#feb24c",
                           high = "#bd0026") +
      coord_flip() +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed") +
      xlab("Gene Ontology") +
      ylab(bquote("-log"[10] ~ " adjusted p-value")) +
      theme_bw() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 12)) +
      ggtitle(title.var)
  }else{

    return(res)

  }
}



