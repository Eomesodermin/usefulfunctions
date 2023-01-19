
#' @title Convert Mouse genes to human genes
#' @author written by radiaj in R bloggers, Oct 14, 2016 https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
#' @usage convert.mouse.to.human(x = c("Gzma", "Gzmb", "Xcl1"))
#'
#'
#' @param genes a character vector of mouse gene names
#' @param use.jax logic variable to convert genes using data from jax website, to use incase connection error with ensemble
#' @return a data frame with mouse and human gene names
#'
#' @examples
#' convert.mouse.to.human(genes = c("Gzma", "Gzmb", "Xcl1"))
#'
#' @export

convert.mouse.to.human <- function(genes, use.jax = FALSE){

  if(use.jax){
    mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")


    output.df = c()

    for(gene.i in genes){
      class_key = (mouse_human_genes %>% dplyr::filter(.data$Symbol == gene.i & .data$Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
      if(!identical(class_key, integer(0)) ){
        human_genes = (mouse_human_genes %>% dplyr::filter(.data$DB.Class.Key == class_key & .data$Common.Organism.Name=="human"))[,"Symbol"]
        for(human_gene in human_genes){
          output.df = append(output.df, human_gene)
        }
      }
    }

    print(head(output.df))

    return(output.df)



  }else{
  human <- biomaRt::useMart("ensembl",
                           dataset = "hsapiens_gene_ensembl")

  mouse <- biomaRt::useMart("ensembl",
                           dataset = "mmusculus_gene_ensembl")

  genes.new <- biomaRt::getLDS(attributes = c("mgi_symbol"),
                            filters = "mgi_symbol",
                            values = genes,
                            mart = mouse,
                            attributesL = c("hgnc_symbol"),
                            martL = human,
                            uniqueRows=T)

  output.df <- genes.new %>%
    dplyr::distinct(.data$MGI.symbol, .keep_all = TRUE) %>%
    dplyr::distinct(.data$HGNC.symbol, .keep_all = TRUE)

  print(head(output.df))

  return(output.df)
  }
}


#' @title Convert Human genes to mouse genes
#' @author written by radiaj in R bloggers, Oct 14, 2016 https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
#' @usage convert.human.to.mouse(x = c("GZMA", "GZMB", "XCL1"))
#'
#'
#' @param genes a character vector of human gene names
#' @param use.jax logic variable to convert genes using data from jax website, to use incase connection error with ensemble
#' @return a data frame with mouse and human gene names
#'
#' @examples
#' convert.human.to.mouse(genes = c("GZMA", "GZMB", "XCL1"))
#'
#' @export

convert.human.to.mouse <- function(genes, use.jax = FALSE){
  if(use.jax){

    mouse_human_genes = utils::read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

    output.df = c()

    for(gene.i in genes){
      class_key = (mouse_human_genes %>% dplyr::filter(.data$Symbol == gene.i & .data$Common.Organism.Name=="human"))[['DB.Class.Key']]
      if(!identical(class_key, integer(0)) ){
        mouse_genes = (mouse_human_genes %>% dplyr::filter(.data$DB.Class.Key == class_key & .data$Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
        for(mouse_gene in mouse_genes){
          output.df = append(output.df, mouse_gene)
        }
      }
    }

    print(head(output.df))

    return(output.df)

  }else{
  human <- biomaRt::useMart("ensembl",
                           dataset = "hsapiens_gene_ensembl")

  mouse <- biomaRt::useMart("ensembl",
                           dataset = "mmusculus_gene_ensembl")

  genes.new <- biomaRt::getLDS(attributes = c("hgnc_symbol"),
                            filters = "hgnc_symbol",
                            values = genes,
                            mart = human,
                            attributesL = c("mgi_symbol"),
                            martL = mouse,
                            uniqueRows=T)

  output.df <- genes.new %>%
    dplyr::distinct(.data$HGNC.symbol, .keep_all = TRUE) %>%
    dplyr::distinct(.data$MGI.symbol, .keep_all = TRUE)

  print(head(output.df))

  return(output.df)
  }

}



