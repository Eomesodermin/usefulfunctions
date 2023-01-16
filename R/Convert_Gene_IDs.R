
#' @title Convert Mouse genes to human genes
#' @author written by radiaj in R bloggers, Oct 14, 2016 https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
#' @usage convert.mouse.to.human(mouse.genes)
#' 
#' @importFrom dplyr distinct
#' @importFrom biomaRt useMart getLDS
#'
#' @param x a character vector of mouse gene names
#' @return a data frame with mouse and human gene names
#'
#' @examples
#' convert.mouse.to.human(c("Gzma", "Gzmb", "Xcl1"))
#'
#' @export

convert.mouse.to.human <- function(x){
  
  human = biomaRt::useMart("ensembl", 
                           dataset = "hsapiens_gene_ensembl")
  
  mouse = biomaRt::useMart("ensembl", 
                           dataset = "mmusculus_gene_ensembl")
  
  genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"),
                            filters = "mgi_symbol",
                            values = x , 
                            mart = mouse,
                            attributesL = c("hgnc_symbol"),
                            martL = human, 
                            uniqueRows=T)
  
  output.df <- genesV2 %>%
    dplyr::distinct(MGI.symbol, .keep_all = TRUE) %>%
    dplyr::distinct(HGNC.symbol, .keep_all = TRUE)
  
  print(head(output.df))
  
  return(output.df)
}


#' @title Convert Human genes to mouse genes
#' @author written by radiaj in R bloggers, Oct 14, 2016 https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
#' @usage convert.human.to.mouse(human.genes)
#' 
#' @importFrom dplyr distinct
#' @importFrom biomaRt useMart getLDS
#'
#' @param x a character vector of human gene names
#' @return a data frame with mouse and human gene names
#'
#' @examples
#' convert.human.to.mouse(c("GZMA", "GZMB", "XCL1"))
#'
#' @export

convert.human.to.mouse <- function(x){
  
  human = biomaRt::useMart("ensembl", 
                           dataset = "hsapiens_gene_ensembl")
  
  mouse = biomaRt::useMart("ensembl", 
                           dataset = "mmusculus_gene_ensembl")
  
  genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol"),
                            filters = "hgnc_symbol",
                            values = x , 
                            mart = human, 
                            attributesL = c("mgi_symbol"),
                            martL = mouse, 
                            uniqueRows=T)
  
  output.df <- genesV2 %>%
    dplyr::distinct(HGNC.symbol, .keep_all = TRUE) %>%
    dplyr::distinct(MGI.symbol, .keep_all = TRUE)
  
  print(head(output.df))
  
  return(output.df)
  
}



