#' @title Get Batlow
#' @description This function returns a vector of colors from the batlow palette
#'
#' @return a vector of colors from the batlow palette
#' @examples
#' # requires scico package
#' # install.packages("scico")
#' # library("scico")
#' batlow.pal <- Get.batlow()
#' @export

Get.batlow <- function(){

  batlow.pal <- scico(100, palette = 'batlow')

}

