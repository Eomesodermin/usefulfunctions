#' @title makeTransparent
#' @description Make a color transparent
#' @author Ricardo Oliveros-Ramos https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
#'
#' @param color a color e.g "red"
#' @param percent a number between 0 and 100 indicating percentage transparent
#' @param name an optional name for the color
#'
#' @return a new color with transparency
#' @examples
#' makeTransparent(color = "red", percent = 50, name = NULL)
#'
#' @export

makeTransparent <- function(color, percent = 50, name = NULL) {

  if(percent<0 | percent>100) stop("percent must be between 0 and 100")

  ## Get RGB values for named color
  rgb.val <- grDevices::col2rgb(color)

  # Make new color using input color as base and alpha set by percent variable
  new.col <- grDevices::rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                            max = 255,
                            alpha = (100 - percent) * 255 / 100,
                            names = name)

  # return new color
  return(new.col)

}




