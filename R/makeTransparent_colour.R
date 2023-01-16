#' @title Make a color transparent
#' @description Make a color transparent
#' @author Ricardo Oliveros-Ramos https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
#' @usage makeTransparent("red", "blue", alpha=0.8)
#' 
#' @param ... a color
#' @param alpha a number between 0 and 1
#' 
#' @return a color
#' @examples
#' makeTransparent("red", alpha=0.5)
#' makeTransparent(rgb(0,0,0,0.5))
#' @export

makeTransparent = function(..., alpha = 0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}

