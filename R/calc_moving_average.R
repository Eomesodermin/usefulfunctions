
#' @title Calculate moving average
#' @description 
#' function to calculate the moving average, was adapted from Riesenberg, S. et al.Nat. Commun. 2015, doi: 10.1038/ncomms9755
#' See also http://www.cookbook-r.com/ free usage under creative commons licence http://creativecommons.org/licenses/by-sa/3.0/legalcode
#' 
#' 
#' @param x A numeric vector
#' @param n The number of values to average over
#' @param centered Whether to center the average
#' @return A vector of the same length as x
#' @examples
#' calc.moving.average(c(1,2,3,4,5), n=3)
#' @export


calc.moving.average <- function(x, n = 1, centered = FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

