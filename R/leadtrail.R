#' Remove leading/training zeros/NAs
#'
#' Remove leading or trailing zeros or NAs from a vector.
#'
#' @param x vector of values to check.
#' @param rm what to remove, can be \code{"zeros"} or \code{"na"}.
#' @param lead If \code{TRUE}, then leading values are removed.
#' @param trail If \code{TRUE}, then trailing values are removed.
#'
#' @return Resulting vector.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @examples
#' x <- c(rep(0,5),rnorm(100),rep(0,5))
#' leadtrail(x)
#'
#' @export leadtrail

leadtrail <- function(x,rm=c("zeros","na"),lead=c(TRUE,FALSE),trail=c(TRUE,FALSE)){

  # Defaults
  rm <- match.arg(rm,c("zeros","na"))
  lead <- lead[1]
  trail <- trail[1]

  # Select what to remove
  if (rm=="zeros"){
    idx <- which(x == 0)
  } else {
    idx <- which(is.na(x))
  }

  n <- length(x)
  l <- length(idx)

  # Handle leading observations
  if (lead==TRUE & l>0){

    if (idx[1]==1){
      d.idx <- diff(idx)
      loc <- which(d.idx > 1)[1]
      if (is.na(loc)){
        loc <- l
      }
      lead.rm <- 1:loc
    } else {
      lead.rm <- NULL
    }

  } else {
    lead.rm <- NULL
  }

  # Handle trailing observations
  if (trail==TRUE & l>0){

    if (tail(idx,1)==n){
      d.idx <- diff(rev(idx))
      loc <- which(d.idx != -1)[1]
      if (is.na(loc)){
        loc <- l
      }
      trail.rm <- (n-loc+1):n
    } else {
      trail.rm <- NULL
    }

  } else {
    trail.rm <- NULL
  }

  keep <- rep(TRUE,n)
  keep[lead.rm] <- FALSE
  keep[trail.rm] <- FALSE

  y <- x[keep]
  return(y)

}

