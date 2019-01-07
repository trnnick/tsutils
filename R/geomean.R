#' Geometric mean
#'
#' Calculate the geometric mean.
#'
#' @param x input data (will be considered as a vector).
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The geometirc mean of the values in \code{x}.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @keywords geomean
#'
#' @examples
#' geomean(c(0.5,1,1.5))
#'
#' @export geomean

geomean <- function(x,na.rm=c(FALSE,TRUE),...){

    na.rm <- na.rm[1]
    gm <- exp(mean(log(x),na.rm=na.rm),...)
    return(gm)

}
