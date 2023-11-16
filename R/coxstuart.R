#' Cox-Stuart test
#'
#' Perform Cox-Stuart test for location or dispersion.
#'
#' @param y input data.
#' @param type type of test. Can be:
#' \itemize{
#' \item \code{"trend"}: test for changes in trend.
#' \item \code{"deviation"}: test for changes in deviation.
#' \item \code{"dispersion"}: test for changes in dispersion (range).
#' }
#' @param alpha significance level.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{H}: hypothesis outcome.
#' \item \code{p.value}: corresponding p-value.
#' \item \code{Htxt}: textual description of the hypothesis outcome.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @keywords htest
#'
#' @examples
#' coxstuart(referrals)
#'
#' @export coxstuart

coxstuart <- function(y,type=c("trend","deviation","dispersion"),alpha=0.05){

  type <- match.arg(type,c("trend","deviation","dispersion"))

  switch(type,
         "trend" = {
           z <- y
         },
         "deviation"={
           z <- split(y)
           z <- apply(z,2,sd)
         },
         "dispersion"={
           z <- split(y)
           z <- apply(z,2,function(x){diff(range(x))})
         })

  # Find number of pairs for comparison
  n <- length(z)
  k <- ceiling(n/2)

  # Create pairs
  idx1 <- 1:(n-k)
  idx2 <- (1+k):n
  pair <- z[idx1] - z[idx2]

  # Calculate statistic
  stat <- min(c(sum(pair>0),sum(pair<0)))

  # P-value
  if (sum(pair)!=0){
    p <- pbinom(stat,k,0.5)
  } else {
    p <- 1
  }

  if (p <= alpha/2){
    H <- 1
    txt <- "H1: There is change (upwards or downwards)."
  } else {
    H <- 0
    txt <- "H0: There is no change present."
  }

  return(list(H=H,p.value=p,Htxt=txt))

}


split <- function(y){
# Helper function

  n <- length(y)

  # Find K according to Cox Stuart guidelines
  if (n < 48){
    k <- 2
  } else if(n < 64){
    k <- 3
  } else if(n < 90){
    k <- 4
  } else {
    k <- 5
  }

  # Split time series to subsamples
  rem <- n %% k
  kk <- floor(n/k)
  idx <- array(1:(n-rem),c(k,kk))
  idx[,(round(kk/2)+1):kk] <- idx[,(round(kk/2)+1):kk] + rem
  z <- array(y[idx], c(k,kk))
  return(z)

}
