#' Centred moving average
#'
#' Calculate the Centred Moving Average (CMA) for time series.
#'
#' @param y input time series. Can be \code{ts} or \code{msts} object.
#' @param ma length of centred moving average. If \code{y} is a \code{ts} object then the default is its frequency. If it is a \code{msts} object the default is the maximum frequency.
#' @param fill if \code{TRUE}, then fill first and last ma/2 observations using exponential smoothing.
#' @param outplot if \code{TRUE}, then output a plot of the time series and the moving average.
#' @param fast if \code{TRUE}, then only a limited set of models are evaluated for CMA extrapolation.
#'
#' @return Centred moving average. If y is a ts object, then cma has the same properties.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' Ord K., Fildes R., Kourentzes N. (2017) \href{http://kourentzes.com/forecasting/2017/10/16/new-forecasting-book-principles-of-business-forecasting-2e/}{Principles of Business Forecasting, 2e}. \emph{Wessex Press Publishing Co.}, p.109.
#'
#' @keywords cma ts
#'
#' @examples
#' cmav(referrals,outplot=TRUE)
#'
#' @importFrom forecast ets
#' @export cmav

cmav <- function(y,ma=NULL,fill=c(TRUE,FALSE),outplot=c(FALSE,TRUE),fast=c(TRUE,FALSE)){

  # Set defaults
  fill <- fill[1]
  outplot <- outplot[1]
  fast <- fast[1]

  # Get MA length
  if (is.null(ma)){
    if (any(class(y) == "ts") || any(class(y) == "msts")){
      ma <- stats::frequency(y)
    } else {
      stop("MA length not defined (y not ts or msts object).")
    }
  }

  # Get bounds for MA and correct length
  n <- length(y)
  mlbounds <- c(floor(ma/2)+1, n-floor(ma/2))
  isodd = ma %% 2
  if (isodd == 0){
    ml <- ma+1
  } else {
    ml <- ma
  }

  # Calculate MA
  # Loop across MA-order to speed up things
  mamat <- matrix(NA,nrow=ml,ncol=(n-ml+1))
  for (i in 1:ml){
    mamat[i,] <- y[i:(n-ml+i)]
  }
  if (isodd == 0){
    mamat[c(1,ml),] <- mamat[c(1,ml),]/2
  }
  mamat <- colSums(mamat)/ma

  # Carry over properties of y
  cma <- y
  cma[] <- NA
  cma[mlbounds[1]:mlbounds[2]] <- mamat

  # Fill MA is requested
  if (fill == TRUE){
    if (fast == FALSE){
      if ((n-mlbounds[2]) >= 1){
        cma[(mlbounds[2]+1):n] <- as.vector(forecast(ets(cma[(mlbounds[1]:mlbounds[2])],
                                                         model="ZZN"),h=(n-mlbounds[2]))$mean)
      }
      if ((mlbounds[1]-1) >= 1){
        cma[1:(mlbounds[1]-1)] <- rev(as.vector(forecast(ets(rev(cma[(mlbounds[1]:mlbounds[2])]),
                                                             model="ZZN"),h=(mlbounds[1]-1))$mean))
      }
    } else {
      fit <- ets(cma[(mlbounds[1]:mlbounds[2])],model="AZN")
      if ((n-mlbounds[2]) >= 1){
        cma[(mlbounds[2]+1):n] <- as.vector(forecast(fit,h=(n-mlbounds[2]))$mean)
      }
      if ((mlbounds[1]-1) >= 1){
        cma[1:(mlbounds[1]-1)] <- rev(as.vector(forecast(ets(rev(cma[mlbounds[1]:mlbounds[2]]),
                                                             fit,use.initial.values=FALSE),h=(mlbounds[1]-1))$mean))
      }
    }
  }

  if (outplot == TRUE){
    graphics::plot(y)
    graphics::lines(cma,col=RColorBrewer::brewer.pal(3,"Set1")[1])
  }

  return(cma)

}
