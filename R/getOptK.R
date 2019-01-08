#' Optimal temporal aggregation level for AR(1), MA(1), ARMA(1,1)
#'
#' Calculate the theoretically optimal temporal aggregation level for AR(1), MA(1) and ARMA(1,1) time series.
#'
#' @param y a time series that must be of either \code{ts} or \code{msts} class.
#' @param m maximum aggregation level.
#' @param type type of data generating process. Can be:
#' \itemize{
#' \item{\code{"ar"}}{: For AR(1) series.}
#' \item{\code{"ma"}}{: For MA(1) series.}
#' \item{\code{"arma"}}{: For ARMA(1,1) series.}
#' }
#'
#' @return Identified optimal temporal aggregation level.
#'
#' @references
#' \itemize{
#' \item{Kourentzes, N., Rostami-Tabar, B., & Barrow, D. K. (2017). \href{http://kourentzes.com/forecasting/2017/04/22/demand-forecasting-by-temporal-aggregation-using-optimal-or-multiple-aggregation-levels/}{Demand forecasting by temporal aggregation: using optimal or multiple aggregation levels?}. Journal of Business Research, 78, 1-9.}
#' \item{Rostami-Tabar, B., Babai, M. Z., Syntetos, A., & Ducq, Y. (2013). Demand forecasting by temporal aggregation. Naval Research Logistics (NRL), 60(6), 479-498.}
#' \item{Rostami-Tabar, B., Babai, M. Z., Syntetos, A., & Ducq, Y. (2014). A note on the forecast performance of temporal aggregation. Naval Research Logistics (NRL), 61(7), 489-500.}
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @keywords htest
#'
#' @examples
#' getOptK(referrals)
#'
#' @export getOptK

getOptK <- function(y,m=12,type=c("ar","ma","arma")){

  # Default
  type <- match.arg(type,c("ar","ma","arma"))
  sigma <- 1

  if (!(any(class(y) == "ts") | any(class(y) == "msts"))){
    stop("Input y must be of class ts or msts.")
  }

  # Get DGP parameter
  if (type == "ma"){
    fit <- tryCatch({
      Arima(y,order=c(0,0,1),method="ML")
    }, error = function(e) {
      Arima(y,order=c(0,0,1))
    })
    t <- fit$coef[1]
    p <- NULL
  } else if (type == "ar"){
    fit <- tryCatch({
      Arima(y,order=c(1,0,0),method="ML")
    }, error = function(e) {
      Arima(y,order=c(1,0,0))
    })
    p <- fit$coef[1]
    t <- NULL
  } else {
    fit <- tryCatch({
      Arima(y,order=c(1,0,1),method="ML")
    }, error = function(e) {
      Arima(y,order=c(1,0,1))
    })
    p <- fit$coef[1]
    t <- fit$coef[2]
  }

  # Aggregate series
  Y <- MAPA::tsaggr(y,1:m)

  # Calculate MSE
  mse <- vector("numeric",m)
  for (i in 1:m){

    y.a <- Y$out[[i]]
    fit.ses <- ets(y.a,model="ANN")
    alpha <- fit.ses$par[1]

    # Get aggregate MSE
    if (type == "ma"){
      # MA
      mse[i] <- ((2*m*(1+t^2))/(2-alpha))*(sigma^2)
    } else if (type == "ar"){
      # AR
      s1.ar <- 0
      s2.ar <- 0
      s3.ar <- 0
      for (k in 1:m) {
        if (k <= m-1){
          s1.ar=s1.ar+(2*(m-k)*(p^(k-1)))
        }
        if (k >= 2){
          s3.ar=s3.ar+((k-1)*(p^(2*m-k)))
        }
        s2.ar=s2.ar+(k*(p^(k-1)))
      }
      part1.ar <- (2*(m+(s1.ar*p)))/((2-alpha)*(1-(p^2)))
      p.prime <- p^m
      part2.ar <- (((2*alpha)*(s2.ar+s3.ar))/((2-alpha)*(1-p.prime+alpha*p.prime)))*(p/(1-(p^2)))
      mse[i] <- (part1.ar-part2.ar)*(sigma^2)
    } else {
      p.prime <- p^m
      # ARMA
      part1.arma11 <- (m*(1-2*p*t+(t^2)))/(1-(p^2))
      s1.arma11 <- 0
      s2.arma11 <- 0
      s3.arma11 <- 0
      for (k in 1:m-1) {
        s2.arma11=s2.arma11+(k*(p^(k-1)))
        if (k <= m-1){
          s1.arma11=s1.arma11+(2*(m-k)*(p^(k-1)))
        }
        if (k >= 2){
          s3.arma11=s3.arma11+((k-1)*(p^(2*m-k)))
        }
      }
      part2.arma11 <- s1.arma11*(((p-t)*(1-p*t))/(1-(p^2)))
      part3.arma11 <- 2*((part1.arma11+part2.arma11)/(2-alpha))
      mse[i] <- (part3.arma11-((((2*alpha)*(s2.arma11+s3.arma11))/((2-alpha)*(1-p.prime+alpha*p.prime)))*(((p-t)*(1-p*t))/(1-p^2))))*(sigma^2)
    }

  }

  k <- which(mse == min(mse))[1]
  return(k)

}
