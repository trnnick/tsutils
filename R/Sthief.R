#' Temporal hierarchy S matrix
#'
#' Calculate the temporal hierarchy summing matrix S for a given time series of seasonal periodicity.
#'
#' @param y input time series (a \code{ts} object) or an integer. 
#' 
#' @return S matrix.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' Athanasopoulos, G., Hyndman, R. J., Kourentzes, N., & Petropoulos, F. (2017). \href{http://kourentzes.com/forecasting/2017/02/27/forecasting-with-temporal-hierarchies-3/}{Forecasting with temporal hierarchies}. European Journal of Operational Research, 262(1), 60-74.
#'
#' @keywords ts
#'
#' @examples
#' Sthief(AirPassengers)
#'
#' @export Sthief


Sthief <- function(y){

  # Get time series frequency
  if (any(class(y) == "ts")){
    f <- frequency(y)    
  } else if(is.numeric(y) && length(y)==1 && (y %% 1 == 0)){
    f <- y
  } else {
    stop("Argument y must be either a ts object or an integer.")
  }


  # Find feasible aggregation levels
  fAll <- f/(1:f)
  fAll <- fAll[fAll %% 1 == 0]
  kmax <- length(fAll)

  # Construct S matrix
  S.thief <- vector("list",kmax)
  for (k in 1:kmax){
    S.thief[[k]] = matrix(0,f/fAll[k],f)
    for (i in 1:(f/fAll[k])){
      S.thief[[k]][i,(1+fAll[k]*(i-1)):(fAll[k]+fAll[k]*(i-1))] <- 1
    }
  }
  S.thief <- do.call("rbind",S.thief)

  return(S.thief)

}

