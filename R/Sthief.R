#' Temporal hierarchy S matrix
#'
#' Calculate the temporal hierarchy summing matrix S for a given time series of seasonal periodicity.
#'
#' @param y input time series (a \code{ts} object) or an integer.
#' @param k optional vector of temporal aggregation levels to include. These must be factors of the seasonal periodicity, i.e., \code{frequency(y)} or \code{y} if it is an integer. The bottom level (\code{k=1}) is always included. If \code{NULL}, all feasible aggregation levels are used.
#'
#' @return S matrix. Rows are ordered from the top level (largest \code{k}) to the bottom level (\code{k=1}). The top level does not need to be a single node, for example \code{Sthief(24,k=c(12,1))} gives two annual nodes over 24 months.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' Athanasopoulos, G., Hyndman, R. J., Kourentzes, N., & Petropoulos, F. (2017). \href{https://kourentzes.com/forecasting/2017/02/27/forecasting-with-temporal-hierarchies-3/}{Forecasting with temporal hierarchies}. European Journal of Operational Research, 262(1), 60-74.
#'
#' @seealso \code{\link{plotThief}}.
#'
#' @keywords ts
#'
#' @examples
#' Sthief(AirPassengers)
#'
#' # Only some aggregation levels: monthly, quarterly, and annual
#' Sthief(12,k=c(12,3))
#'
#' # Two years of monthly data, with half-yearly, annual and 2-yearly levels
#' S <- Sthief(24,k=c(24,12,6))
#' plotThief(S,freq="monthly")
#'
#' @export Sthief


Sthief <- function(y,k=NULL){

  # Get time series frequency
  if (any(class(y) == "ts")){
    f <- frequency(y)
  } else if(is.numeric(y) && length(y)==1 && (y %% 1 == 0)){
    f <- y
  } else {
    stop("Argument y must be either a ts object or an integer.")
  }
  if (f < 1 || f %% 1 != 0){
    stop("The seasonal periodicity must be a positive integer.")
  }

  # Find feasible aggregation levels
  fAll <- f/(1:f)
  fAll <- fAll[fAll %% 1 == 0]

  # Keep only the requested ones, plus the bottom level
  if (!is.null(k)){
    if (!is.numeric(k) || any(!(k %in% fAll))){
      stop(paste0("Argument k must contain factors of the seasonal periodicity (",f,"): ",
                  paste(sort(fAll),collapse=", "),"."))
    }
    fAll <- sort(unique(c(k,1)),decreasing=TRUE)
  }
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
