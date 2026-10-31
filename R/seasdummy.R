#' Create seasonal dummy variables.
#'
#' Create binary or trigonometric seasonal dummies.
#'
#' @param n number of observations to create.
#' @param m seasonal periodicity. If \code{NULL} it will take the information from the provided time series (\code{y} argument). See notes.
#' @param y this is an optional time series input that can be used to get seasonal periodicity (\code{m}) and the start point.
#' @param type type of seasonal dummies to create.
#' \itemize{
#' \item \code{"bin}: binary dummies.
#' \item \code{"trg}: trigonometric dummies. See notes.
#' }
#' @param full If full is \code{TRUE}, then keeps the m-th dummy that is co-linear to the rest. See notes.
#'
#' @note If the seasonal periodicity is fractional then the the type will be overriden to trigonometric and only two seasonal dummies with be produced. One cosine and one sine.
#'
#' @return An array with seasonal dummies, where rows correspond observations and columns to dummy variables.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
#'
#' @examples
#' seasdummy(24,12)
#'
#' @export seasdummy

seasdummy <- function(n,m=NULL,y=NULL,type=c("bin","trg"),full=c(FALSE,TRUE)){
  # Create deterministic seasonality dummies

  # Default
  type <- match.arg(type,c("bin","trg"))
  full <- full[1]

  # Get time series information if available
  if (!is.null(y)){
    if (is.null(m)){
      m <- frequency(y)
    }
    start <- start(y)
    # Deal with fractional seaosnalities
    isdd <-length(start)==2
    if (isdd){
      start <- start[2]
    } else {
      start <- start %% 1
    }
  } else {
    start <- 1
    isdd <- TRUE
  }

  if (is.null(m)){
    stop("Seasonality not specified.")
  }

  # Create dummies
  if ((m %% 1) == 0){
    if (type == "bin"){
      # Seasonal period of each observation, starting from period start
      p <- ((start - 1) + seq_len(n) - 1) %% m + 1
      x <- diag(1,m)[p,,drop=FALSE]
    } else { # trg
      x <- array(0,c(n,m))
      # Shift time for start, so that every harmonic has the correct phase
      t <- 1:n + (start - isdd)
      for (i in 1:(m/2)){
        x[,1+(i-1)*2] <- cos((2*t*pi*i)/m)
        x[,2+(i-1)*2] <- sin((2*t*pi*i)/m)
      }
    }
    # Trim co-linear dummy
    if (full==FALSE){
      x <- x[,1:(m-1),drop=FALSE]
    }
  } else {
    # Fractional seasonality
    x <- array(0,c(n,2))
    t <- 1:n
    x[,1] <- cos((2*t*pi)/m+(2*pi*(start-isdd))) # Added phase shift for start
    x[,2] <- sin((2*t*pi)/m+(2*pi*(start-isdd)))
  }

  return(x)

}
