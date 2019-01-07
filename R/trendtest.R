#' Test a time series for trend
#'
#' Test a time series for trend by either fitting exponential smoothing models and comparing then using the AICc, or by using the non-parametric Cox-Stuart test. The tests can be augmented by using multiple temporal aggregation.
#'
#' @param y a time series that must be of either \code{ts} or \code{msts} class.
#' @param extract if \code{TRUE} then the centred moving average of the time series is calculated and the test is performed on that. Otherwise, the test is performed on the raw data.
#' @param type type of test. Can be:
#' \itemize{
#' \item{\code{"aicc"}}{: test by comparing the AICc of exponential smoothing models. See details.}
#' \item{\code{"cs"}}{: test by using the Cox-Stuart test. See details.}
#' }
#' @param mta If \code{TRUE} augment testing by using Multiple Temporal Aggregation.
#'
#' @return The function returns \code{TRUE} when there is evidence of trend and \code{FALSE} otherwise.
#'
#' @details All tests are performed at 5% significance level. When exponential smoothing is used the following three models are compared: ETS(A,N,N), ETS(A,A,N), ETS(A,Ad,N).
#'
#' @references The multiple temporal aggregation follows the construction approach suggested by Kourentzes, N., Petropoulos, F., & Trapero, J. R. (2014). \href{http://kourentzes.com/forecasting/2014/04/19/improving-forecasting-by-estimating-time-series-structural-components-across-multiple-frequencies/}{Improving forecasting by estimating time series structural components across multiple frequencies}. International Journal of Forecasting, 30(2), 291-302.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @keywords htest
#'
#' @examples
#' trendtest(referrals,TRUE)
#'
#' @export trendtest

trendtest <- function(y,extract=c("FALSE","TRUE"),type=c("aicc","cs"),mta=c(FALSE,TRUE)){

    # Defaults
    extract <- extract[1]
    type <- match.arg(type,c("aicc","cs"))
    mta <- mta[1]

    if (!(any(class(y) == "ts") | any(class(y) == "msts"))){
        stop("Input y must be of class ts or msts.")
    }

    # Extract trend
    if (extract == TRUE){
        y.in <- cmav(y)
    } else {
        y.in <- y
    }

    # MTA if needed
    if (mta == TRUE){
        f <- frequency(y)
        yaggr <- MAPA::tsaggr(y.in,1:f)$out
        # Remove levels with very few observations, but always keep dissagregate
        idx <- unlist(lapply(yaggr,function(x){length(x)>=5}))
        idx[1] <- TRUE
        yaggr <- yaggr[idx]
    } else {
        yaggr <- MAPA::tsaggr(y.in,1)$out
    }

    if (type == "cs"){
        H <- unlist(lapply(yaggr,test.cs))
    } else {
        H <- unlist(lapply(yaggr,test.aicc))
    }

    H <- median(H)>0.5

    return(H)

}

test.cs <- function(y){
    H <- coxstuart(y)$H
    return(H)
}

test.aicc <- function(y){
    # Construct and AICc test
    aicc1 <- tryCatch({ets(y,model="ANN")$aicc}, error = function(e){10^5})
    aicc2 <- tryCatch({ets(y,model="AAN",damped=FALSE)$aicc}, error = function(e){10^5})
    aicc3 <- tryCatch({ets(y,model="AAN",damped=TRUE)$aicc}, error = function(e){10^5})
    aicc <- c(aicc1,aicc2,aicc3)
    H.aicc <- which(aicc == min(aicc))[1] != 1
    return(H.aicc)
}
