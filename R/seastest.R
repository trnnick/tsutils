#' Seasonality test
#'
#' Test a time series for deterministic seasonality using the Friedman test on the distributions of the (detrended) seasonal periods. The series can be automatically detrended as needed. This is the test reported by \code{\link{seasplot}}.
#'
#' @param y input time series. Can be \code{ts} object.
#' @param m seasonal period. If \code{y} is a \code{ts} object then the default is its frequency.
#' @param s starting period in the season. If \code{y} is a \code{ts} object then this is picked up from \code{y}.
#' @param trend if \code{TRUE}, then presence of trend is assumed and removed. If \code{FALSE} no trend is assumed. Use \code{NULL} to identify automatically.
#' @param alpha significance level for statistical tests.
#' @param decomposition type of seasonal decomposition. This can be \code{"multiplicative"}, \code{"additive"} or \code{"auto"}. If \code{y} contains non-positive values then this is forced to \code{"additive"}.
#' @param cma input precalculated level/trend for the analysis. If \code{trend=TRUE} it is removed from \code{y} without testing, otherwise the presence of trend is tested on \code{cma} (overriding \code{trend=FALSE}).
#'
#' @details
#' The series is first detrended (if needed) using a centred moving average (see \code{\link{cmav}}), with the presence of trend tested using the Cox-Stuart test (see \code{\link{coxstuart}}) when \code{trend=NULL}. The detrended series is then arranged in a matrix, with one row per season and one column per seasonal period, and the Friedman test is used to test whether the seasonal periods come from the same distribution. Only complete seasons are used for the test, so at least two complete seasons are needed.
#'
#' @return An object of class \code{seastest} containing:
#' \itemize{
#' \item \code{is.seasonal}: \code{TRUE}/\code{FALSE} result of the seasonality test. \code{NA} if the test could not be performed.
#' \item \code{statistic}: Friedman test statistic.
#' \item \code{pvalue}: p-value of the seasonality test.
#' \item \code{alpha}: significance level used.
#' \item \code{season}: matrix of (detrended) seasonal elements.
#' \item \code{trend}: CMA estimate (using \code{\link{cmav}}) or \code{NULL} if no trend was removed.
#' \item \code{trend.exist}: \code{TRUE}/\code{FALSE} results of trend test, or \code{NULL} if not tested.
#' \item \code{trend.pval}: p-value of trend test (Cox-Stuart), or \code{NULL} if not tested.
#' \item \code{decomposition}: type of decomposition used. \code{"none"} if no trend was removed.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{seasplot}}, \code{\link{mseastest}}, \code{\link{trendtest}}.
#'
#' @keywords ts htest
#'
#' @examples
#' seastest(referrals)
#' seastest(referrals,trend=TRUE,decomposition="auto")$pvalue
#'
#' @export seastest

seastest <- function(y,m=NULL,s=NULL,trend=NULL,alpha=0.05,
                     decomposition=c("multiplicative","additive","auto"),cma=NULL){

  sm <- seasmat(y=y,m=m,s=s,trend=trend,alpha=alpha,decomposition=decomposition,cma=cma)
  test <- seasfriedman(sm$season,alpha=alpha)

  if (is.na(test$pvalue)){
    warning("Seasonality not tested. At least two complete seasons and m > 1 are needed.")
  }

  out <- list(is.seasonal=test$is.seasonal,statistic=test$statistic,pvalue=test$pvalue,
              alpha=alpha,season=sm$season,trend=sm$trend,trend.exist=sm$trend.exist,
              trend.pval=sm$trend.pval,decomposition=sm$decomposition)
  out <- structure(out,class="seastest")
  return(out)

}

#' @export
#' @method summary seastest
summary.seastest <- function(object,...){
  print(object)
}

#' @export
#' @method print seastest
print.seastest <- function(x,...){
  cat("Seasonality test (Friedman) on the seasonal periods\n")
  if (x$decomposition == "none"){
    cat("Series not detrended.\n")
  } else {
    cat(paste0("Series detrended (",x$decomposition,").\n"))
  }
  if (!is.null(x$trend.exist)){
    cat(paste0("Evidence of trend: ",x$trend.exist,"  (pval: ",round(x$trend.pval,3),")\n"))
  }
  if (is.na(x$pvalue)){
    cat("Seasonality not tested (not enough complete seasons).\n")
  } else {
    cat(paste0("Evidence of seasonality: ",x$is.seasonal,"  (statistic: ",round(x$statistic,3),
               ", pval: ",round(x$pvalue,3),", alpha: ",x$alpha,")\n"))
  }
  invisible(x)
}
