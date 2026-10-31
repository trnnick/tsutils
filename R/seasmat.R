# Internal helpers shared by seasplot and seastest. Not exported.

# Detrend a series (if needed) and reorganise it into a seasonal matrix.
#
# Handles the common input processing of seasplot and seastest: identifying m and s,
# validating labels and cma, testing for trend (Cox-Stuart on the CMA), choosing the
# decomposition and removing the trend, and finally arranging the (detrended) series
# in a matrix with one row per season and one column per seasonal period. Incomplete
# first and last seasons are padded with NA.
#
# Returns a list with:
#   season: matrix of (detrended) seasonal elements (ns x m).
#   m, s: seasonal period and starting period used.
#   labels: labels of the seasonal periods (column names of season).
#   trend: CMA used for detrending, or NULL if no trend was removed.
#   trend.exist, trend.pval: results of the trend test, or NULL if not tested.
#   decomposition: "multiplicative", "additive" or "none" (no trend removed).
#' @noRd
seasmat <- function(y,m=NULL,s=NULL,trend=NULL,alpha=0.05,
                    decomposition=c("multiplicative","additive","auto"),
                    cma=NULL,labels=NULL){

  decomposition <- match.arg(decomposition,c("multiplicative","additive","auto"))
  is.ts <- inherits(y,"ts")

  # Get m (seasonality)
  if (is.null(m)){
    if (is.ts){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }
  if (length(m) != 1 || m < 1 || m %% 1 != 0){
    stop("m must be a single positive integer.")
  }

  # Get starting period in seasonality if available
  if (is.null(s)){
    if (is.ts){
      s <- start(y)[2]
      # Temporal aggregation can mess-up s, so override if needed
      if (is.na(s)){s <- 1}
    } else {
      s <- 1
    }
  }
  s <- ((s - 1) %% m) + 1

  # Make sure that labels argument is fine
  if (!is.null(labels)){
    if (length(labels) < m){
      labels <- NULL
      warning("Incorrect number of labels. It must be equal to m. Defaulting to labels=NULL.")
    } else {
      labels <- labels[1:m]
    }
  }
  if (is.null(labels)){
    labels <- paste(1:m)
  }

  n <- length(y)

  if (min(y,na.rm=TRUE) <= 0){
    decomposition <- "additive"
  }

  # If cma is given, it is removed without testing when trend=TRUE, otherwise the trend is tested on it
  if (!is.null(cma)){
    if (!isTRUE(trend)){
      trend <- NULL
    }
    if (n != length(cma)){
      stop("Length of series and cma do not match.")
    }
  }

  # Calculate CMA
  if ((is.null(trend) || isTRUE(trend)) && is.null(cma)){
    if (anyNA(y)){
      stop("y contains missing values, so the trend cannot be estimated. Use trend=FALSE, provide cma, or remove leading/trailing NAs with leadtrail.",call.=FALSE)
    }
    cma <- cmav(y=y,ma=m,fill=TRUE,outplot=FALSE)
  }

  # Test for changes in the CMA (trend)
  if (is.null(trend)){
    trend.pval <- coxstuart(cma)$p.value
    trend.exist <- trend.pval <= alpha/2
    trend <- trend.exist
  } else {
    trend.exist <- NULL
    trend.pval <- NULL
  }

  # Remove trend
  if (isTRUE(trend)){
    if (decomposition == "auto"){
      # mseastest needs a few seasons, fall back to additive if it cannot be run
      is.mult <- tryCatch(mseastest(y,m,cma=cma)$is.multiplicative,
                          error=function(e){
                            warning("Could not test for multiplicative seasonality (",conditionMessage(e),"). Using additive decomposition.",call.=FALSE)
                            FALSE})
      decomposition <- if (isTRUE(is.mult)) "multiplicative" else "additive"
    }
    if (decomposition == "multiplicative"){
      ynt <- y/cma
    } else {
      ynt <- y - cma
    }
  } else {
    decomposition <- "none"
    ynt <- y
    cma <- NULL
  }

  # Fill with NA start and end of season and arrange in matrix
  ks <- s - 1
  ke <- (m - ((n + ks) %% m)) %% m
  ynt <- c(rep(NA,times=ks),as.vector(ynt),rep(NA,times=ke))
  ns <- length(ynt)/m
  ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
  colnames(ynt) <- labels
  rownames(ynt) <- paste0("s",1:ns)

  return(list(season=ynt,m=m,s=s,labels=labels,trend=cma,
              trend.exist=trend.exist,trend.pval=trend.pval,decomposition=decomposition))

}

# Friedman test on a seasonal matrix (as produced by seasmat).
#
# Only complete seasons (rows without NA) are used. If fewer than two complete seasons
# are available, or m == 1, the test is not performed and NA is returned. If the
# elements are tied within every season (e.g. a constant series) the test statistic is
# undefined, in which case no seasonality is reported (p-value of 1).
#' @noRd
seasfriedman <- function(season,alpha=0.05){

  season <- season[complete.cases(season),,drop=FALSE]

  if (ncol(season) < 2 || nrow(season) < 2){
    return(list(is.seasonal=NA,statistic=NA,pvalue=NA,nseas=nrow(season)))
  }

  test <- tryCatch(friedman.test(season),error=function(e){NULL})
  if (is.null(test) || !is.finite(test$p.value)){
    statistic <- 0
    pvalue <- 1
  } else {
    statistic <- unname(test$statistic)
    pvalue <- test$p.value
  }

  return(list(is.seasonal=pvalue <= alpha,statistic=statistic,pvalue=pvalue,nseas=nrow(season)))

}
