#' Theta method
#'
#' Estimate Theta method.
#'
#' @param y input time series. Can be \code{ts} object.
#' @param m seasonal period. If \code{y} is a \code{ts} object then the default is its frequency.
#' @param sign.level significance level for trend and seasonality tests.
#' @param cost0 cost function of theta0 line. Can be:
#' \itemize{
#' \item \code{"MSE"}: mean squared error.
#' \item \code{"MdSE"}: median squared error.
#' \item \code{"MAE"}: mean absolute error.
#' \item \code{"MdAE"}: median absolute error.
#' }
#' @param cost2 cost function of theta2 line. Same options as \code{cost0}.
#' @param costs cost function of seasonal element. Same options as \code{cost0}.
#' @param season type of seasonal decomposition. This can be \code{"multiplicative"}, \code{"additive"} or \code{"auto"}. If \code{y} contains non-positive values then this is forced to \code{"additive"}. With \code{"auto"} the type is chosen with \code{\link{mseastest}}. Whether seasonality is present at all is always tested (see details).
#' @param cma input precalculated level/trend for the analysis. Use \code{NULL} to estimate internally.
#' @param outliers provide vector of location of observations that are considered outliers (see \code{\link{residout}}). These will be considered in the estimation of theta0. For no outliers use \code{NULL}.
#' @param multiplicative deprecated, use \code{season} instead. Accepts the same values as \code{season}. It is kept for backwards compatibility, gives a warning when used, and will be removed in a future version.
#'
#' @return An object of class \code{theta}, containing:
#' \itemize{
#' \item \code{"method"}: "Theta".
#' \item \code{"y"}: the input time series.
#' \item \code{"m"}: seasonal periods.
#' \item \code{"exist"}: Statistical testing results, \code{exist[1]} is the result for trend, \code{exist[2]} is for season.
#' \item \code{"multiplicative"}: If \code{TRUE}, then seasonality is modelled multiplicatively.
#' \item \code{"theta0"}: fitted theta0 line values.
#' \item \code{"theta2"}: fitted theta2 line values.
#' \item \code{"season"}: fitted season values.
#' \item \code{"x.out"}: modelled outliers.
#' \item \code{"cost"}: cost functions for theta0, theta2 and season components.
#' \item \code{"a"}: SES parameters of theta2.
#' \item \code{"b"}: regression parameters of theta0.
#' \item \code{"p"}: coefficients of outliers from theta0 and theta2 estimation.
#' \item \code{"g"}: pure seasonal exponential smoothing parameters.
#' \item \code{"fitted"}: fitted values.
#' \item \code{"residuals"}: in-sample residuals.
#' \item \code{"MSE"}: in-sample Mean Squared Error.
#' \item \code{"outliers"}: location of outliers provided in \code{outliers}.
#' }
#'
#' @details This implementation of the Theta method tests automatically for seasonality and trend. Seasonal decomposition can be done either additively or multiplicatively and the seasonality is treated as a pure seasonal model. The various Theta components can be optimised using different cost functions. The originally proposed Theta method always assumed multiplicative seasonality and presence of trend, while all theta lines were optimised using MSE and seasonality was estimated using classical decomposition.
#'
#' The argument \code{season} replaces the argument \code{multiplicative} of earlier versions. Calls that use \code{multiplicative} still work, but give a deprecation warning. The output element \code{multiplicative} is unchanged.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' \itemize{
#' \item The original Theta method was proposed by: Assimakopoulos, V., & Nikolopoulos, K. (2000). The theta model: a decomposition approach to forecasting. International journal of forecasting, 16(4), 521-530. See details in how the implementation here differs.
#' \item The THieF forecasting methodology used for \code{theta.thief} is proposed by: Athanasopoulos, G., Hyndman, R. J., Kourentzes, N., & Petropoulos, F. (2017). Forecasting with temporal hierarchies. European Journal of Operational Research, 262(1), 60-74.
#' }
#'
#' @keywords ts
#'
#' @examples
#' fit <- theta(referrals)
#' plot(fit)
#'
#' @export theta

theta <- function(y,m=NULL,sign.level=0.05,
                  cost0=c("MSE","MdSE","MAE","MdAE"),
                  cost2=c("MSE","MdSE","MAE","MdAE"),
                  costs=c("MSE","MdSE","MAE","MdAE"),
                  season=c("multiplicative","additive","auto"),cma=NULL,
                  outliers=NULL,multiplicative=NULL){

  # Defaults
  cost0 <- match.arg(cost0,c("MSE","MdSE","MAE","MdAE"))
  cost2 <- match.arg(cost2,c("MSE","MdSE","MAE","MdAE"))
  costs <- match.arg(costs,c("MSE","MdSE","MAE","MdAE"))

  # Deprecated argument multiplicative, replaced by season
  if (!is.null(multiplicative)){
    if (missing(season)){
      warning("Argument multiplicative is deprecated, use season instead.",call.=FALSE)
      season <- multiplicative
    } else {
      warning("Argument multiplicative is deprecated and ignored, as season is provided.",call.=FALSE)
    }
  }
  season <- match.arg(season,c("multiplicative","additive","auto"))

  n <- length(y)

  # Detrend and organise in seasonal matrix. This also identifies m, calculates the CMA if
  # not provided and chooses the type of seasonality (forced to additive for non-positive data)
  sm <- seasmat(y=y,m=m,trend=TRUE,decomposition=season,cma=cma)
  m <- sm$m
  s <- sm$s
  cma <- sm$trend
  ynt <- sm$season
  multiplicative <- sm$decomposition == "multiplicative"

  # Test for trend
  trend.exist <- coxstuart(cma)$p.value <= sign.level/2

  # Test for seasonality, not seasonal if it cannot be tested
  season.exist <- isTRUE(seasfriedman(ynt,alpha=sign.level)$is.seasonal)

  # If seasonality exist then decompose
  if (season.exist == TRUE){
    # y.des <- y/rep(ynt,ceiling(n/m)+1)[1:n]
    y.des <- cma
  } else {
    y.des <- y
  }

  # Create theta lines

  # Include in theta0 outliers if any provided
  if (is.null(outliers)){
    X.out <- NULL
    n.out <- 0
  } else {
    n.out <- length(outliers)
    #     X.out <- array(0,c(n,n.out))
    #     X.out[outliers+n*(0:(n.out-1))] <- 1
    # Through the MA the outlier is spread across m observations, create dummies to account for that
    m.half <- floor((m + (m+1) %% 2)/2)
    X.out <- array(0,c(n,n.out))
    for (i in 1:n.out){
      X.out[max(1,outliers[i]-m.half):min(n,outliers[i]+m.half),i] <- 1
    }
  }
  # Include trend component in theta0
  if (trend.exist == TRUE){
    X.trend <- matrix(c(1:n), nrow=n, ncol=1)
  } else {
    X.trend <- NULL
  }
  X <- cbind(matrix(1,nrow=n,ncol=1),X.trend,X.out)

  #   if (trend.exist == TRUE || !is.null(outliers)){
  if ((n.out + trend.exist)>=1){
    b0 <- solve(t(X)%*%X)%*%t(X)%*%y.des # Initialise theta0 parameters
    b <- opt.trnd(y.des,X,cost0,b0)      # Optimise theta0
  } else {
    # If no trend then theta0 is just a mean
    b <- mean(y.des)
  }

  # Create theta0 line without outliers
  # 0*y.des To take ts object properties
  theta0 <- X[,1:(1+trend.exist)]%*%matrix(b[1:(1+trend.exist)],ncol=1) + 0*y.des

  # Estimate SES
  theta2 <- 2*y.des - theta0              # Construct theta2
  a0 <- rbind(0.1,theta2[1],rep(0,n.out)) # Initialise theta2 parameters
  a <- opt.ses(theta2,cost2,a0,2,X.out)   # Optimise theta2

  # In-sample fit
  in.theta0 <- theta0
  in.theta2 <- fun.ses(theta2,a,X.out)$ins
  if (!is.null(X.out)){
    # Remove outlier from fit - the complete outlier will be modelled afterwards
    in.theta2 <- in.theta2 - X.out %*% matrix(a[3:(2+n.out)])
  }
  in.fit <- (in.theta0 + in.theta2)/2

  # Separate theta0 and theta2 parameters into a, b and outliers (p[,1:2])
  if (n.out > 0){
    p <- matrix(b[(2+trend.exist):(1+trend.exist+n.out)],nrow=n.out)
    p <- cbind(p,a[3:(2+n.out)])
    colnames(p) <- c('Theta0','Theta2')
  } else {
    p <- NULL
  }
  if (trend.exist == FALSE){
    b <- rbind(b[1], 0)
  } else {
    b <- b[1:2]
  }
  b <- matrix(b,ncol=1)
  a <- matrix(a[1:2],ncol=1)

  # Reseasonalise
  if (season.exist == TRUE){
    # Seasonality is modelled with a pure seasonal smoothing
    sout <- opt.sfit(ynt,costs,n,m,s,y,in.fit,multiplicative,outliers)
    season <- sout$season
    g <- sout$g
    if (n.out > 0){
      p <- cbind(p,matrix(sout$p,ncol=1,dimnames=list(NULL,'Season')))
    }
    # sout$in.season includes the outlier
    if (multiplicative == TRUE){
      in.fit <- in.fit * sout$in.season
    } else {
      in.fit <- in.fit + sout$in.season
    }
  } else {
    g <- NULL
    season <- NULL
  }

  # Prepare output
  exist <- rbind(trend.exist,season.exist)
  rownames(exist) <- c("Trend","Season")
  rownames(a) <- c("Alpha","Initial level")
  rownames(b) <- c("Intercept","Slope")
  if (season.exist==TRUE){
    g <- matrix(g,nrow=1,ncol=m+1)
    colnames(g) <- c("Gamma",paste("s",1:m,sep=""))
  }
  costf <- rbind(cost0,cost2,costs)
  rownames(costf) <- c("Theta0","Theta2","Seasonal")
  if (any(class(y) == "ts")){
    in.fit <- ts(in.fit,frequency=m,end=end(y))
    theta0 <- ts(theta0,frequency=m,end=end(y))
    theta2 <- ts(theta2,frequency=m,end=end(y))
  }

  return(structure(list("method"="Theta","y"=y,"m"=m,"exist"=exist,
                        "multiplicative"=multiplicative,
                        "theta0"=theta0,"theta2"=theta2,
                        "season"=season,"x.out"=X.out,
                        "cost"=costf,"a"=a,"b"=b,"p"=p,"g"=g,
                        "fitted"=in.fit,"residuals"=y-in.fit,
                        "MSE"=mean((y-in.fit)^2),"outliers"=outliers)
                   ,class="theta"))

}

#' @rdname theta
#'
#' @title forecast.theta
#' @description Forecast with fitted Theta method.
#'
#' @param object object of class \code{theta}.
#' @param h forecast horizon. If \code{h} is \code{NULL}, then the horizon is set equal to the the seasonal frequency.
#' @param ... additional arguments passed to functions.
#'
#' @importFrom forecast forecast
#'
#' @examples
#' forecast.theta(fit,h=12) # Or simply use forecast(fit)
#'
#' @export
#' @export forecast.theta

forecast.theta <- function(object,h=NULL,...){
  # Produce forecasts with Theta

  a <- object$a
  b <- object$b
  p <- object$p
  m <- object$m
  n <- length(object$y)
  X.out <- object$x.out
  theta2 <- object$theta2
  y <- object$y

  if (is.null(h)){
    h <- m
  }

  # Forecast theta line
  frc.theta0 <- b[1] + b[2]*((n+1):(n+h))
  if (!is.null(X.out)){
    n.out <- ncol(X.out)
    a.frc <- rbind(a,array(p[,2],c(n.out,1)))
  } else {
    a.frc <- a
  }
  frc.theta2 <- fun.ses(theta2,a.frc,X.out)$outs * rep(1,h)
  frc <- (frc.theta0 + frc.theta2)/2

  # Include seasonality
  if (object$exist[2] == TRUE){
    season <- rep(object$season, h %/% m + 1)[1:h]
    if (object$multiplicative == TRUE){
      frc <- frc * season
    } else {
      frc <- frc + season
    }
  } else {
    season <- NULL
  }

  # Convert to ts
  if (any(class(y) == "ts")){
    fstart <- tsp(y)[2] + deltat(y)
    frc <- ts(frc,start=fstart,frequency=m)
    frc.theta0 <- ts(frc.theta0,start=fstart,frequency=m)
    frc.theta2 <- ts(frc.theta2,start=fstart,frequency=m)
  }

  return(structure(list("method"=class(object),"mean"=frc,
                        "frc.theta0"=frc.theta0,"frc.theta2"=frc.theta2,
                        "frc.season"=season,"x"=y,
                        "fitted"=object$fitted,"residuals"=object$residuals),
                   class="forecast"))

}

#' @rdname theta
#'
#' @title plot.theta
#' @description Produce a plot of the fitted Theta method.
#'
#' @param x object of class \code{theta}.
#' @param thetalines if \code{TRUE}, then theta lines are included in the plot.
#' @param legend if \code{TRUE} a legend is added to the plot. Use \code{FALSE} to omit it, or a keyword such as \code{"bottomleft"} to change its position.
#'
#' @details For \code{plot.theta}, if \code{h > 0} the forecasts are added to the plot, together with the forecasts of the theta lines if \code{thetalines=TRUE}. Outliers provided in \code{outliers} are marked on the data. Arguments in \code{...} control the plot: \code{col} sets the colours of the data, the fit and forecast, theta0, theta2 and the outliers, in this order (recycled as needed); \code{lwd} and \code{lty} the lines; everything else, such as \code{main}, \code{xlab}, \code{ylab}, \code{xlim} or \code{ylim}, is passed to \code{plot}. The default colours are from the \code{Set1} palette of \code{\link[RColorBrewer]{brewer.pal}}.
#'
#' @examples
#' plot(fit,h=24)
#' plot(fit,thetalines=FALSE,h=12,main="Referrals",legend="bottomright")
#'
#' @export
#' @method plot theta

plot.theta <- function(x,thetalines=c(TRUE,FALSE),h=0,legend=TRUE,...){
  # Produce in-sample fit plot, with forecasts if h > 0

  thetalines <- isTRUE(as.logical(thetalines[1]))
  if (is.null(h)){
    h <- x$m
  }
  n <- length(x$y)
  y <- as.vector(x$y)

  # Time index, including forecast horizon
  if (inherits(x$y,"ts")){
    tt <- as.vector(time(x$y))
    dt <- deltat(x$y)
    xlab <- "Time"
  } else {
    tt <- 1:n
    dt <- 1
    xlab <- "Period"
  }
  tf <- tt[n] + seq_len(h)*dt
  if (h > 0){
    frc <- forecast.theta(x,h=h)
  }

  # Split user arguments between the plotted lines and the plot frame
  args <- list(...)
  anames <- names(args)
  if (is.null(anames)){anames <- rep("",length(args))}
  el.names <- c("col","lwd","lty")
  el <- args[anames %in% el.names]
  args <- args[!(anames %in% el.names)]
  cmp <- RColorBrewer::brewer.pal(9,"Set1")
  cmp <- if (is.null(el$col)) c("black",cmp[c(2,1,3,5)]) else rep_len(el$col,5)
  lwd <- if (is.null(el$lwd)) 1 else el$lwd
  lty <- if (is.null(el$lty)) 1 else el$lty

  # Default limits of plot, with some space for the legend
  draw.legend <- !identical(legend,FALSE) && !is.null(legend)
  yy <- c(y,as.vector(x$fitted))
  if (thetalines){
    yy <- c(yy,as.vector(x$theta0),as.vector(x$theta2))
  }
  if (h > 0){
    yy <- c(yy,as.vector(frc$mean))
    if (thetalines){
      yy <- c(yy,as.vector(frc$frc.theta0),as.vector(frc$frc.theta2))
    }
  }
  yy <- range(yy,na.rm=TRUE)
  yy <- yy + c(-0.04,if (draw.legend) 0.12 else 0.04)*diff(yy)
  xx <- c(tt[1],if (h > 0) tf[h] else tt[n])

  # Allow user to override plot defaults
  if (is.null(args$main)){
    args$main <- paste("Theta method",thetadesc(x))
  }
  if (is.null(args$xlab)){
    args$xlab <- xlab
  }
  if (is.null(args$ylab)){
    args$ylab <- ""
  }
  if (is.null(args$xlim)){
    args$xlim <- xx
  }
  if (is.null(args$ylim)){
    args$ylim <- yy
  }
  if (is.null(args$xaxs) && h > 0){
    args$xaxs <- "i"
  }
  do.call(plot,c(list(x=NA,y=NA),args))

  # Forecast region
  if (h > 0){
    usr <- par("usr")
    rect(tt[n],usr[3],usr[2],usr[4],col="gray93",border=NA)
    box()
  }

  # Theta lines
  ftype <- if (h == 1) "p" else "l"
  if (thetalines){
    lines(tt,as.vector(x$theta0),col=cmp[3],lwd=lwd,lty=lty)
    lines(tt,as.vector(x$theta2),col=cmp[4],lwd=lwd,lty=lty)
    if (h > 0){
      lines(tf,as.vector(frc$frc.theta0),col=cmp[3],lwd=lwd,lty=lty,type=ftype,pch=20)
      lines(tf,as.vector(frc$frc.theta2),col=cmp[4],lwd=lwd,lty=lty,type=ftype,pch=20)
    }
  }

  # Data, fit and forecasts
  lines(tt,y,col=cmp[1],lwd=lwd,lty=lty)
  lines(tt,as.vector(x$fitted),col=cmp[2],lwd=lwd,lty=lty)
  if (h > 0){
    lines(c(tt[n],tf),c(as.vector(x$fitted)[n],as.vector(frc$mean)),col=cmp[2],lwd=2*lwd,lty=lty)
  }

  # Outliers
  has.out <- !is.null(x$outliers) && length(x$outliers) > 0
  if (has.out){
    points(tt[x$outliers],y[x$outliers],col=cmp[5],pch=19)
  }

  # Legend
  if (draw.legend){
    lgd <- c("Data",if (h > 0) "Fit & forecast" else "Fit")
    lcol <- cmp[1:2]
    if (thetalines){
      lgd <- c(lgd,"Theta 0","Theta 2")
      lcol <- c(lcol,cmp[3:4])
    }
    lty.l <- rep(lty,length(lgd))
    pch.l <- rep(NA,length(lgd))
    if (has.out){
      lgd <- c(lgd,"Outliers")
      lcol <- c(lcol,cmp[5])
      lty.l <- c(lty.l,NA)
      pch.l <- c(pch.l,19)
    }
    legend(if (is.character(legend)) legend[1] else "topleft",lgd,col=lcol,lty=lty.l,
           pch=pch.l,lwd=lwd,bty="n",cex=0.8,horiz=TRUE,text.width=NA)
  }

  invisible(NULL)

}

#' @export
#' @method summary theta
summary.theta <- function(object,...){
  print(object)
}

#' @export
#' @method print theta
print.theta <- function(x,...){

  mdl <- thetadesc(x)

  writeLines(paste("Theta method",mdl))
  writeLines("")
  writeLines("Parameters:")
  writeLines(paste("    Theta 0, intercept =",round(x$b[1],3)))
  if (x$exist[1]){
    writeLines(paste("    Theta 0, slope     =",round(x$b[2],3)))
  }
  writeLines(paste("    Theta 2, alpha     =",round(x$a[1],3)))
  if (x$exist[2]){
    writeLines(paste("    Season, gamma      =",round(x$g[1],3)))
  }
  writeLines("")
  writeLines("Initial states:")
  writeLines(paste("    Theta 2            =",round(x$a[2],3)))
  if (x$exist[2]){
    writeLines(paste("    Season             =",
                     paste(round(x$g[2:(x$m+1)],3),collapse=", ")))
  }
  writeLines("")
  writeLines(paste0("RMSE: ",round(sqrt(x$MSE),3)))
  invisible(x)

}

# Internal: short description of the fitted Theta model, used by print and plot
thetadesc <- function(x){
  if (all(x$exist == FALSE)){
    return("")
  }
  mdl <- NULL
  if (x$exist[1]){
    mdl <- "trend"
  }
  if (x$exist[2]){
    stp <- if (x$multiplicative) "multiplicative" else "additive"
    mdl <- c(mdl,paste0(stp," season (",x$m,")"))
  }
  return(paste0("(",paste(mdl,collapse=", "),")"))
}

#' @describeIn theta Wrapper function to use Theta with \code{\link[thief]{thief}}.
#'
#' @examples
#' \dontrun{
#' library(thief)
#' thief(referrals,forecastfunction=theta.thief)
#' }
#'
#' @export theta.thief
theta.thief <- function(y,h=NULL,...){

  # Remove level input from ellipsis
  ellipsis.args <- list(...)
  ellipsis.args$level <- NULL
  ellipsis.args$y <- y

  # Fit network
  fit <- do.call(theta,ellipsis.args)

  # Default h
  if (is.null(h)){
    h <- frequency(y)
  }

  # Forecast
  out <- forecast(fit,h)

  return(out)

}

opt.sfit <- function(ynt,costs,n,m,s,y,in.fit,multiplicative,outliers){
  # Optimise pure seasonal model and predict out-of-sample seasonality
  # Initial seasonal states refer to the first m observations, so align the seasonal
  # means of ynt (columns are seasonal periods) to the starting period s
  s.init <- unname(colMeans(ynt,na.rm=TRUE))[((s-1):(s+m-2)) %% m + 1]
  if (is.null(outliers)){
    g0 <- c(0.001,s.init)                         # Initialise seasonal model
    season.sample <- matrix(t(ynt),ncol=1)        # Transform back to vector
    season.sample <- season.sample[!is.na(season.sample)]
    X.out <- NULL
    n.out <- 0
  } else {
    n.out <- length(outliers)
    g0 <- c(0.001,s.init,rep(0,n.out))            # Initialise seasonal model
    if (multiplicative == TRUE){
      season.sample <- as.numeric(y / in.fit)
    } else {
      season.sample <- as.numeric(y - in.fit)
    }
    X.out <- array(0,c(n,n.out))
    X.out[outliers+n*(0:(n.out-1))] <- 1
  }
  opt <- optim(par=g0, cost.sfit, method = "Nelder-Mead", season.sample=season.sample,
               cost=costs, n=n, m=m, X.out = X.out, control = list(maxit = 2000))
  g <- opt$par
  sfit <- fun.sfit(season.sample,g,n,m,X.out)
  out.season <- sfit$outs
  in.season <- sfit$ins
  # Size of outliers
  if (n.out > 0){
    p <- g[(m+2):(m+1+n.out)]
  } else {
    p <- NULL
  }
  return(list("season"=out.season,"in.season"=in.season,"g"=g[1:(1+m)],"p"=p))
}

fun.sfit <- function(season.sample,g,n,m,X.out){
  # Fit pure seasonal model
  s.init <- g[2:(m+1)]
  season.fit <- c(s.init,rep(NA,n))
  for (i in 1:n){
    season.fit[i+m] <- season.fit[i] + g[1]*(season.sample[i] - season.fit[i])
  }
  if (!is.null(X.out)){
    n.out <- length(X.out[1,])
    g.out <- matrix(g[(m+2):(m+1+n.out)],ncol=1)
    X.out <- rbind(X.out, array(0,c(m,n.out)))
    season.fit <- season.fit + X.out %*% g.out
  }
  return(list("ins"=season.fit[1:n],"outs"=season.fit[(n+1):(n+m)]))
}

cost.sfit <- function(g,season.sample,cost,n,m,X.out){
  # Cost function of pure seasonal model
  err <- season.sample-fun.sfit(season.sample,g,n,m,X.out)$ins
  err <- cost.err(err,cost,NULL)
  if (g[1]<0 | g[1]>1){
    err <- 9*10^99
  }
  return(err)
}

fun.ses <- function(line,a,X.out=NULL){
  # Fit SES model on theta line
  n <- length(line)
  ses <- matrix(NA,nrow=n+1,ncol=1)
  ses[1] <- a[2]
  for (i in 2:(n+1)){
    ses[i] <- a[1]*line[i-1] + (1-a[1])*ses[i-1]
  }
  # If X.out is !null then model outliers
  if (!is.null(X.out)){
    n.out <- length(X.out[1,])
    a.out <- matrix(a[3:(n.out+2)],ncol=1)
    X.out <- rbind(X.out, array(0,c(1,n.out)))
    ses <- ses + X.out %*% a.out
  }
  return(list("ins"=ses[1:n],"outs"=ses[n+1]))
}

opt.ses <- function(line,cost,a0,theta,X.out=NULL){
  # Optimise SES on theta
  opt <- optim(par=a0, cost.ses, method = "Nelder-Mead", line=line, cost=cost,
               theta=theta, X.out=X.out, control = list(maxit = 2000))
  a <- opt$par
  return(a)
}

cost.ses <- function(a,line,cost,theta=2,X.out=NULL){
  # Cost function for SES optimisation
  err <- line-fun.ses(line,a,X.out)$ins
  err <- cost.err(err,cost,theta)
  if (!a[1]<0.99 | !a[1]>0.01){
    err <- 9*10^99
  }
  return(err)
}

opt.trnd <- function(y,X,cost,b0){
  # Optimise theta line 0
  opt <- optim(par=b0, cost.trnd, method = "Nelder-Mead", y=y, X=X,
               cost=cost, control = list(maxit = 2000))
  b <- opt$par
  return(b)
}

cost.trnd <- function(b,y,X,cost){
  # Theta 0 cost function
  err <- y-X%*%b
  err <- cost.err(err,cost,NULL)
  return(err)
}

cost.err <- function(err,cost,theta=NULL){
  # Cost calculation
  if (cost == "MAE"){
    err <- mean(abs(err))
  }
  if (cost == "MdAE"){
    err <- median(abs(err))
  }
  if (cost == "MSE"){
    err <- mean((err)^2)
  }
  if (cost == "MdSE"){
    err <- median((err)^2)
  }
  if (cost == "MTE"){
    err <- mean(abs((err)^theta))
  }
  if (cost == "MdTE"){
    err <- median(abs((err)^theta))
  }
  return(err)
}
