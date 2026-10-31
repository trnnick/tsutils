#' Classical time series decomposition
#'
#' Perform classical time series decomposition. Optionally, the seasonal component is only estimated when there is evidence of seasonality.
#'
#' @param y input time series. Can be \code{ts} object.
#' @param m seasonal period. If \code{y} is a \code{ts} object then the default is its frequency.
#' @param s starting period in the season. If \code{y} is a \code{ts} object then this is picked up from \code{y}.
#' @param trend vector of the level/trend of \code{y}. Use \code{NULL} to estimate internally.
#' @param outplot if \code{TRUE}, then provide a plot of the decomposed components.
#' @param decomposition type of decomposition. This can be \code{"multiplicative"}, \code{"additive"} or \code{"auto"}. If \code{y} contains non-positive values then this is forced to \code{"additive"}.
#' @param h forecast horizon for seasonal component.
#' @param type calculation for seasonal component:
#' \itemize{
#' \item \code{"mean"}: the mean of each seasonal period.
#' \item \code{"median"}: the median of each seasonal period.
#' \item \code{"pure.seasonal"}: estimate using a pure seasonal model.
#' }
#' @param w percentage or number of observations to winsorise in the calculation of mean seasonal indices. If w>1 then it is the number of observations, otherwise it is a percentage. If \code{type != "mean"} then this is ignored.
#' @param seastest if \code{TRUE}, then the detrended series is tested for seasonality with \code{\link{seastest}} and the seasonal component is only estimated if there is evidence of seasonality. Otherwise, the seasonal component is constant (1 for multiplicative and 0 for additive decomposition). If \code{FALSE} the seasonal component is always estimated.
#' @param alpha significance level for the seasonality test. Used only if \code{seastest=TRUE}.
#' @param ... additional arguments to control the plot:
#' \itemize{
#' \item \code{main}: title of the plot.
#' \item \code{xlab}: label of the time axis.
#' \item \code{ylab}: labels of the four panels.
#' \item \code{col}: colours of the data and components, the reconstructed series, and the seasonal forecast, recycled as needed.
#' \item \code{lwd}: line width.
#' \item \code{las}, \code{cex.axis}, \code{cex.lab}: axis label orientation, size of axis annotation and of panel labels.
#' }
#'
#' @return A list containing:
#' \itemize{
#' \item \code{trend}: trend component.
#' \item \code{season}: season component.
#' \item \code{irregular}: irregular component.
#' \item \code{f.season}: forecasted seasonal component if \code{h>0}.
#' \item \code{g}: pure seasonal model parameters.
#' \item \code{decomposition}: type of decomposition used.
#' \item \code{season.exist}: \code{TRUE}/\code{FALSE} result of the seasonality test, or \code{NULL} if not tested.
#' \item \code{season.pval}: p-value of the seasonality test, or \code{NULL} if not tested.
#' }
#'
#' @references
#' Ord K., Fildes R., Kourentzes N. (2017) \href{https://kourentzes.com/forecasting/2017/10/16/new-forecasting-book-principles-of-business-forecasting-2e/}{Principles of Business Forecasting, 2e}. \emph{Wessex Press Publishing Co.}, p.106-111.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{seastest}}, \code{\link{seasplot}}, \code{\link{cmav}}.
#'
#' @keywords ts
#'
#' @examples
#' decomp(referrals)
#' decomp(referrals,outplot=TRUE,h=12,seastest=TRUE)
#'
#' @export decomp

decomp <- function(y,m=NULL,s=NULL,trend=NULL,outplot=c(FALSE,TRUE),
                   decomposition=c("multiplicative","additive","auto"),
                   h=0,type=c("mean","median","pure.seasonal"),w=NULL,
                   seastest=c(FALSE,TRUE),alpha=0.05,...)
{

  # Defaults
  outplot <- isTRUE(as.logical(outplot[1]))
  decomposition <- match.arg(decomposition,c("multiplicative","additive","auto"))
  type <- match.arg(type,c("mean","median","pure.seasonal"))

  n <- length(y)
  if (!is.null(trend) && n != length(trend)){
    stop("Length of series and trend input do not match.")
  }

  # Detrend and organise in seasonal matrix
  sm <- seasmat(y=y,m=m,s=s,trend=TRUE,decomposition=decomposition,cma=trend)
  ynt <- sm$season
  m <- sm$m
  s <- sm$s
  trend <- sm$trend
  decomposition <- sm$decomposition
  colnames(ynt) <- paste0("p",1:m)
  ns <- nrow(ynt)
  ks <- s - 1
  ke <- ns*m - n - ks

  # Test for seasonality
  season.exist <- season.pval <- NULL
  if (isTRUE(seastest)){
    test <- seasfriedman(ynt,alpha=alpha)
    if (is.na(test$pvalue)){
      warning("Seasonality could not be tested (at least two complete seasons and m > 1 are needed). The seasonal component is estimated.")
    } else {
      season.exist <- test$is.seasonal
      season.pval <- test$pvalue
    }
  }

  # Estimate seasonal component
  g <- NULL
  if (isFALSE(season.exist) || type != "pure.seasonal"){
    if (isFALSE(season.exist)){
      # No evidence of seasonality, seasonal component is neutral
      season <- rep(if (decomposition == "multiplicative") 1 else 0, m)
    } else if (type == "mean"){
      # Calculate the seasonality as the overall mean
      if (!is.null(w)){
        ynt <- colWins(ynt,w)
      }
      season <- colMeans(ynt,na.rm=TRUE)
    } else {
      # Calculate the seasonality as the overall median
      season <- unname(apply(ynt,2,median,na.rm=TRUE))
    }
    i.season <- rep(season,ns)[(ks+1):(ns*m-ke)] + y*0
    f.season <- season[((ks+n+seq_len(h)-1) %% m) + 1]
  } else {
    # Seasonality is modelled with a pure seasonal smoothing
    sout <- decomp.opt.sfit(ynt=ynt,costs="MSE",n=n,m=m,s=s)
    g <- sout$g
    season.sample <- as.vector(t(ynt))
    season.sample <- season.sample[!is.na(season.sample)]
    i.season <- as.vector(decomp.fun.sfit(season.sample,g,n,m)$ins) + y*0
    f.season <- rep(sout$season,length.out=h)
  }
  f.season <- as.vector(f.season)
  if (h > 0){
    if (inherits(y,"ts")){
      f.season <- ts(f.season,start=tsp(y)[2]+deltat(y),frequency=m)
    }
  } else {
    f.season <- NULL
  }

  # Irregular component
  if (decomposition == "multiplicative"){
    resid <- y - (trend*i.season)
  } else {
    resid <- y - (trend+i.season)
  }

  # Produce plots
  if (outplot){
    decompdraw(y=y,trend=trend,i.season=i.season,resid=resid,f.season=f.season,h=h,
               decomposition=decomposition,season.exist=season.exist,season.pval=season.pval,...)
  }

  return(list(trend=trend,season=i.season,irregular=resid,f.season=f.season,g=g,
              decomposition=decomposition,season.exist=season.exist,season.pval=season.pval))

}

# Internal: plot of decomposed components for decomp
#' @noRd
decompdraw <- function(y,trend,i.season,resid,f.season,h,decomposition,
                       season.exist,season.pval,...){

  n <- length(y)
  is.mult <- decomposition == "multiplicative"

  # User options
  dots <- list(...)
  known <- c("main","xlab","ylab","col","lwd","las","cex.axis","cex.lab")
  dnames <- names(dots)
  if (is.null(dnames)){dnames <- rep("",length(dots))}
  if (any(!(dnames %in% known))){
    warning("Ignoring unsupported plotting arguments: ",
            paste(dnames[!(dnames %in% known)],collapse=", "),call.=FALSE)
  }
  getopt <- function(nm,default){
    if (is.null(dots[[nm]])) default else dots[[nm]]
  }
  cmp <- rep_len(getopt("col",c("black",brewer.pal(3,"Set1")[1:2])),3)
  lwd <- getopt("lwd",1)
  las <- getopt("las",0)
  cex.axis <- getopt("cex.axis",0.8)
  cex.lab <- getopt("cex.lab",0.8)
  panel <- getopt("ylab",c("Data","Trend","Season","Irregular"))
  if (length(panel) != 4){
    warning("ylab must contain four labels, one for each panel. Using defaults.",call.=FALSE)
    panel <- c("Data","Trend","Season","Irregular")
  }
  main <- getopt("main",paste0("Classical decomposition (",decomposition,")"))
  xlab <- getopt("xlab",if (inherits(y,"ts")) "Time" else "Period")

  # Time index, including forecast horizon
  if (inherits(y,"ts")){
    tt <- as.vector(time(y))
    dt <- deltat(y)
  } else {
    tt <- 1:n
    dt <- 1
  }
  tf <- tt[n] + seq_len(h)*dt
  xlim <- c(tt[1],if (h > 0) tf[h] else tt[n])

  # Plot ranges
  fit <- if (is.mult) trend*i.season else trend+i.season
  yrange <- function(x,ref=NULL){
    r <- range(c(x,ref),na.rm=TRUE)
    if (diff(r) == 0){
      r <- r + c(-1,1)*max(abs(r[1])*0.1,0.1)
    }
    r + 0.04*c(-1,1)*diff(r)
  }
  ref.s <- if (is.mult) 1 else 0
  yl.y <- yrange(c(as.vector(y),as.vector(fit),as.vector(trend)))
  yl.s <- yrange(c(as.vector(i.season),as.vector(f.season)),ref.s)
  yl.i <- c(-1,1)*max(diff(yl.y)/2,1.04*abs(as.vector(resid)),na.rm=TRUE)

  # Plot layout, restoring user settings on exit
  def.par <- par(no.readonly=TRUE)
  on.exit(par(def.par))
  layout(matrix(1:4,ncol=1),heights=c(1.4,1,1,1))
  par(mar=c(0,4.5,0,1),oma=c(3.5,0,3,0))

  panel.frame <- function(ylim,label,last=FALSE){
    plot(NA,NA,xlim=xlim,ylim=ylim,xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i")
    if (h > 0){
      rect(tt[n],ylim[1]-diff(ylim),xlim[2]+dt,ylim[2]+diff(ylim),col="gray93",border=NA)
    }
    abline(v=axTicks(1),col="gray90")
    # Drop labels too close to the panel edges, as they collide with the neighbouring panels
    at <- axTicks(2)
    usr <- par("usr")[3:4]
    at <- at[at >= usr[1] + 0.06*diff(usr) & at <= usr[2] - 0.06*diff(usr)]
    axis(2,at=at,labels=format(at,scientific=FALSE,trim=TRUE,drop0trailing=TRUE),
         cex.axis=cex.axis,las=las)
    mtext(label,side=2,line=3,cex=cex.lab)
    if (last){
      axis(1,cex.axis=cex.axis,las=las)
    }
    box()
  }

  # Data and reconstructed series
  panel.frame(yl.y,panel[1])
  lines(tt,as.vector(y),col=cmp[1],lwd=lwd)
  lines(tt,as.vector(fit),col=cmp[2],lwd=lwd)
  legend("topleft",c("Data","Reconstructed"),lty=1,lwd=lwd,col=cmp[1:2],
         cex=cex.lab,bty="n",horiz=TRUE)

  # Trend
  panel.frame(yl.y,panel[2])
  lines(tt,as.vector(trend),col=cmp[1],lwd=lwd)

  # Season
  panel.frame(yl.s,panel[3])
  abline(h=ref.s,lty=2,col="grey")
  lines(tt,as.vector(i.season),col=cmp[1],lwd=lwd)
  if (h > 0){
    lines(c(tt[n],tf),c(i.season[n],as.vector(f.season)),col=cmp[3],lwd=lwd)
  }
  if (!is.null(season.exist)){
    legend("topleft",paste0(if (season.exist) "Seasonal" else "Nonseasonal",
                            " (p-val: ",round(season.pval,3),")"),bty="n",cex=cex.lab)
  }

  # Irregular
  panel.frame(yl.i,panel[4],last=TRUE)
  abline(h=0,lty=2,col="grey")
  lines(tt,as.vector(resid),col=cmp[1],lwd=lwd)
  legend("topleft",paste("RMSE:",format(signif(sqrt(mean(resid^2,na.rm=TRUE)),4))),
         bty="n",cex=cex.lab)

  # Titles
  mtext(main,side=3,outer=TRUE,line=1,font=2)
  mtext(xlab,side=1,outer=TRUE,line=2.3,cex=cex.lab)

  invisible(NULL)

}

decomp.opt.sfit <- function(ynt,costs,n,m,s=1){
  # Optimise pure seasonal model and predict out-of-sample seasonality
  # Initialise seasonal model, with initial states aligned to the first observation (period s)
  g0 <- c(0.001,colMeans(ynt,na.rm=TRUE)[((s-1):(s+m-2)) %% m + 1])
  season.sample <- matrix(t(ynt),ncol=1)        # Transform back to vector
  season.sample <- season.sample[!is.na(season.sample)]
  opt <- optim(par=g0, decomp.cost.sfit, method = "Nelder-Mead", season.sample=season.sample,
               cost=costs, n=n, m=m, control = list(maxit = 2000))
  g <- opt$par
  season <- decomp.fun.sfit(season.sample,g,n,m)$outs
  return(list("season"=season,"g"=g))
}

decomp.fun.sfit <- function(season.sample,g,n,m){
  # Fit pure seasonal model
  s.init <- g[2:(m+1)]
  season.fit <- c(s.init,rep(NA,n))
  for (i in 1:n){
    season.fit[i+m] <- season.fit[i] + g[1]*(season.sample[i] - season.fit[i])
  }
  return(list("ins"=season.fit[1:n],"outs"=season.fit[(n+1):(n+m)]))
}

decomp.cost.sfit <- function(g,season.sample,cost,n,m){
  # Cost function of pure seasonal model
  err <- season.sample-decomp.fun.sfit(season.sample,g,n,m)$ins
  err <- decomp.cost.err(err,cost)
  if (g[1]<0 | g[1]>1){
    err <- 9*10^99
  }
  return(err)
}

decomp.cost.err <- function(err,cost){
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
  return(err)
}
