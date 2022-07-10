#' Classical time series decomposition
#'
#' Perform classical time series decomposition.
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
#' \item{\code{"mean"}}{: the mean of each seasonal period.}
#' \item{\code{"median"}}{: the median of each seasonal period.}
#' \item{\code{"pure.seasonal"}}{: estimate using a pure seasonal model.}
#' }
#' @param w percentage or number of observations to winsorise in the calculation of mean seasonal indices. If w>1 then it is the number of observations, otherwise it is a percentage. If \code{type != "mean"} then this is ignored.
#'
#' @return A list containing:
#' \itemize{
#' \item{\code{trend}}{: trend component.}
#' \item{\code{season}}{: season component.}
#' \item{\code{irregular}}{: irregular component.}
#' \item{\code{f.season}}{: forecasted seasonal component if \code{h>0}.}
#' \item{\code{g}}{: pure seasonal model parameters.}
#' }
#'
#' @references
#' Ord K., Fildes R., Kourentzes N. (2017) \href{https://kourentzes.com/forecasting/2017/10/16/new-forecasting-book-principles-of-business-forecasting-2e/}{Principles of Business Forecasting, 2e}. \emph{Wessex Press Publishing Co.}, p.106-111.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @keywords ts
#'
#' @examples
#' decomp(referrals)
#'
#' @export decomp

decomp <- function(y,m=NULL,s=NULL,trend=NULL,outplot=c(FALSE,TRUE),
                   decomposition=c("multiplicative","additive","auto"),
                   h=0,type=c("mean","median","pure.seasonal"),w=NULL)
{

  # Defaults
  outplot <- outplot[1]
  decomposition <- match.arg(decomposition,c("multiplicative","additive","auto"))
  type <- match.arg(type,c("mean","median","pure.seasonal"))

  # Get m (seasonality)
  if (is.null(m)){
    if (any(class(y) == "ts")){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }

  # Get starting period in seasonality if available
  if (is.null(s)){
    if (any(class(y) == "ts")){
      s <- start(y)
      s <- s[2]
    } else {
      s <- 1
    }
  }

  n <- length(y)

  if (min(y)<=0){
    decomposition <- "additive"
  }

  # If trend is not given then calculate CMA
  if (is.null(trend)){
    trend <- cmav(y=y,ma=m,fill=TRUE,outplot=FALSE,fast=TRUE)
  } else {
    if (n != length(trend)){
      stop("Length of series and trend input do not match.")
    }
  }

  if (decomposition == "auto"){
    if (mseastest(y,m)$is.multiplicative){
      decomposition <- "multiplicative"
    } else {
      decomposition <- "additive"
    }
  }

  if (decomposition == "multiplicative"){
    ynt <- y/trend
  } else {
    ynt <- y - trend
  }

  if (outplot != FALSE){
    yminmax.s <- range(ynt)
    yminmax.s <- yminmax.s + 0.04*c(-1,1)*diff(yminmax.s)
    yminmax.y <- range(y)
    yminmax.y <- yminmax.y + 0.04*c(-1,1)*diff(yminmax.y)
  }

  # Fill with NA start and end of season
  k <- m - (n %% m)
  ks <- s-1
  ke <- m - ((n+ks) %% m)
  ynt <- c(rep(NA,times=ks),as.vector(ynt),rep(NA,times=ke))
  ns <- length(ynt)/m
  ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
  colnames(ynt) <- paste("p",1:m,sep="")
  rownames(ynt) <- paste("s",1:ns,sep="")

  # If h>0 then produce forecasts of seasonality
  g <- NULL
  if (type=="mean"){
    # Calculate the seasonality as the overall mean
    if (!is.null(w)){
      ynt <- colWins(ynt,w)
    }
    season <- colMeans(ynt, na.rm=TRUE)
    if (h>0){
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h)])
    } else {
      f.season <- NULL
    }
    i.season <- rep(season,ns)
    i.season <- i.season[(ks+1):(ns*m-ke)]+y*0
  }
  if (type=="median"){
    # Calculate the seasonality as the overall median
    season <- array(NA,c(1,m))
    for (si in 1:m){
      season[si] <- median(ynt[,si], na.rm=TRUE)
    }
    if (h>0){
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h)])
    } else {
      f.season <- NULL
    }
    i.season <- rep(season,ns)
    i.season <- i.season[(ks+1):(ns*m-ke)] + y*0
  }
  if (type=="pure.seasonal"){
    # Seasonality is modelled with a pure seasonal smoothing
    sout <- decomp.opt.sfit(ynt=ynt,costs="MSE",n=n,m=m)
    g <- sout$g
    if (h>0){
      season <- sout$season
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-k+1):(m-k+h)])
      f.season <- rep(season, h %/% m + 1)[1:h]
    } else {
      f.season <- NULL
    }
    season.sample <- matrix(t(ynt),ncol=1)
    season.sample <- season.sample[!is.na(season.sample)]
    i.season <- as.vector(decomp.fun.sfit(season.sample,g,n,m)$ins) + y*0
  }

  # Convert f.season to ts object
  if (any(class(y) == "ts") && h>0){
    f.season <- ts(f.season,start=tsp(y)[2]+deltat(y),frequency=m)
  }

  if (decomposition == "multiplicative"){
    resid <- y - (trend*i.season)
  } else {
    resid <- y - (trend+i.season)
  }

  # Produce plots
  if (outplot == TRUE){
    # Write down the default values of par
    def.par <- par(no.readonly = TRUE);
    par(mfrow=c(4,1),mar=c(0,2,0,0),oma=c(2,2,2,2))

    # Series
    plot(1:n,as.vector(y),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),ylim=yminmax.y,xaxs = "i")
    if (decomposition == "multiplicative"){
      lines(1:n,trend*i.season,col="red",lty=1)
    } else {
      lines(1:n,trend+i.season,col="red",lty=1)
    }
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.y,yminmax.y[2:1]),border=NA,col="gray93")
    }
    mtext("Data",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    legend("topleft",c("Data","Reconstructed"),lty=c(1,1),lwd=c(1,1),col=c("black","red"),cex=0.8,bty="n",horiz=TRUE)

    # Trend
    plot(1:n,as.vector(trend),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.y,xaxs="i")
    mtext("Trend",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.y,yminmax.y[2:1]),border=NA,col="gray93")
    }

    # Season
    plot(1:n,as.vector(i.season),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.s,xaxs="i")
    mtext("Season",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.s,yminmax.s[2:1]),border=NA,col="gray93")
      lines((n+1):(n+h),f.season,col="blue")
    }
    if (decomposition == "multiplicative"){
      lines(1:(n+h),rep(1,n+h),lty=2,col="grey")
    } else {
      lines(1:(n+h),rep(0,n+h),lty=2,col="grey")
    }

    # Irregular
    yminmax.i = yminmax.y-mean(yminmax.y)
    plot(1:n,as.vector(resid),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.i,xaxs="i")
    mtext("Irregular",side=2,cex=0.8,padj=-2.5)
    text(1,yminmax.i[2],paste("RMSE:",round(sqrt(mean(resid^2)),2)),pos=4,cex=0.8)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.i,yminmax.i[2:1]),border=NA,col="gray93")
    }
    lines(1:n,rep(0,n),lty=2,col="grey")

    par(def.par)
  }

  return(list(trend=trend,season=i.season,irregular=resid,
              f.season=f.season,g=g))

}


decomp.opt.sfit <- function(ynt,costs,n,m){
  # Optimise pure seasonal model and predict out-of-sample seasonality
  g0 <- c(0.001,colMeans(ynt,na.rm=TRUE))       # Initialise seasonal model
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
