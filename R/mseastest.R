#' Multiplicative seasonality test
#'
#' correlation based multiplicative seasonality test.
#'
#' @param y input time series. Can be \code{ts} object.
#' @param m seasonal period. If \code{y} is a \code{ts} object then the default is its frequency.
#' @param type type of correlation
#' \itemize{
#' \item{\code{"pearson"}}
#' \item{\code{"spearman"}}
#' \item{\code{"kendall"}}
#' }
#' @param cma input precalculated level/trend for the analysis. Use \code{NULL} to calculate internally.
#' @param sn number of seasonal periods of decreasing magnitude to consider for the test.
#' @param alpha significance level.
#' @param outplot type of output plot:
#' \itemize{
#' \item \code{0}: none.
#' \item \code{1}: scatter plot.
#' \item \code{2}: time series plot.
#' }
#'
#' @return A list with the following:
#' \itemize{
#' \item \code{is.multiplicative}: if \code{TRUE} the test found evidence of multiplicative seasonality.
#' \item \code{statistic}: test statistic.
#' \item \code{pvalue}: p-value of the test.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @keywords ts htest internal
#'
#' @examples
#' mseastest(referrals)
#'
#' @export mseastest

mseastest <- function(y,m=NULL,type=c("pearson","spearman","kendall"),cma=NULL,
                      sn=1, alpha=0.05,outplot=c(0,1,2))
{

  # Defaults
  type <- match.arg(type,c("pearson","spearman","kendall"))
  outplot <- outplot[1]

  # Get m (seasonality)
  if (is.null(m)){
    if (any(class(y) == "ts")){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }

  # Test sn
  if (sn>m){
    sn <- m
  }

  # Calculate CMA if not provided
  if (is.null(cma)){
    cma <- cmav(y,ma=m)
  }

  # Remove trend from time series
  seas <- (y - cma)

  # Convert to seasonal matrix
  n <- length(y)
  k <- m - (n %% m)
  seas <- c(as.vector(seas),rep(NA,times=k))
  ns <- length(seas)/m
  seas <- matrix(seas,nrow=ns,ncol=m,byrow=TRUE)
  colnames(seas) <- paste("m",1:m,sep="")

  # Find size of seasonality and adjust for direction (of mean size)
  mag <- colMeans(seas, na.rm = TRUE)
  idx <- order(abs(mag),decreasing=TRUE)
  ssign <- (mag<0)*-2+1
  seas <- seas*matrix(rep(ssign,ns),nrow=ns,ncol=m,byrow=TRUE)
  idx <- idx[1:sn]

  # Measure correlation for sn seasonal periods
  cor.size <- array(NA,c(sn,1))
  cor.pvalue <- array(NA,c(sn,1))
  for (i in 1:sn){
    X <- seas[,idx[i]]
    X <- X[!is.na(X)]
    Y <- cma[seq(idx[i],n,m)]
    Y <- Y[!is.na(X)]
    if (length(X)<3){
      stop("Not enough seasons to test for multiplicative seasonality.")
    }
    test <- cor.test(X,Y,method=type,alternative="greater")
    cor.size[i] <- test$estimate
    cor.pvalue[i] <- test$p.value
  }

  # Aggregate cases
  p <- median(cor.pvalue)
  c <- median(cor.size)

  # Reset p-value for negative correlations
  if (c <= 0){
    p <- 1
  }

  # Do test
  if (p <= alpha){
    is.multiplicative <- TRUE
    H <- "Multiplicative"
  } else {
    is.multiplicative <- FALSE
    H <- "Additive"
  }

  # Plot if requested
  if (outplot == 1){
    plot.title = paste(H, " seasonality (pval: ",round(p,3),")",sep="")
    xmin <- as.vector(seas[,idx[1:sn]])
    xmin <- xmin[!is.na(xmin)]
    xmax <- max(xmin)
    xmin <- min(xmin)
    xminmax <- c(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin))
    ymin <- as.vector(mapply(Vectorize(seq),idx[1:sn],idx[1:sn]+m*ns,m))
    ymin <- cma[ymin[ymin<n]]
    ymax <- max(ymin)
    ymin <- min(ymin)
    yminmax <- c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    cmp <- rainbow(sn)
    X <- seas[,idx[1]]
    X <- X[!is.na(X)]
    Y <- cma[seq(idx[1],n,m)]
    Y <- Y[!is.na(X)]
    plot(X, Y, col=cmp[1], pch=20, xlab="Seasonal value", ylab="Level value",
         main=plot.title, xlim=xminmax, ylim=yminmax)
    if (sn>1){
      text(X[1],Y[1],"s1",col=cmp[1], pos=3, cex=0.7)
      for (i in 2:sn){
        X <- seas[,idx[i]]
        X <- X[!is.na(X)]
        Y <- cma[seq(idx[i],n,m)]
        Y <- Y[!is.na(X)]
        points(X,Y,col=cmp[i],pch=20)
        text(X[1],Y[1],paste("s",i,sep=""),col=cmp[i], pos=3, cex=0.7)
      }
    }
  }

  if (outplot == 2){
    cmp <- rainbow(sn)
    plot.title = paste(H, " seasonality (pval: ",round(p,3),")",sep="")
    plot(1:n,as.vector(y),type="l",main=plot.title,ylab="",xlab="Period")
    lines(as.vector(cma),lty=1,col="black",lwd=2)
    midx <- seq(idx[1],n,m)
    points(midx,y[midx],col=cmp[1],pch=16)
    if (sn>1){
      if (y[midx[1]]>cma[midx[1]]){
        text(midx[1],y[midx[1]],"s1",col=cmp[1], pos=3, cex=0.7)
      }
      for (i in 2:sn){
        midx <- seq(idx[i],n,m)
        points(midx,(y)[midx],col=cmp[i],pch=16)
        if (y[midx[1]]>cma[midx[1]]){
          text(midx[1],y[midx[1]],paste("s",i,sep=""),col=cmp[i], pos=3, cex=0.7)
        } else {
          text(midx[1],y[midx[1]],paste("s",i,sep=""),col=cmp[i], pos=1, cex=0.7)
        }
      }
    }
  }

  return(list("is.multiplicative"=is.multiplicative,"statistic"=c,"pvalue"=p))

}
