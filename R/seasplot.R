#' Seasonal plots with simplistic trend/season tests
#'
#' Construct seasonal plots of various styles for a given time series. The series can automatically detrended as needed.
#'
#' @param y input time series. Can be \code{ts} object.
#' @param m seasonal period. If \code{y} is a \code{ts} object then the default is its frequency.
#' @param s starting period in the season. If \code{y} is a \code{ts} object then this is picked up from \code{y}.
#' @param trend if \code{TRUE}, then presence of trend is assumed and removed. If \code{FALSE} no trend is assumed. Use \code{NULL} to identify automatically.
#' @param colour single colour override for plots.
#' @param alpha significance level for statistical tests.
#' @param outplot type of seasonal plot
#' \itemize{
#' \item \code{0}: none.
#' \item \code{1}: seasonal diagram.
#' \item \code{2}: seasonal boxplots.
#' \item \code{3}: seasonal subseries.
#' \item \code{4}: seasonal distribution.
#' \item \code{5}: seasonal density.
#' }
#' @param decomposition type of seasonal decomposition. This can be \code{"multiplicative"}, \code{"additive"} or \code{"auto"}. If \code{y} contains non-positive values then this is forced to \code{"additive"}.
#' @param cma input precalculated level/trend for the analysis. This overrides \code{trend=NULL}.
#' @param labels external labels for the seasonal periods. Use \code{NULL} for none. If \code{length(labels) < m}, then this input is ignored.
#' @param ... additional arguments passed to plotting functions. For example, use \code{main=""} to replace the title.
#'
#' @return An object of class \code{seasexpl} containing:
#' \itemize{
#' \item \code{season}: matrix of (detrended) seasonal elements.
#' \item \code{season.exist}: \code{TRUE}/\code{FALSE} results of seasonality test.
#' \item \code{season.pval}: p-value of seasonality test (Friedman test).
#' \item \code{trend}: CMA estimate (using \code{\link{cmav}}) or \code{NULL} if \code{trend=FALSE}.
#' \item \code{trend.exist}: \code{TRUE}/\code{FALSE} results of trend test.
#' \item \code{trend.pval}: p-value of trend test (Cox-Stuart).
#' \item \code{decomposition}: type of decomposition used.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @keywords ts
#'
#' @examples
#' seasplot(referrals,outplot=1)
#'
#' @export seasplot

seasplot <- function(y,m=NULL,s=NULL,trend=NULL,colour=NULL,alpha=0.05,
                     outplot=c(1,0,2,3,4,5),decomposition=c("multiplicative","additive","auto"),
                     cma=NULL,labels=NULL,...)
{

  # Defaults
  decomposition <- match.arg(decomposition,c("multiplicative","additive","auto"))
  outplot <- outplot[1]

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
      # Temporal aggregation can mess-up s, so override if needed
      if (is.na(s)){s<-1}
    } else {
      s <- 1
    }
  }

  # Make sure that labels argument is fine
  if (!is.null(labels)){
    if (length(labels) < m){
      labels  <- NULL
      warning("Incorrect number of labels. It must be equal to m. Defaulting to labels=NULL.")
    } else {
      labels <- labels[1:m]
    }
  }

  if (is.null(labels)){
    labels <- paste(1:m)
  }

  n <- length(y)

  if (min(y)<=0){
    decomposition <- "additive"
  }

  # Override trend if cma is given
  if (!is.null(cma)){
    trend <- NULL
    if (n != length(cma)){
      stop("Length of series and cma do not match.")
    }
  }

  # Calculate CMA
  if ((is.null(trend) || trend == TRUE) && (is.null(cma))){
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

  if (trend == TRUE){
    is.mult <- mseastest(y,m)$is.multiplicative
    if (decomposition == "auto"){
      if (is.mult == TRUE){
        decomposition <- "multiplicative"
      } else {
        decomposition <- "additive"
      }
    }
    if (decomposition == "multiplicative"){
      ynt <- y/cma
    } else {
      ynt <- y - cma
    }
    title.trend <- "(Detrended)"
  } else {
    if (decomposition == "auto"){
      decomposition <- "none"
    }
    ynt <- y
    title.trend <- ""
    cma <- NULL
  }

  ymin <- min(ynt)
  ymax <- max(ynt)
  ymid <- median(ynt)

  # Fill with NA start and end of season
  k <- m - (n %% m)
  ks <- s-1
  ke <- m - ((n+ks) %% m)
  if (ke == m){ke <- 0}
  if (ks == m){ks <- 0}
  ynt <- c(rep(NA,times=ks),as.vector(ynt),rep(NA,times=ke))
  ns <- length(ynt)/m
  ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
  colnames(ynt) <- labels
  rownames(ynt) <- paste("s",1:ns,sep="")

  # Check seasonality with Friedman
  if (m>1 && (length(y)/m)>=2){
      # Check if ynt is a constant, if it is do not get it through friedman, as it fails
      if (diff(range(colSums(ynt,na.rm=TRUE))) <= .Machine$double.eps ^ 0.5){
          season.pval <- 1
      } else {
          season.pval <- friedman.test(ynt)$p.value
      }
      season.exist <- season.pval <= alpha
    if (season.exist==TRUE){
      title.season <- "Seasonal"
    } else {
      title.season <- "Nonseasonal"
    }
  } else {
    season.pval <- NULL
    season.exist <- NULL
    title.season <- "Nonseasonal"
  }

  # Produce plots
  if (outplot != 0){
    yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
    yminmax <- c(-1,1)*max(abs(ymid-yminmax))+ymid
    if (is.null(season.pval)){
      plottitle <- paste(title.trend, "\n", title.season,sep="")
    } else {
      plottitle <- paste(title.trend, "\n", title.season,
                   " (p-val: ",round(season.pval,3),")",sep="")
    }
    # Allow user to override plot defaults
    args <- list(...)
    if (!("main" %in% names(args))){
      addtitle <- TRUE
    } else {
      addtitle <- FALSE
    }
    if (!("xlab" %in% names(args))){
      args$xlab <- "Period"
    }
    if (!("ylab" %in% names(args))){
      args$ylab <- ""
    }
    if (!("yaxs" %in% names(args))){
      args$yaxs <- "i"
    }
    if (!("ylim" %in% names(args))){
      args$ylim <- yminmax
    }
    # Remaining defaults
    args$x <- args$y <- NA
   }

  if (outplot == 1){
    # Conventional seasonal diagramme
    if (is.null(colour)){
      cmp <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu")[4:8])(ns)
    } else {
      cmp <- rep(colour,times=ns)
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    if (addtitle){
      args$main <- paste0("Seasonal plot ",main=plottitle)
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    for (i in 1:ns){
      lines(ynt[i,],type="l",col=cmp[i])
    }
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
    legend("topleft",c("Oldest","Newest"),col=c(cmp[1],cmp[ns]),lty=1,bty="n",lwd=2,cex=0.7)
    axis(1,at=1:m,labels=labels)
  }
  if (outplot == 2){
    # Seasonal boxplots
    if (is.null(colour)){
      cmp <- RColorBrewer::brewer.pal(3,"Set3")[1]
    } else {
      cmp <- colour
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (addtitle){
      args$main <- paste0("Seasonal boxplot ",main=plottitle)
    }
    args$x <- ynt
    args$col <- cmp
    # Produce plot
    do.call(boxplot,args)
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
  }
  if (outplot == 3){
    # Subseries plots
    if (is.null(colour)){
      cmp <- RColorBrewer::brewer.pal(3,"Set1")
    } else {
      cmp <- colour
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m*ns)
    }
    if (addtitle){
      args$main <- paste0("Seasonal subseries ",main=plottitle)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    lines(c(1,ns),median(ynt[,1],na.rm=TRUE)*c(1,1),col=cmp[1],lwd=2)
    for (i in 1:m){
      lines((1+(i-1)*ns):(ns+(i-1)*ns),ynt[,i],type="o",col=cmp[2],pch=20,cex=0.75)
      lines(c((1+(i-1)*ns),(ns+(i-1)*ns)),median(ynt[,i],na.rm=TRUE)*c(1,1),col=cmp[1],lwd=2)
      if (i < m){
        lines(c(1,1)*i*ns+0.5,yminmax,col="gray")
      }
    }
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    axis(1,at=seq(0.5+(ns/2),0.5+m*ns-ns/2,ns),labels=labels)
  }
  if (outplot == 4){
    # Seasonal distribution
    if (is.null(colour)){
      cmp <- "royalblue"
    }
    qntl <- matrix(NA,nrow=9,ncol=m)
    for (i in 1:m){
      qntl[,i] <- quantile(ynt[!is.na(ynt[,i]),i], c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1))
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (addtitle){
      args$main <- paste0("Seasonal distribution ",main=plottitle)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    polygon(c(1:m,rev(1:m)),c(qntl[9,],rev(qntl[1,])),col=gray(0.8),border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[8,],rev(qntl[2,])),col="lightblue",border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[7,],rev(qntl[3,])),col="skyblue",border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[6,],rev(qntl[4,])),col="skyblue",border=NA)
    lines(1:m,qntl[5,],col=cmp,lwd=2)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    legend("topleft",c("Median","25%-75%","10%-90%","MinMax"),col=c(cmp,"skyblue","lightblue",gray(0.8)),lty=1,bty="n",lwd=2,cex=0.7)
    axis(1,at=1:m,labels=labels)
    box()
  }
  if (outplot == 5){
    dnst <- matrix(NA,nrow=m,ncol=512)
    for (i in 1:m){
       tmp <- density(ynt[!is.na(ynt[,i]),i], bw = "SJ", n = 512, from = yminmax[1], to = yminmax[2])
       dnst[i,] <- tmp$y
       if (i == 1){
        llc <- tmp$x
       }
    }
    cmp <- c(rgb(0,0,0,0), colorRampPalette((RColorBrewer::brewer.pal(9,"Blues")[3:9]))(100))
    # Additional plot options
    if (addtitle){
      args$main <- paste0("Seasonal density ",main=plottitle)
    }
    args$xaxt <- "n"
    args$col <- cmp
    args$x <- 1:m
    args$y <- llc
    args$z <- dnst
    # Produce plot
    do.call(image,args)
    lines(colMeans(ynt,na.rm=TRUE),type="o",lty=1,bg=RColorBrewer::brewer.pal(3,"Set1")[1],pch=21,cex=1.1,lwd=2)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    box()
    axis(1,at=1:m,labels=labels)
  }

  out <- list(season=ynt,season.exist=season.exist,season.pval=season.pval,
              trend=cma,trend.exist=trend.exist,trend.pval=trend.pval,decomposition=decomposition)
  out <- structure(out,class="seasexpl")
  return(out)

}

#' @export
#' @method summary seasexpl
summary.seasexpl <- function(object,...){
  print(object)
}

#' @export
#' @method print seasexpl
print.seasexpl <- function(x,...){
  cat("Results of statistical testing\n")
  if (!is.null(x$trend.exist)){
    cat(paste0("Evidence of trend: ",x$trend.exist, "  (pval: ",round(x$trend.pval,3),")\n"))
  } else {
    cat("Presence of trend not tested.\n")
  }
  cat(paste0("Evidence of seasonality: ",x$season.exist, "  (pval: ",round(x$season.pval,3),")\n"))
}
