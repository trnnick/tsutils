#' Seasonal plots with simplistic trend/season tests
#'
#' Construct seasonal plots of various styles for a given time series. The series can be automatically detrended as needed.
#'
#' @param y input time series. Can be \code{ts} object.
#' @param m seasonal period. If \code{y} is a \code{ts} object then the default is its frequency.
#' @param s starting period in the season. If \code{y} is a \code{ts} object then this is picked up from \code{y}.
#' @param trend if \code{TRUE}, then presence of trend is assumed and removed. If \code{FALSE} no trend is assumed. Use \code{NULL} to identify automatically.
#' @param colour colour override for plots. A single colour, or a vector of colours, used as follows:
#' \itemize{
#' \item \code{outplot=1}: colours of the seasons, from oldest to newest. Two or more colours (other than one per season) are interpolated into a gradient.
#' \item \code{outplot=2}: colours of the boxes, recycled across seasonal periods.
#' \item \code{outplot=3}: colours of the medians and of the subseries, in this order. A single colour is used for both.
#' \item \code{outplot=4}: colour of the median. Optionally, three more colours for the 25\%-75\%, 10\%-90\% and min-max bands. Otherwise, the bands are lighter shades of the median colour.
#' \item \code{outplot=5}: colour of the highest density. Two or more colours are interpolated into the density palette.
#' }
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
#' @param cma input precalculated level/trend for the analysis. If \code{trend=TRUE} it is removed from \code{y} without testing, otherwise the presence of trend is tested on \code{cma} (overriding \code{trend=FALSE}).
#' @param labels external labels for the seasonal periods. Use \code{NULL} for none. If \code{length(labels) < m}, then this input is ignored.
#' @param legend if \code{TRUE} a legend is added to the plot, where applicable (\code{outplot=1} and \code{4}). Use \code{FALSE} to omit it, or a keyword such as \code{"topright"} to change its position.
#' @param ... additional arguments passed to plotting functions. Arguments that refer to the plotted data, \code{col}, \code{lwd}, \code{lty}, \code{pch}, \code{cex}, \code{bg} and \code{type}, are applied to the seasonal lines, boxes or points (\code{col} acts as \code{colour}). Axis arguments (\code{las}, \code{cex.axis}, \code{col.axis}, \code{font.axis}, \code{tck}, \code{tcl}, \code{mgp}) are also applied to the seasonal period axis. Everything else, such as \code{main}, \code{xlab}, \code{ylab}, \code{xlim}, \code{ylim}, \code{xaxt="n"} or \code{axes=FALSE}, is passed to the function that sets up the plot (\code{plot}, \code{boxplot} or \code{image}).
#'
#' @details
#' The seasonality test is the same as in \code{\link{seastest}}, which can be used to run it without producing a plot.
#'
#' @return An object of class \code{seasexpl} containing:
#' \itemize{
#' \item \code{season}: matrix of (detrended) seasonal elements.
#' \item \code{season.exist}: \code{TRUE}/\code{FALSE} results of seasonality test, or \code{NULL} if not tested.
#' \item \code{season.pval}: p-value of seasonality test (Friedman test), or \code{NULL} if not tested.
#' \item \code{trend}: CMA estimate (using \code{\link{cmav}}) or \code{NULL} if no trend was removed.
#' \item \code{trend.exist}: \code{TRUE}/\code{FALSE} results of trend test.
#' \item \code{trend.pval}: p-value of trend test (Cox-Stuart).
#' \item \code{decomposition}: type of decomposition used. \code{"none"} if no trend was removed.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{seastest}}, \code{\link{decomp}}.
#'
#' @keywords ts
#'
#' @examples
#' seasplot(referrals,outplot=1)
#' seasplot(referrals,outplot=1,colour=c("grey","red"),lwd=2,legend="bottomright")
#' seasplot(referrals,outplot=3,labels=month.abb,las=2,main="Referrals",ylab="Index")
#'
#' @export seasplot

seasplot <- function(y,m=NULL,s=NULL,trend=NULL,colour=NULL,alpha=0.05,
                     outplot=c(1,0,2,3,4,5),decomposition=c("multiplicative","additive","auto"),
                     cma=NULL,labels=NULL,legend=TRUE,...)
{

  outplot <- outplot[1]
  if (!(outplot %in% 0:5)){
    stop("outplot must be one of 0, 1, 2, 3, 4 or 5.")
  }

  # Detrend and organise in seasonal matrix
  sm <- seasmat(y=y,m=m,s=s,trend=trend,alpha=alpha,decomposition=decomposition,
                cma=cma,labels=labels)
  ynt <- sm$season
  m <- sm$m
  labels <- sm$labels
  ns <- nrow(ynt)

  # Test for seasonality
  test <- seasfriedman(ynt,alpha=alpha)
  if (is.na(test$pvalue)){
    season.pval <- season.exist <- NULL
  } else {
    season.pval <- test$pvalue
    season.exist <- test$is.seasonal
  }

  # Produce plots
  if (outplot != 0){
    seasdraw(ynt=ynt,m=m,ns=ns,labels=labels,outplot=outplot,colour=colour,
             legend=legend,detrended=(sm$decomposition != "none"),
             season.exist=season.exist,season.pval=season.pval,...)
  }

  out <- list(season=ynt,season.exist=season.exist,season.pval=season.pval,
              trend=sm$trend,trend.exist=sm$trend.exist,trend.pval=sm$trend.pval,
              decomposition=sm$decomposition)
  out <- structure(out,class="seasexpl")
  return(out)

}

# Internal: draw the seasonal plots for seasplot
#' @noRd
seasdraw <- function(ynt,m,ns,labels,outplot,colour,legend,detrended,
                     season.exist,season.pval,...){

  # Plot range and title
  ymin <- min(ynt,na.rm=TRUE)
  ymax <- max(ynt,na.rm=TRUE)
  ymid <- median(ynt,na.rm=TRUE)
  yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
  yminmax <- c(-1,1)*max(abs(ymid-yminmax))+ymid
  if (diff(yminmax) == 0){
    yminmax <- ymid + c(-1,1)*max(1,abs(ymid))*0.1
  }
  title.trend <- if (detrended) "(Detrended)" else ""
  title.season <- if (isTRUE(season.exist)) "Seasonal" else "Nonseasonal"
  if (is.null(season.pval)){
    plottitle <- paste0(title.trend,"\n",title.season)
  } else {
    plottitle <- paste0(title.trend,"\n",title.season," (p-val: ",round(season.pval,3),")")
  }

  # Split user arguments between the plotted data elements, the period axis and the plot frame
  dots <- list(...)
  dnames <- names(dots)
  if (is.null(dnames)){dnames <- rep("",length(dots))}
  el.names <- c("col","lwd","lty","pch","cex","bg","type")
  ax.names <- c("las","cex.axis","col.axis","font.axis","tck","tcl","mgp")
  el <- dots[dnames %in% el.names]
  args <- dots[!(dnames %in% el.names)]
  ax <- dots[dnames %in% ax.names]
  if (is.null(colour) && !is.null(el$col)){
    colour <- el$col
  }
  el$col <- NULL
  draw.xaxis <- !identical(args$xaxt,"n") && !identical(args$axes,FALSE)
  draw.legend <- !identical(legend,FALSE) && !is.null(legend)
  legend.pos <- if (is.character(legend)) legend[1] else "topleft"

  # Fill in defaults for anything the user did not set
  setdef <- function(lst,def){
    for (nm in names(def)){
      if (is.null(lst[[nm]])){
        lst[nm] <- def[nm]
      }
    }
    return(lst)
  }
  ptype <- c("Seasonal plot ","Seasonal boxplot ","Seasonal subseries ",
             "Seasonal distribution ","Seasonal density ")[outplot]
  args <- setdef(args,list(main=paste0(ptype,plottitle),xlab="Period",ylab="",
                           yaxs="i",ylim=yminmax))
  yminmax <- args$ylim
  xaxis <- function(at){
    if (draw.xaxis){
      do.call(axis,c(list(side=1,at=at,labels=labels),ax))
    }
  }
  refline <- function(){
    abline(h=ymid,col="black",lty=2)
  }

  if (outplot == 1){
    # Conventional seasonal diagramme
    if (is.null(colour)){
      cmp <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu")[4:8])(ns)
    } else if (length(colour) == 1 || length(colour) == ns){
      cmp <- rep_len(colour,ns)
    } else {
      cmp <- colorRampPalette(colour)(ns)
    }
    el <- setdef(el,list(type="l",lty=1,lwd=1))
    args <- setdef(args,list(xlim=c(1,m),xaxs="i"))
    args$xaxt <- "n"
    do.call(plot,c(list(x=NA,y=NA),args))
    for (i in 1:ns){
      do.call(lines,c(list(x=1:m,y=ynt[i,],col=cmp[i]),el))
    }
    refline()
    if (draw.legend){
      legend(legend.pos,c("Oldest","Newest"),col=c(cmp[1],cmp[ns]),lty=el$lty,
             lwd=max(2,el$lwd),bty="n",cex=0.7)
    }
    xaxis(1:m)
  }

  if (outplot == 2){
    # Seasonal boxplots
    if (is.null(colour)){
      colour <- RColorBrewer::brewer.pal(3,"Set3")[1]
    }
    el$type <- NULL
    args <- setdef(args,list(xlim=c(0.5,m+0.5)))
    if (!draw.xaxis){
      args$xaxt <- "n"
    }
    do.call(boxplot,c(list(x=ynt,col=rep_len(colour,m)),el,args))
    refline()
  }

  if (outplot == 3){
    # Subseries plots
    if (is.null(colour)){
      cmp <- RColorBrewer::brewer.pal(3,"Set1")[1:2]
    } else {
      cmp <- rep_len(colour,2)
    }
    el <- setdef(el,list(type="o",pch=20,cex=0.75,lwd=1))
    args <- setdef(args,list(xlim=c(1,m*ns),xaxs="i"))
    args$xaxt <- "n"
    do.call(plot,c(list(x=NA,y=NA),args))
    if (m > 1){
      abline(v=(1:(m-1))*ns+0.5,col="gray")
    }
    for (i in 1:m){
      idx <- (1+(i-1)*ns):(ns+(i-1)*ns)
      do.call(lines,c(list(x=idx,y=ynt[,i],col=cmp[2]),el))
      lines(range(idx),median(ynt[,i],na.rm=TRUE)*c(1,1),col=cmp[1],lwd=2*el$lwd)
    }
    refline()
    xaxis(seq(0.5+(ns/2),0.5+m*ns-ns/2,ns))
  }

  if (outplot == 4){
    # Seasonal distribution
    if (is.null(colour)){
      cmp <- c("royalblue","skyblue","lightblue",gray(0.8))
    } else if (length(colour) >= 4){
      cmp <- colour[1:4]
    } else {
      cmp <- c(colour[1],colour.blend(colour[1],c(0.5,0.3,0.15)))
    }
    qntl <- apply(ynt,2,quantile,probs=c(0,0.1,0.25,0.5,0.75,0.9,1),na.rm=TRUE)
    el <- setdef(el,list(type="l",lty=1,lwd=2))
    args <- setdef(args,list(xlim=c(1,m),xaxs="i"))
    args$xaxt <- "n"
    do.call(plot,c(list(x=NA,y=NA),args))
    for (j in 1:3){
      polygon(c(1:m,m:1),c(qntl[8-j,],rev(qntl[j,])),col=cmp[5-j],border=NA)
    }
    do.call(lines,c(list(x=1:m,y=qntl[4,],col=cmp[1]),el))
    refline()
    if (draw.legend){
      legend(legend.pos,c("Median","25%-75%","10%-90%","MinMax"),col=cmp,lty=1,
             bty="n",lwd=c(el$lwd,6,6,6),cex=0.7)
    }
    xaxis(1:m)
    box()
  }

  if (outplot == 5){
    # Seasonal density
    dnst <- matrix(0,nrow=m,ncol=512)
    llc <- seq(yminmax[1],yminmax[2],length.out=512)
    for (i in 1:m){
      x <- ynt[!is.na(ynt[,i]),i]
      if (length(x) >= 2){
        tmp <- tryCatch(density(x,bw="SJ",n=512,from=yminmax[1],to=yminmax[2]),
                        error=function(e){
                          tryCatch(density(x,n=512,from=yminmax[1],to=yminmax[2]),
                                   error=function(e){NULL})})
        if (!is.null(tmp)){
          dnst[i,] <- tmp$y
        }
      }
    }
    if (is.null(colour)){
      cmp <- RColorBrewer::brewer.pal(9,"Blues")[3:9]
    } else if (length(colour) == 1){
      cmp <- c(colour.blend(colour,0.2),colour)
    } else {
      cmp <- colour
    }
    cmp <- c(rgb(0,0,0,0),colorRampPalette(cmp)(100))
    el <- setdef(el,list(type="o",lty=1,pch=21,cex=1.1,lwd=2,
                         bg=RColorBrewer::brewer.pal(3,"Set1")[1]))
    args$xaxt <- "n"
    do.call(image,c(list(x=1:m,y=llc,z=dnst,col=cmp),args))
    do.call(lines,c(list(x=1:m,y=colMeans(ynt,na.rm=TRUE),col="black"),el))
    refline()
    box()
    xaxis(1:m)
  }

  invisible(NULL)

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
  if (!is.null(x$season.exist)){
    cat(paste0("Evidence of seasonality: ",x$season.exist, "  (pval: ",round(x$season.pval,3),")\n"))
  } else {
    cat("Presence of seasonality not tested (not enough complete seasons).\n")
  }
  invisible(x)
}
