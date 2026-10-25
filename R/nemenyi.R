#' Nonparametric multiple comparisons (Nemenyi test)
#'
#' Perform nonparametric multiple comparisons, across columns, using the Friedman and the post-hoc Nemenyi tests.
#'
#' @param data an array that includes values to be compared for several treatments (in columns) for several observations (rows), of size n x k. For example, if these are forecast errors, different methods should be in columns and errors for different time series or forecast origins in rows.
#' @param conf.level the confidence level used for the comparison. Default is 0.95.
#' @param sort if \code{TRUE}, then function sorts the outputted values of mean ranks. If plots are request, this is forced to \code{TRUE}.
#' @param plottype type of plot to produce:
#' \itemize{
#'  \item \code{"none"}: no plot.
#'  \item \code{"mcb"}: \emph{Multiple Comparison with the Best} style plot.
#'  \item \code{"vmcb"}: vertical \emph{MCB} plot.
#'  \item \code{"line"}: summarised \emph{line} plot. Each line connects a group of treatments that are all not significantly different.
#'  \item \code{"vline"}: vertical \emph{line} plot.
#'  \item \code{"matrix"}: \emph{matrix} visualisation. Both axes are in mean rank order and each row/column is a tested treatment (dark cell on the diagonal). Each off-diagonal cell gives the comparison of a pair of treatments; shaded cells indicate no evidence of significant difference. As the comparison is symmetric, the upper triangle mirrors the lower, so each row can be read on its own as the outcome for the tested treatment.
#'  \item \code{"lmatrix"}: lower triangular \emph{matrix} visualisation. Both axes are in mean rank order, each row/column is a tested treatment (dark cell on the diagonal), and only the lower triangle of the matrix is shown, giving the comparison for each pair once. Shaded cells indicate no evidence of significant difference.
#'  \item \code{"ranges"}: \emph{anchored} plot. One row per tested treatment on a common mean rank axis. The shaded band is the tested treatment's mean rank +/- the critical distance; treatments inside the band (filled points) are not significantly different from it, treatments outside (open points) are.
#'  }
#' @param select highlight selected treatment (column). Number 1 to k. Use NULL for no highlighting.
#' @param labels optional labels for models. If NULL column names of \code{data} will be used.
#' @param group if \code{TRUE}, show the group membership with a compact letter display: treatments that share a letter are not significantly different. For \code{"ranges"}, \code{"vline"} and \code{"matrix"} the letters are shown on the right axis, for \code{"line"} on the top axis, for \code{"lmatrix"} in the empty upper triangle, next to each diagonal cell, and for \code{"mcb"} and \code{"vmcb"} they are appended to the labels. If \code{NULL} (default), they are shown for all plots except \code{"mcb"} and \code{"vmcb"}.
#' @param line margin lines (as in the \code{line} argument of \code{title} and \code{mtext}) where the plot annotation is placed. A named vector with elements:
#' \itemize{
#'  \item \code{main}: plot title.
#'  \item \code{xlab}: x-axis title.
#'  \item \code{ylab}: y-axis title.
#'  \item \code{xaxis}: x-axis labels (treatment labels or numbers).
#'  \item \code{yaxis}: y-axis labels (treatment labels or numbers).
#'  \item \code{group}: compact letter display, see \code{group}. For \code{"lmatrix"} this is the distance from the right edge of the diagonal cell.
#'  }
#' Elements that are \code{NA} are placed automatically: axis labels and letters at \code{par("mgp")[2]} (letters at 0.3 for \code{"lmatrix"}), axis titles outside the axis labels and the title above the letters when these are on the top axis (otherwise in its default position). Only the elements to change need to be given, e.g. \code{line=c(main=0.5, ylab=10)}; an unnamed vector is assigned in the order above. Margins are enlarged to fit the requested lines.
#' @param ... additional arguments passed to the \code{plot} function. Axis annotation arguments (\code{cex.axis}, \code{font.axis}, \code{col.axis}, \code{cex.lab}, \code{font.lab}, \code{col.lab}, \code{cex.main}, \code{font.main}, \code{col.main}, \code{mgp}) are also used for the axes and titles the function draws.
#'
#' @details
#' Plots respect the graphical parameters currently set with \code{par()} (for example \code{mfrow}, \code{cex}, \code{mgp}, \code{fg}, \code{bg}, \code{family} or \code{las} for numerical axes). The only parameter changed is \code{mar}: margins that are too small for the treatment labels, the group letters or the lines requested in \code{line} are enlarged, using the actual width of the labels, while all other margins are kept as set by the user. The user's \code{mar} is restored when the function exits. As the plot region depends on the margins, low-level additions to the plot after the call (e.g. \code{abline}) may not align; set \code{par(mar=...)} large enough for the labels beforehand to avoid this.
#'
#' @return Return object of class \code{nemenyi} and contains:
#' \itemize{
#' \item \code{means}: mean rank of each treatment.
#' \item \code{intervals}: intervals within there is no evidence of significance difference according to the Nemenyi test at requested confidence level.
#' \item \code{fpval}: Friedman test p-value.
#' \item \code{fH}: Friedman test hypothesis outcome.
#' \item \code{cd}: Nemenyi critical distance. Output \code{intervals} is calculate as \code{means} +/- \code{cd}.
#' \item \code{conf.level}: confidence level used for testing.
#' \item \code{k}: number of treatments (columns).
#' \item \code{n}: number of observations (rows).
#' \item \code{group}: group membership as a compact letter display; treatments sharing a letter are not significantly different.
#' \item \code{diff}: k x k matrix of absolute differences in mean ranks. Pairs with a difference smaller than \code{cd} are not significantly different.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com},
#' @author Ivan Svetunkov, \email{ivan@svetunkov.ru}.
#'
#' @references
#' \itemize{
#' \item The tests are detailed by Hollander, M., Wolfe, D.A. and Chicken, E. (2014) Nonparametric Statistical Methods. 3rd Edition, John Wiley & Sons, Inc., New York.
#' \item The \emph{line} plot is introduced \href{https://kourentzes.com/forecasting/wp-content/uploads/2014/04/ISF2012_Tests_Kourentzes.pdf}{here} and a first example of its use, along with a short description is provided by Kourentzes, N. (2013). \href{https://kourentzes.com/forecasting/2013/04/19/intermittent-demand-forecasts-with-neural-networks/}{Intermittent demand forecasts with neural networks}. International Journal of Production Economics, 143(1), 198-206.
#' \item The \emph{matrix} plot is introduced by Kourentzes, N., & Athanasopoulos, G. (2018). Cross-temporal coherent forecasts for Australian tourism (No. 24/18). Monash University, Department of Econometrics and Business Statistics.
#' \item The \emph{MCB} plot is described by Koning, A. J., Franses, P. H., Hibon, M., & Stekler, H. O. (2005). The M3 competition: Statistical tests of the results. International Journal of Forecasting, 21(3), 397-409.
#' }
#'
#' @keywords htest
#'
#' @examples
#' x <- matrix( rnorm(50*4,mean=0,sd=1), 50, 4)
#' x[,2] <- x[,2]+1
#' x[,3] <- x[,3]+0.7
#' x[,4] <- x[,4]+0.5
#' colnames(x) <- c("Method A","Method B","Method C - long name","Method D")
#' nemenyi(x,conf.level=0.95,plottype="vline")
#' nemenyi(x,plottype="matrix")
#' nemenyi(x,plottype="lmatrix")
#' nemenyi(x,plottype="ranges")
#' nemenyi(x,plottype="line",ylab="Groups",line=c(main=0.5,group=0.2))
#'
#' @export nemenyi

nemenyi <- function(data, conf.level=0.95, sort=c(TRUE,FALSE),
                    plottype=c("vline","none","mcb","vmcb","line","matrix","lmatrix","ranges"),
                    select=NULL, labels=NULL, group=NULL,
                    line=c("main"=NA,"xlab"=NA,"ylab"=NA,"xaxis"=NA,"yaxis"=NA,"group"=NA), ...){

  # Default
  sort <- sort[1]
  plottype <- match.arg(plottype,c("vline","none","mcb","vmcb","line","matrix","lmatrix","ranges"))

  # Check data
  if (length(dim(data)) != 2){
    stop("Data must be organised as methods in columns and observations in rows.")
  }
  data <- as.matrix(data)
  data <- na.exclude(data)
  rows.number <- nrow(data)
  cols.number <- ncol(data)
  if (cols.number < 2){
    stop("At least two treatments (columns) are needed.")
  }

  # Check select argument
  if (!is.null(select) && ((select < 1) || (select > cols.number))){
    select <- NULL
  }

  # If plot is asked, always sort the results
  if (plottype != "none"){
    sort <- TRUE
  }

  # Group letters shown by default for all plots, except MCB where they are appended to the labels
  if (is.null(group)){
    group <- !(plottype %in% c("mcb","vmcb"))
  }

  # Lines for plot annotation
  line <- nemenyi.line(line)

  # Checks for labels
  if (is.null(labels)){
    labels <- colnames(data)
    if (is.null(labels)){
      labels <- 1:cols.number
    }
  } else {
    labels <- labels[1:cols.number]
  }
  labels <- as.character(labels)

  # First run Friedman test. If insignificant then ignore Nemenyi (Weaker)
  fried.pval <- stats::friedman.test(data)$p.value
  if (fried.pval <= 1-conf.level){
    fried.H <- "Ha: Different" # At least one method is different
  } else {
    fried.H <- "H0: Identical" # No evidence of differences between methods
  }
  fried.sig <- fried.pval <= 1-conf.level

  # Nemenyi critical distance and bounds of intervals
  # Nikos: Use the formulas for mean ranks, not sums
  r.stat <- stats::qtukey(conf.level,cols.number,Inf)*sqrt((cols.number*(cols.number+1))/(12*rows.number)) # Means
  # r.stat <- qtukey(conf.level,cols.number,Inf)*sqrt((rows.number*cols.number*(cols.number+1))/(12)) # Sums
  # NSM3::cWNMT(0.95, cols.number, rows.number, method="Asymptotic") # This is the sums result according to Hollander, Wolfe & Chicken, 2014, p. 332 (7.27)

  # Rank methods for each time series
  ranks.matrix <- t(apply(data,1,function(x){rank(x,na.last="keep",ties.method="average")}))

  # Calculate mean rank values
  ranks.means <- colMeans(ranks.matrix)
  names(ranks.means) <- labels

  # Sort means and labels
  if(sort==TRUE){
    order.idx <- order(ranks.means)
  } else {
    order.idx <- 1:cols.number
  }
  ranks.means <- ranks.means[order.idx]
  labels <- labels[order.idx]
  if (!is.null(select)){
    select <- which(order.idx == select)
  }

  # Calculate intervals for each of the methods
  # The comparison is abs(diff(ranks)) < r.stat, otherwise different groups
  ranks.intervals <- rbind(ranks.means - r.stat,ranks.means + r.stat)

  # Pairwise absolute differences and outcome
  ranks.diff <- abs(outer(ranks.means,ranks.means,"-"))
  dimnames(ranks.diff) <- list(labels,labels)
  ranks.same <- ranks.diff < r.stat

  # Compact letter display
  ranks.cld <- nemenyi.cld(ranks.means,r.stat)

  # Produce plots
  # For all plots
  if (plottype != "none"){
    args <- list(...)

    # Create title for plots. It is drawn separately, so that its line can be controlled
    if (is.null(args[["main"]])){
      main <- paste0("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ")\nCritical distance: ", format(round(r.stat,3),nsmall=3))
    } else {
      main <- args[["main"]]
    }
    args$main <- NULL

    # Remaining defaults
    if (is.null(args[["xaxs"]])){
      args$xaxs <- "i"
    }
    if (is.null(args[["yaxs"]])){
      args$yaxs <- "i"
    }

    # Axis titles are drawn separately, so that they can be placed outside the labels
    xlab <- if (is.null(args[["xlab"]])) "" else args[["xlab"]]
    ylab <- if (is.null(args[["ylab"]])) "" else args[["ylab"]]
    args$xlab <- args$ylab <- ""
    args$x <- args$y <- NA
    args$axes <- FALSE

    # Axis annotation settings: arguments in ... take precedence over par()
    ax <- nemenyi.axpar(args)

    # Keep user margins, all other graphical parameters are not modified
    mar.user <- graphics::par("mar")
    on.exit(graphics::par(mar=mar.user), add=TRUE)

    # Colours. Light fills are opaque, mixed with the background, so that they
    # are kept in devices that do not support semi-transparency (e.g. postscript)
    fg <- graphics::par("fg")
    col.best <- RColorBrewer::brewer.pal(3,"Set1")[1]
    col.ns <- RColorBrewer::brewer.pal(3,"Set1")[2]

    # Prepare method labels, with letters for MCB plots if requested, and add mean rank to them
    # For all other plots the letters are shown on a separate axis
    if (group && (plottype %in% c("mcb","vmcb"))){
      lbl <- paste0(labels," [",ranks.cld,"]")
      gtxt <- NULL
    } else {
      lbl <- labels
      gtxt <- if (group) ranks.cld else NULL
    }
    lbl.rank <- paste0(lbl," - ",sprintf("%1.2f",round(ranks.means,2)))
  }

  # MCB style plots
  if ((plottype == "mcb")|(plottype == "vmcb")){

    # Set colours
    if (fried.sig){pcol <- col.ns} else {pcol <- "gray"}

    # Find min max
    mnmx <- range(ranks.means) + c(-0.5,0.5)*r.stat
    mnmx <- mnmx + diff(mnmx)*0.04*c(-1,1)

    # Set plot - horizontal or vertical
    if (plottype == "mcb"){ # Horizontal
      if (is.null(list(...)[["ylab"]])){
        ylab <- "Mean ranks"
      }
      if(is.null(args[["xlim"]])){
        args$xlim <- c(0,cols.number+1)
      }
      if(is.null(args[["ylim"]])){
        args$ylim <- mnmx
      }
      lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,xtxt=lbl.rank)
    } else { # Vertical
      if (is.null(list(...)[["xlab"]])){
        xlab <- "Mean ranks"
      }
      if(is.null(args[["ylim"]])){
        args$ylim <- c(0,cols.number+1)
      }
      if(is.null(args[["xlim"]])){
        args$xlim <- mnmx
      }
      lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,ytxt=lbl.rank)
    }
    graphics::par(mar=lay$mar)

    if (is.null(select)){
      sel <- 1
    } else {
      sel <- select
    }

    # Use do.call to use manipulated ellipsis (...)
    do.call(graphics::plot,args)
    # Plot rest
    if (plottype == "mcb"){
      # Intervals for best method
      graphics::polygon(c(0,rep(cols.number+1,2),0),rep(ranks.means[sel],4)+r.stat/2*c(1,1,-1,-1),col="gray90",border=NA)
      # Ranks
      graphics::points(1:cols.number,ranks.means,pch=20,lwd=3,col=fg)
      nemenyi.axis(1,lay,ax,at=1:cols.number,labels=lbl.rank,las=2)
      nemenyi.axis(2,lay,ax)
      # Intervals for all methods
      for (i in 1:cols.number){
        graphics::lines(rep(i,times=2),ranks.means[i]+c(-1,1)*0.5*r.stat, type="o", lwd=1, col=pcol,pch=20)
      }
      # Highlight identical
      idx <- ranks.same[sel,]
      graphics::points((1:cols.number)[idx],ranks.means[idx],pch=20,lwd=3,col=col.best)
    } else { # vmcb
      # Intervals for best method
      graphics::polygon(rep(ranks.means[sel],4)+r.stat/2*c(1,1,-1,-1),c(0,rep(cols.number+1,2),0),col="gray90",border=NA)
      # Ranks
      graphics::points(ranks.means,1:cols.number,pch=20,lwd=3,col=fg)
      nemenyi.axis(2,lay,ax,at=1:cols.number,labels=lbl.rank,las=2)
      nemenyi.axis(1,lay,ax)
      # Intervals for all methods
      for (i in 1:cols.number){
        graphics::lines(ranks.means[i]+c(-1,1)*0.5*r.stat, rep(i,times=2), type="o", lwd=1, col=pcol,pch=20)
      }
      # Highlight identical
      idx <- ranks.same[sel,]
      graphics::points(ranks.means[idx],(1:cols.number)[idx],pch=20,lwd=3,col=col.best)
    }
    graphics::box(which="plot")
    nemenyi.titles(xlab,ylab,main,lay,ax)

  }

  # Matrix visualisations, full or lower triangle
  if ((plottype == "matrix")|(plottype == "lmatrix")){

    # Colours
    if (fried.sig){
      fill.tested <- fg
      fill.ns <- colour.blend(col.ns,0.5)
    } else {
      fill.tested <- "gray50"
      fill.ns <- colour.blend("gray50",0.4)
    }

    # Margins for both label axes, and letters on the right
    # For lmatrix the letters are next to the diagonal, so only the last one
    # (bottom row) reaches the right margin
    if (plottype == "lmatrix"){
      if (is.na(line["group"])){
        line["group"] <- 0.3
      }
      if (!is.null(gtxt)){
        gtxt <- gtxt[cols.number]
      }
    }
    lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,xtxt=lbl,ytxt=lbl.rank,gtxt=gtxt,gside=4)
    graphics::par(mar=lay$mar)

    if(is.null(args[["xlim"]])){
      args$xlim <- c(0.5,cols.number+0.5)
    }
    if(is.null(args[["ylim"]])){
      args$ylim <- c(0.5,cols.number+0.5)
    }
    do.call(graphics::plot,args)

    # Size of numbers in cells, so that they fit in a cell
    cex.cell <- min(1, 0.8/graphics::strwidth("0.00"), 0.6/graphics::strheight("0.00"))

    # Best method at the top. Full matrix, or lower triangle only
    # Cells are not clipped, so that borders on the edge of the plot region have the
    # same width as the inner ones. Axes are drawn without a line, only ticks, so that
    # they do not thicken the outer borders
    for (i in 1:cols.number){
      y <- cols.number - i + 1
      cols.plot <- if (plottype == "matrix") 1:cols.number else 1:i
      for (j in cols.plot){
        if (i == j){
          fcol <- if (!is.null(select) && (select == i)) col.best else fill.tested
          graphics::rect(j-0.5,y-0.5,j+0.5,y+0.5,col=fcol,border=fg,xpd=NA)
        } else {
          fcol <- if (ranks.same[i,j]) fill.ns else NA
          graphics::rect(j-0.5,y-0.5,j+0.5,y+0.5,col=fcol,border=fg,xpd=NA)
          # Add rank differences
          # Makes the figure cluttered. I need to figure out a better way to show this
          # if (cex.cell >= 0.4){
          #   graphics::text(j,y,sprintf("%1.2f",ranks.diff[i,j]),cex=cex.cell,col=fg)
          # }
        }
      }
    }

    nemenyi.axis(1,lay,ax,at=1:cols.number,labels=lbl,las=2,lwd=0,lwd.ticks=graphics::par("lwd"))
    nemenyi.axis(2,lay,ax,at=cols.number:1,labels=lbl.rank,las=2,lwd=0,lwd.ticks=graphics::par("lwd"))
    if (group){
      if (plottype == "matrix"){
        nemenyi.axis(4,lay,ax,at=cols.number:1,labels=ranks.cld,las=2,tick=FALSE)
      } else {
        # Letters in the empty upper triangle, right of each diagonal cell
        # Convert margin lines to user x units
        usr <- graphics::par("usr")
        dx <- lay$line[["group"]]*graphics::par("csi")*graphics::par("mex") /
          (graphics::par("pin")[1]/diff(usr[1:2]))
        graphics::text((1:cols.number)+0.5+dx, cols.number:1, ranks.cld, adj=c(0,0.5), xpd=NA,
                       cex=ax$cex.axis, font=ax$font.axis, col=ax$col.axis, family=ax$family)
      }
    }
    nemenyi.titles(xlab,ylab,main,lay,ax)

  }

  # Anchored ranges plot
  if (plottype == "ranges"){

    # Colours
    if (fried.sig){
      band.col <- colour.blend(col.ns,0.25)
      pcol <- col.ns
    } else {
      band.col <- colour.blend("gray50",0.25)
      pcol <- "gray50"
    }

    if (is.null(list(...)[["xlab"]])){
      xlab <- "Mean ranks"
    }

    # Margins for labels, and letters on the right
    lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,ytxt=lbl.rank,gtxt=gtxt,gside=4)
    graphics::par(mar=lay$mar)

    if(is.null(args[["xlim"]])){
      mnmx <- range(ranks.means) + c(-1,1)*r.stat
      args$xlim <- mnmx + diff(mnmx)*0.04*c(-1,1)
    }
    if(is.null(args[["ylim"]])){
      args$ylim <- c(0.5,cols.number+0.5)
    }
    do.call(graphics::plot,args)

    # Bands of no significant difference from the tested method
    for (i in 1:cols.number){
      y <- cols.number - i + 1
      graphics::rect(ranks.means[i]-r.stat,y-0.3,ranks.means[i]+r.stat,y+0.3,col=band.col,border=NA)
    }

    # Guides at the mean rank of each method, drawn over the opaque bands
    graphics::abline(v=ranks.means,col="gray80",lty=3)

    for (i in 1:cols.number){
      y <- cols.number - i + 1
      same <- ranks.same[i,]
      same[i] <- FALSE
      diffr <- !ranks.same[i,]
      if (any(diffr)){
        graphics::points(ranks.means[diffr],rep(y,sum(diffr)),pch=1,col="gray50")
      }
      if (any(same)){
        graphics::points(ranks.means[same],rep(y,sum(same)),pch=19,col=pcol)
      }
      graphics::points(ranks.means[i],y,pch=19,cex=1.4,
                       col=if (!is.null(select) && (select == i)) col.best else fg)
    }

    nemenyi.axis(1,lay,ax)
    nemenyi.axis(2,lay,ax,at=cols.number:1,labels=lbl.rank,las=2)
    if (group){
      nemenyi.axis(4,lay,ax,at=cols.number:1,labels=ranks.cld,las=2,tick=FALSE)
    }
    graphics::box(which="plot")
    nemenyi.titles(xlab,ylab,main,lay,ax)

  }

  # Line style plot (as in ISF reference)
  if ((plottype == "line")|(plottype == "vline")){

    # Find groups: each line connects treatments that are all pairwise not
    # significantly different (same groups as the compact letter display)
    rline <- nemenyi.groups(ranks.means,r.stat)
    # Remove single member groups
    rline <- rline[rline[,1] != rline[,2],,drop=FALSE]
    k <- nrow(rline)
    # Choose colour depending on Friedman test result
    cmp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(k)
    if (!fried.sig){cmp <- rep("gray",times=k)}

    if (plottype == "line"){
      if(is.null(args[["xlim"]])){
        args$xlim <- c(1,cols.number)
      }
      if(is.null(args[["ylim"]])){
        args$ylim <- c(0,k+1)
      }
      # Labels at the bottom, letters at the top
      lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,xtxt=lbl.rank,gtxt=gtxt,gside=3)
    } else { # vline
      if(is.null(args[["xlim"]])){
        args$xlim <- c(0,k+1)
      }
      if(is.null(args[["ylim"]])){
        args$ylim <- c(1,cols.number)
      }
      # Labels on the left, letters on the right
      lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,ytxt=lbl.rank,gtxt=gtxt,gside=4)
    }
    graphics::par(mar=lay$mar)

    do.call(graphics::plot,args)

    if (plottype == "line"){
      graphics::points(1:cols.number,rep(0,cols.number),pch=20,lwd=4,col=fg)
      if (k>0){
        for (i in 1:k){
          graphics::lines(rline[i,],c(i,i), col=cmp[i], lwd = 4)
          graphics::lines(rep(rline[i,1],times=2),c(0,i), col="gray", lty = 2)
          graphics::lines(rep(rline[i,2],times=2),c(0,i), col="gray", lty = 2)
        }
      }
      nemenyi.axis(1,lay,ax,at=1:cols.number,labels=lbl.rank,las=2)
      if (group){
        nemenyi.axis(3,lay,ax,at=1:cols.number,labels=ranks.cld,las=1,tick=FALSE)
      }
      if (!is.null(select)){
        graphics::points(select,0,pch=20,col=col.best,cex=2)
      }
    } else { # vline
      graphics::points(rep(0,cols.number),1:cols.number,pch=20,lwd=4,col=fg)
      if (k>0){
        for (i in 1:k){
          graphics::lines(c(i,i), rline[i,], col=cmp[i], lwd = 4)
          graphics::lines(c(0,i), rep(rline[i,1],times=2), col="gray", lty = 2)
          graphics::lines(c(0,i), rep(rline[i,2],times=2), col="gray", lty = 2)
        }
      }
      nemenyi.axis(2,lay,ax,at=1:cols.number,labels=lbl.rank,las=2)
      if (group){
        nemenyi.axis(4,lay,ax,at=1:cols.number,labels=ranks.cld,las=2,tick=FALSE)
      }
      if (!is.null(select)){
        graphics::points(0,select,pch=20,col=col.best,cex=2)
      }
    }
    nemenyi.titles(xlab,ylab,main,lay,ax)

  }

  return(structure(list("means"=ranks.means,"intervals"=ranks.intervals,"fpval"=fried.pval,"fH"=fried.H,
                        "cd"=r.stat,"conf.level"=conf.level,"k"=cols.number,"n"=rows.number,
                        "group"=ranks.cld,"diff"=ranks.diff),class="nemenyi"))

}

# Groups of treatments that are not significantly different from each other
# As the test is on distances of mean ranks, each group is a run of consecutive
# treatments in rank order. Only maximal groups are kept: every pair within a group
# is not significantly different. Returns start and end positions in rank order.
nemenyi.groups <- function(means, cd){
  m <- sort(means)
  k <- length(m)
  # Last treatment (in rank order) that is not different from treatment i
  hi <- sapply(1:k, function(i){max(which(m - m[i] < cd))})
  # Keep only groups that are not contained in the previous one
  keep <- c(TRUE, hi[-1] > hi[-k])
  cbind("start"=which(keep), "end"=hi[keep])
}

# Compact letter display for the Nemenyi test
# Treatments sharing a letter are not significantly different
nemenyi.cld <- function(means, cd){
  k <- length(means)
  idx <- order(means)
  grp <- nemenyi.groups(means, cd)
  ng <- nrow(grp)
  if (ng <= 52){
    let <- c(letters,LETTERS)[1:ng]
  } else {
    let <- paste0(".",1:ng)
  }
  out <- rep("",k)
  for (g in 1:ng){
    out[grp[g,1]:grp[g,2]] <- paste0(out[grp[g,1]:grp[g,2]],let[g])
  }
  res <- rep("",k)
  res[idx] <- out
  names(res) <- names(means)
  return(res)
}

# Lines of plot annotation. Accepts a named (in any order, partial names allowed)
# or unnamed (in the default order) vector or list. Missing elements are NA (automatic)
nemenyi.line <- function(line){
  nms <- c("main","xlab","ylab","xaxis","yaxis","group")
  out <- stats::setNames(rep(NA_real_,length(nms)),nms)
  if (is.null(line)){
    return(out)
  }
  line <- unlist(line)
  if (length(line) > length(nms)){
    stop(paste0("Argument line can have up to ", length(nms), " elements: ", paste(nms,collapse=", "), "."))
  }
  nm <- names(line)
  if (is.null(nm) || all(nm == "")){
    out[seq_along(line)] <- as.numeric(line)
  } else {
    if (any(nm == "")){
      stop("Elements of argument line must be either all named or all unnamed.")
    }
    idx <- pmatch(nm, nms)
    if (any(is.na(idx))){
      stop(paste0("Unknown element(s) in argument line: ", paste(nm[is.na(idx)],collapse=", "),
                  ". Use: ", paste(nms,collapse=", "), "."))
    }
    out[idx] <- as.numeric(line)
  }
  out
}

# Axis annotation settings, arguments passed by the user take precedence over par()
nemenyi.axpar <- function(args){
  getp <- function(nm){
    if (!is.null(args[[nm]])){
      args[[nm]]
    } else {
      graphics::par(nm)
    }
  }
  list(cex.axis=getp("cex.axis"), font.axis=getp("font.axis"), col.axis=getp("col.axis"),
       cex.lab=getp("cex.lab"), font.lab=getp("font.lab"), col.lab=getp("col.lab"),
       cex.main=getp("cex.main"), font.main=getp("font.main"), col.main=getp("col.main"),
       mgp=getp("mgp"), family=getp("family"))
}

# Width of text drawn perpendicular to an axis, in margin lines
nemenyi.textlines <- function(txt, ax){
  w <- max(graphics::strwidth(txt, units="inches", cex=ax$cex.axis,
                              font=ax$font.axis, family=ax$family))
  w / (graphics::par("csi")*graphics::par("mex"))
}

# Is there a title to draw, and how many lines does it span
nemenyi.nlines <- function(txt){
  if (is.null(txt)){
    0
  } else if (is.character(txt)){
    length(unlist(strsplit(txt[nzchar(txt)],"\n")))
  } else {
    1
  }
}

# Resolve the margin line of each annotation element and enlarge the margins
# that cannot fit them. Other margins are kept as set by the user.
#   xtxt, ytxt: treatment labels drawn perpendicular on sides 1 and 2, NULL for numerical axes
#   gtxt, gside: group letters and the side they are drawn on (3 parallel, 4 perpendicular)
nemenyi.layout <- function(mar, line, ax, xlab, ylab, main,
                           xtxt=NULL, ytxt=NULL, gtxt=NULL, gside=4){
  pos <- line
  # Axis labels and letters default to the usual axis label line
  for (nm in c("xaxis","yaxis","group")){
    if (is.na(pos[nm])){
      pos[nm] <- ax$mgp[2]
    }
  }

  # Sides 1 and 2: axis labels, and axis title outside them
  sides <- list(list(side=1, axis="xaxis", lab="xlab", txt=xtxt, title=xlab),
                list(side=2, axis="yaxis", lab="ylab", txt=ytxt, title=ylab))
  for (s in sides){
    if (!is.null(s$txt)){
      # Treatment labels, perpendicular to the axis
      need <- pos[[s$axis]] + nemenyi.textlines(s$txt,ax) + 0.3
      auto <- need
      gap <- 0.3
    } else {
      # Numerical axis, keep the usual distance between labels and title
      need <- pos[[s$axis]] + ax$cex.axis
      auto <- pos[[s$axis]] + ax$mgp[1] - ax$mgp[2]
      gap <- 0
    }
    if (nemenyi.nlines(s$title) > 0){
      if (is.na(pos[s$lab])){
        pos[s$lab] <- auto
      }
      need <- max(need, pos[[s$lab]] + ax$cex.lab + gap)
    }
    mar[s$side] <- max(mar[s$side], need)
  }

  # Group letters
  top <- NA
  if (!is.null(gtxt)){
    if (gside == 3){
      need <- pos[["group"]] + ax$cex.axis + 0.3
      top <- need
    } else {
      need <- pos[["group"]] + nemenyi.textlines(gtxt,ax) + 0.3
    }
    mar[gside] <- max(mar[gside], need)
  }

  # Main title: above the letters if these are on top, otherwise default position
  nmain <- nemenyi.nlines(main)
  if (nmain > 0){
    if (is.na(pos["main"])){
      pos["main"] <- top
    }
    if (!is.na(pos["main"])){
      mar[3] <- max(mar[3], pos[["main"]] + nmain*ax$cex.main + 0.3)
    }
  }

  list(mar=mar, line=pos)
}

# Draw an axis with the user's annotation settings, at the resolved line
nemenyi.axis <- function(side, lay, ax, ...){
  pos <- lay$line[[c("xaxis","yaxis","group","group")[side]]]
  graphics::axis(side, cex.axis=ax$cex.axis, font.axis=ax$font.axis,
                 col.axis=ax$col.axis, mgp=c(ax$mgp[1],pos,ax$mgp[3]), family=ax$family, ...)
}

# Draw main and axis titles, at the resolved lines
nemenyi.titles <- function(xlab, ylab, main, lay, ax){
  pos <- lay$line
  if (nemenyi.nlines(xlab) > 0){
    graphics::title(xlab=xlab, line=pos[["xlab"]], cex.lab=ax$cex.lab,
                    font.lab=ax$font.lab, col.lab=ax$col.lab, mgp=ax$mgp, family=ax$family)
  }
  if (nemenyi.nlines(ylab) > 0){
    graphics::title(ylab=ylab, line=pos[["ylab"]], cex.lab=ax$cex.lab,
                    font.lab=ax$font.lab, col.lab=ax$col.lab, mgp=ax$mgp, family=ax$family)
  }
  if (nemenyi.nlines(main) > 0){
    graphics::title(main=main, line=pos[["main"]], cex.main=ax$cex.main,
                    font.main=ax$font.main, col.main=ax$col.main, family=ax$family)
  }
}

#' @export
#' @method summary nemenyi
summary.nemenyi <- function(object,...){
  print(object)
}

#' @export
#' @method print nemenyi
print.nemenyi <- function(x,...){

  writeLines("Friedman and Nemenyi Tests")
  writeLines(paste0("The confidence level is ", x$conf.level*100, "% (significance level ", (1-x$conf.level)*100, "%)"))
  writeLines(paste0("Number of observations is ", x$n, " and number of methods is ", x$k))
  writeLines(paste0("Friedman test p-value: ", format(round(x$fpval,4),nsmall=4) , " - ", x$fH))
  writeLines(paste0("Critical distance: ", format(round(x$cd,4),nsmall=4)))
  if (!is.null(x$group)){
    writeLines("")
    out <- data.frame("Mean rank"=format(round(x$means,4),nsmall=4), "Group"=x$group, check.names=FALSE)
    rownames(out) <- names(x$means)
    print(out)
  }

}
