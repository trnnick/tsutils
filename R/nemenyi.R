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
#'  \item \code{"matrix"}: \emph{matrix} visualisation. Both axes are in mean rank order, each row/column is a tested treatment (dark cell on the diagonal), and the lower triangle gives the absolute difference in mean ranks for each pair. Shaded cells indicate no evidence of significant difference.
#'  \item \code{"ranges"}: \emph{anchored} plot. One row per tested treatment on a common mean rank axis. The shaded band is the tested treatment's mean rank +/- the critical distance; treatments inside the band (filled points) are not significantly different from it, treatments outside (open points) are.
#'  }
#' @param select highlight selected treatment (column). Number 1 to k. Use NULL for no highlighting.
#' @param labels optional labels for models. If NULL column names of \code{data} will be used.
#' @param cld if \code{TRUE}, show the compact letter display: treatments that share a letter are not significantly different. For \code{"ranges"} the letters are shown on the right axis, for the other plots they are appended to the labels. If \code{NULL} (default), they are shown only for \code{"ranges"}.
#' @param ... additional arguments passed to the \code{plot} function. Axis annotation arguments (\code{cex.axis}, \code{font.axis}, \code{col.axis}, \code{cex.lab}, \code{font.lab}, \code{col.lab}, \code{mgp}) are also used for the axes the function draws.
#'
#' @details
#' Plots respect the graphical parameters currently set with \code{par()} (for example \code{mfrow}, \code{cex}, \code{mgp}, \code{fg}, \code{bg}, \code{family} or \code{las} for numerical axes). The only parameter changed is \code{mar}: margins that are too small for the treatment labels are enlarged, using the actual width of the labels, while all other margins are kept as set by the user. The user's \code{mar} is restored when the function exits. As the plot region depends on the margins, low-level additions to the plot after the call (e.g. \code{abline}) may not align; set \code{par(mar=...)} large enough for the labels beforehand to avoid this.
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
#' \item \code{cld}: compact letter display; treatments sharing a letter are not significantly different.
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
#' nemenyi(x,plottype="ranges")
#'
#' @export nemenyi

nemenyi <- function(data, conf.level=0.95, sort=c(TRUE,FALSE),
                    plottype=c("vline","none","mcb","vmcb","line","matrix","ranges"),
                    select=NULL, labels=NULL, cld=NULL, ...){

  # Default
  sort <- sort[1]
  plottype <- match.arg(plottype,c("vline","none","mcb","vmcb","line","matrix","ranges"))

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

  # Compact letter display shown by default only for ranges plot
  if (is.null(cld)){
    cld <- plottype == "ranges"
  }

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

    # Create title for plots
    if (is.null(args[["main"]])){
      args$main <- paste0("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ")\nCritical distance: ", format(round(r.stat,3),nsmall=3))
    }

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

    # Colours (semi-transparent fills work on any background)
    fg <- graphics::par("fg")
    col.best <- RColorBrewer::brewer.pal(3,"Set1")[1]
    col.ns <- RColorBrewer::brewer.pal(3,"Set1")[2]

    # Prepare method labels, with letters if requested, and add mean rank to them
    if (cld && plottype != "ranges"){
      lbl <- paste0(labels," [",ranks.cld,"]")
    } else {
      lbl <- labels
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
      mm <- nemenyi.margins(mar.user,1,list(lbl.rank),list(xlab),ax)
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
      mm <- nemenyi.margins(mar.user,2,list(lbl.rank),list(ylab),ax)
    }
    graphics::par(mar=mm$mar)

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
      nemenyi.axis(1,at=1:cols.number,labels=lbl.rank,ax=ax,las=2)
      nemenyi.axis(2,ax=ax)
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
      nemenyi.axis(2,at=1:cols.number,labels=lbl.rank,ax=ax,las=2)
      nemenyi.axis(1,ax=ax)
      # Intervals for all methods
      for (i in 1:cols.number){
        graphics::lines(ranks.means[i]+c(-1,1)*0.5*r.stat, rep(i,times=2), type="o", lwd=1, col=pcol,pch=20)
      }
      # Highlight identical
      idx <- ranks.same[sel,]
      graphics::points(ranks.means[idx],(1:cols.number)[idx],pch=20,lwd=3,col=col.best)
    }
    graphics::box(which="plot")
    nemenyi.titles(xlab,ylab,mm$title.line,ax)

  }

  # Matrix visualisation
  if (plottype == "matrix"){

    # Colours
    if (fried.sig){
      fill.tested <- fg
      fill.ns <- grDevices::adjustcolor(col.ns,alpha.f=0.5)
    } else {
      fill.tested <- "gray50"
      fill.ns <- grDevices::adjustcolor("gray50",alpha.f=0.4)
    }

    # Margins for both label axes
    mm <- nemenyi.margins(mar.user,c(1,2),list(lbl,lbl.rank),list(xlab,ylab),ax)
    graphics::par(mar=mm$mar)

    if(is.null(args[["xlim"]])){
      args$xlim <- c(0.5,cols.number+0.5)
    }
    if(is.null(args[["ylim"]])){
      args$ylim <- c(0.5,cols.number+0.5)
    }
    do.call(graphics::plot,args)

    # Size of numbers in cells, so that they fit in a cell
    cex.cell <- min(1, 0.8/graphics::strwidth("0.00"), 0.6/graphics::strheight("0.00"))

    # Best method at the top, lower triangle only
    for (i in 1:cols.number){
      y <- cols.number - i + 1
      for (j in 1:i){
        if (i == j){
          fcol <- if (!is.null(select) && (select == i)) col.best else fill.tested
          graphics::rect(j-0.5,y-0.5,j+0.5,y+0.5,col=fcol,border=fg)
        } else {
          fcol <- if (ranks.same[i,j]) fill.ns else NA
          graphics::rect(j-0.5,y-0.5,j+0.5,y+0.5,col=fcol,border=fg)
          # Add rank differences
          # Makes the figure cluttered. I need to figure out a better way to show this
          # if (cex.cell >= 0.4){
          #   graphics::text(j,y,sprintf("%1.2f",ranks.diff[i,j]),cex=cex.cell,col=fg)
          # }
        }
      }
    }

    nemenyi.axis(1,at=1:cols.number,labels=lbl,ax=ax,las=2)
    nemenyi.axis(2,at=cols.number:1,labels=lbl.rank,ax=ax,las=2)
    nemenyi.titles(xlab,ylab,mm$title.line,ax)

  }

  # Anchored ranges plot
  if (plottype == "ranges"){

    # Colours
    if (fried.sig){
      band.col <- grDevices::adjustcolor(col.ns,alpha.f=0.25)
      pcol <- col.ns
    } else {
      band.col <- grDevices::adjustcolor("gray50",alpha.f=0.25)
      pcol <- "gray50"
    }

    if (is.null(list(...)[["xlab"]])){
      xlab <- "Mean ranks"
    }

    # Margins for labels, and letters on the right
    if (cld){
      mm <- nemenyi.margins(mar.user,c(2,4),list(lbl.rank,ranks.cld),list(ylab,""),ax)
    } else {
      mm <- nemenyi.margins(mar.user,2,list(lbl.rank),list(ylab),ax)
    }
    graphics::par(mar=mm$mar)

    if(is.null(args[["xlim"]])){
      mnmx <- range(ranks.means) + c(-1,1)*r.stat
      args$xlim <- mnmx + diff(mnmx)*0.04*c(-1,1)
    }
    if(is.null(args[["ylim"]])){
      args$ylim <- c(0.5,cols.number+0.5)
    }
    do.call(graphics::plot,args)

    # Guides at the mean rank of each method
    graphics::abline(v=ranks.means,col="gray80",lty=3)

    for (i in 1:cols.number){
      y <- cols.number - i + 1
      # Band of no significant difference from the tested method
      graphics::rect(ranks.means[i]-r.stat,y-0.3,ranks.means[i]+r.stat,y+0.3,col=band.col,border=NA)
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

    nemenyi.axis(1,ax=ax)
    nemenyi.axis(2,at=cols.number:1,labels=lbl.rank,ax=ax,las=2)
    if (cld){
      nemenyi.axis(4,at=cols.number:1,labels=ranks.cld,ax=ax,las=2,tick=FALSE)
    }
    graphics::box(which="plot")
    nemenyi.titles(xlab,ylab,mm$title.line,ax)

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
      mm <- nemenyi.margins(mar.user,1,list(lbl.rank),list(xlab),ax)
    } else { # vline
      if(is.null(args[["xlim"]])){
        args$xlim <- c(0,k+1)
      }
      if(is.null(args[["ylim"]])){
        args$ylim <- c(1,cols.number)
      }
      mm <- nemenyi.margins(mar.user,2,list(lbl.rank),list(ylab),ax)
    }
    graphics::par(mar=mm$mar)

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
      nemenyi.axis(1,at=1:cols.number,labels=lbl.rank,ax=ax,las=2)
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
      nemenyi.axis(2,at=1:cols.number,labels=lbl.rank,ax=ax,las=2)
      if (!is.null(select)){
        graphics::points(0,select,pch=20,col=col.best,cex=2)
      }
    }
    nemenyi.titles(xlab,ylab,mm$title.line,ax)

  }

  return(structure(list("means"=ranks.means,"intervals"=ranks.intervals,"fpval"=fried.pval,"fH"=fried.H,
                        "cd"=r.stat,"conf.level"=conf.level,"k"=cols.number,"n"=rows.number,
                        "cld"=ranks.cld,"diff"=ranks.diff),class="nemenyi"))

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
       mgp=getp("mgp"), family=getp("family"))
}

# Width of text drawn perpendicular to an axis, in margin lines
nemenyi.textlines <- function(txt, ax){
  w <- max(graphics::strwidth(txt, units="inches", cex=ax$cex.axis,
                              font=ax$font.axis, family=ax$family))
  w / (graphics::par("csi")*graphics::par("mex"))
}

# Enlarge margins that cannot fit axis labels (and axis title outside them)
# Other margins are kept as set by the user
nemenyi.margins <- function(mar, sides, txt, titles, ax){
  title.line <- rep(NA,4)
  for (s in seq_along(sides)){
    side <- sides[s]
    need <- ax$mgp[2] + nemenyi.textlines(txt[[s]],ax) + 0.3
    if (nzchar(titles[[s]])){
      title.line[side] <- need
      need <- need + ax$cex.lab + 0.3
    }
    mar[side] <- max(mar[side], need)
  }
  list(mar=mar, title.line=title.line)
}

# Draw an axis with the user's annotation settings
nemenyi.axis <- function(side, ax, ...){
  graphics::axis(side, cex.axis=ax$cex.axis, font.axis=ax$font.axis,
                 col.axis=ax$col.axis, mgp=ax$mgp, family=ax$family, ...)
}

# Draw axis titles, outside the labels when needed
nemenyi.titles <- function(xlab, ylab, title.line, ax){
  if (nzchar(xlab)){
    graphics::title(xlab=xlab, line=title.line[1], cex.lab=ax$cex.lab,
                    font.lab=ax$font.lab, col.lab=ax$col.lab, mgp=ax$mgp, family=ax$family)
  }
  if (nzchar(ylab)){
    graphics::title(ylab=ylab, line=title.line[2], cex.lab=ax$cex.lab,
                    font.lab=ax$font.lab, col.lab=ax$col.lab, mgp=ax$mgp, family=ax$family)
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
  if (!is.null(x$cld)){
    writeLines("")
    out <- data.frame("Mean rank"=format(round(x$means,4),nsmall=4), "Group"=x$cld, check.names=FALSE)
    rownames(out) <- names(x$means)
    print(out)
  }

}
