#' Sensitivity of the Nemenyi test to the confidence level
#'
#' Shows how the grouping of treatments by the Nemenyi test changes with the confidence level. As the test uses a single critical distance for all pairs, each pair of treatments changes from significantly different to not significantly different at exactly one confidence level, given by its pairwise p-value. The groupings therefore change only at these levels and are constant in between, so all distinct groupings can be found exactly.
#'
#' @param data an array that includes values to be compared for several treatments (in columns) for several observations (rows), of size n x k. See \code{\link{nemenyi}}.
#' @param conf.level confidence levels to show. If \code{NULL} (default), one level is shown for each distinct grouping within \code{conf.range}, labelled with the range of confidence levels that give that grouping.
#' @param conf.range range of confidence levels searched for distinct groupings, when \code{conf.level} is \code{NULL}.
#' @param ref optional reference confidence level, e.g. \code{0.95}. If given, the column of this level is outlined in the \code{"groups"} and \code{"letters"} plots, and the matching significance level is marked on the colour key of the \code{"pvalue"} plot. If \code{NULL} (default), nothing is highlighted.
#' @param plottype type of plot to produce:
#' \itemize{
#'  \item \code{"none"}: no plot.
#'  \item \code{"groups"}: \emph{group stability} plot. Columns are confidence levels and rows are treatments in mean rank order, best at the top. In each column, every vertical bar connects a group of treatments that are all not significantly different, as in the \emph{line} plot of \code{\link{nemenyi}}.
#'  \item \code{"letters"}: as \code{"groups"}, but showing the compact letter display in each column: treatments that share a letter are not significantly different.
#'  \item \code{"pvalue"}: lower triangular matrix of the pairwise Nemenyi p-values. The shading increases with the p-value, on a log scale; shaded cells indicate little evidence of significant difference. A pair is not significantly different at a confidence level above 1 - p-value.
#'  }
#' Changes from the previous column are highlighted. In \code{"groups"} these are the treatments whose set of not significantly different treatments changed, i.e. the members of the pairs that changed outcome, which mark the edges of the groups that changed. In \code{"letters"} these are all the letters that changed, which also includes treatments whose letters were relabelled because other groups changed. In both plots, columns where the Friedman test is not significant are shaded grey, and, if \code{ref} is given, the column of that confidence level is outlined.
#' @param labels optional labels for models. If NULL column names of \code{data} will be used.
#' @param line margin lines where the plot annotation is placed, see \code{\link{nemenyi}}. For \code{"pvalue"}, \code{group} is the distance of the colour key from the plot.
#' @param ... additional arguments passed to the \code{plot} function, see \code{\link{nemenyi}}.
#'
#' @details
#' The pairwise p-value of treatments i and j is \code{1 - ptukey(|R_i - R_j| / sqrt(k(k+1)/(12n)), k, Inf)}, where \code{R} are the mean ranks. The pair is not significantly different at confidence levels above \code{1 - p}, which is the same outcome as \code{\link{nemenyi}} at that level. Graphical parameters and colours are handled as in \code{\link{nemenyi}}.
#'
#' @return Return object of class \code{nemenyi.sens} and contains:
#' \itemize{
#' \item \code{means}: mean rank of each treatment.
#' \item \code{pval}: k x k matrix of pairwise Nemenyi p-values.
#' \item \code{pairs}: table of all pairs, ordered by p-value, with the absolute difference in mean ranks, the p-value and the confidence level above which the pair is not significantly different.
#' \item \code{fpval}: Friedman test p-value. The Friedman test is significant for confidence levels below \code{1 - fpval}.
#' \item \code{levels}: confidence levels shown.
#' \item \code{range}: range of confidence levels that each shown level represents. For user specified levels this is the level itself.
#' \item \code{cd}: Nemenyi critical distance at each shown level.
#' \item \code{group}: k x levels matrix of compact letter displays.
#' \item \code{changed}: list of two k x levels logical matrices, \code{TRUE} where a treatment changed from the previous level: \code{groups}, if its set of not significantly different treatments changed (highlighted in the \code{"groups"} plot), and \code{letters}, if its letters changed (highlighted in the \code{"letters"} plot).
#' \item \code{k}: number of treatments (columns).
#' \item \code{n}: number of observations (rows).
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
#'
#' @seealso \code{\link{nemenyi}}
#'
#' @keywords htest
#'
#' @examples
#' x <- matrix( rnorm(50*4,mean=0,sd=1), 50, 4)
#' x[,2] <- x[,2]+1
#' x[,3] <- x[,3]+0.7
#' x[,4] <- x[,4]+0.5
#' colnames(x) <- c("Method A","Method B","Method C - long name","Method D")
#' nemenyi.sens(x)
#' nemenyi.sens(x,plottype="letters")
#' nemenyi.sens(x,plottype="pvalue",ref=0.95)
#' nemenyi.sens(x,conf.level=c(0.8,0.9,0.95,0.99),ref=0.95)
#'
#' @export nemenyi.sens

nemenyi.sens <- function(data, conf.level=NULL, conf.range=c(0.80,0.999), ref=NULL,
                         plottype=c("groups","letters","pvalue","none"), labels=NULL,
                         line=c("main"=NA,"xlab"=NA,"ylab"=NA,"xaxis"=NA,"yaxis"=NA,"group"=NA), ...){

  plottype <- match.arg(plottype,c("groups","letters","pvalue","none"))
  line <- nemenyi.line(line)

  # Check reference level
  if (!is.null(ref) && ((length(ref) != 1) || (ref <= 0) || (ref >= 1))){
    stop("ref must be NULL or a single confidence level between 0 and 1.")
  }

  # Mean ranks and Friedman test, these do not depend on the confidence level
  res <- nemenyi(data, plottype="none", labels=labels)
  means <- res$means
  k <- res$k
  n <- res$n
  lbl <- names(means)
  se <- sqrt((k*(k+1))/(12*n))

  # Pairwise p-values: a pair is not significantly different for conf.level > 1 - p
  pval <- 1 - stats::ptukey(res$diff/se, k, Inf)
  dimnames(pval) <- dimnames(res$diff)
  ut <- which(upper.tri(pval), arr.ind=TRUE)
  pairs <- data.frame("Treatment 1"=lbl[ut[,1]], "Treatment 2"=lbl[ut[,2]],
                      "Difference"=res$diff[ut], "p-value"=pval[ut], "Same above"=1-pval[ut],
                      check.names=FALSE)
  pairs <- pairs[order(pairs[["p-value"]]),]
  rownames(pairs) <- NULL

  # Confidence levels to show
  if (is.null(conf.level)){
    if ((length(conf.range) != 2) || any(conf.range <= 0) || any(conf.range >= 1) || (conf.range[1] >= conf.range[2])){
      stop("conf.range must be two increasing values between 0 and 1.")
    }
    # Groupings change only at the critical levels of the pairs, and the Friedman
    # outcome at its own critical level
    crit <- sort(unique(c(pairs[["Same above"]], 1 - res$fpval)))
    crit <- crit[(crit > conf.range[1]) & (crit < conf.range[2])]
    edges <- c(conf.range[1], crit, conf.range[2])
    lo <- edges[-length(edges)]
    hi <- edges[-1]
    levels <- (lo + hi)/2
    e <- nemenyi.pct(edges)
    lab <- paste0(e[-length(e)],"-",e[-1])
    ref.col <- if (is.null(ref)) integer(0) else which((lo < ref) & (ref <= hi))
  } else {
    if (any(conf.level <= 0) || any(conf.level >= 1)){
      stop("conf.level must be between 0 and 1.")
    }
    levels <- lo <- hi <- sort(unique(conf.level))
    lab <- nemenyi.pct(levels)
    ref.col <- if (is.null(ref)) integer(0) else which(abs(levels - ref) < 1e-10)
  }
  L <- length(levels)

  # Groups at each level
  cd <- stats::qtukey(levels, k, Inf)*se
  grp <- matrix(sapply(cd, function(x){nemenyi.cld(means, x)}), nrow=k, dimnames=list(lbl, lab))
  # Changes from the previous level
  # groups: treatments whose set of not significantly different treatments changed
  # letters: treatments whose letters changed, including relabelling
  changed <- list("groups"=matrix(FALSE, k, L, dimnames=list(lbl, lab)),
                  "letters"=matrix(FALSE, k, L, dimnames=list(lbl, lab)))
  if (L > 1){
    for (l in 2:L){
      changed$groups[,l] <- rowSums((res$diff < cd[l]) != (res$diff < cd[l-1])) > 0
      changed$letters[,l] <- grp[,l] != grp[,l-1]
    }
  }
  # Friedman test not significant at these levels
  fried.ns <- levels > 1 - res$fpval

  # Plots
  if (plottype != "none"){
    args <- list(...)

    if (is.null(args[["main"]])){
      main <- paste0("Friedman: ", format(round(res$fpval,3),nsmall=3), "\n",
                     if (plottype == "pvalue") "Nemenyi pairwise p-values" else "Nemenyi groups across confidence levels")
    } else {
      main <- args[["main"]]
    }
    args$main <- NULL
    if (is.null(args[["xaxs"]])){
      args$xaxs <- "i"
    }
    if (is.null(args[["yaxs"]])){
      args$yaxs <- "i"
    }
    xlab <- if (is.null(args[["xlab"]])) "" else args[["xlab"]]
    ylab <- if (is.null(args[["ylab"]])) "" else args[["ylab"]]
    if ((plottype != "pvalue") && is.null(args[["xlab"]])){
      xlab <- "Confidence level (%)"
    }
    args$xlab <- args$ylab <- ""
    args$x <- args$y <- NA
    args$axes <- FALSE

    ax <- nemenyi.axpar(args)
    mar.user <- graphics::par("mar")
    on.exit(graphics::par(mar=mar.user), add=TRUE)

    fg <- graphics::par("fg")
    col.best <- RColorBrewer::brewer.pal(3,"Set1")[1]
    col.ns <- RColorBrewer::brewer.pal(3,"Set1")[2]
    lbl.rank <- paste0(lbl," - ",sprintf("%1.2f",round(means,2)))
    if (is.null(args[["ylim"]])){
      args$ylim <- c(0.5,k+0.5)
    }
  }

  # Group stability and letters plots
  if ((plottype == "groups") | (plottype == "letters")){

    lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,xtxt=lab,ytxt=lbl.rank)
    graphics::par(mar=lay$mar)
    if (is.null(args[["xlim"]])){
      args$xlim <- c(0.5,L+0.5)
    }
    do.call(graphics::plot,args)

    # Columns where Friedman is not significant and reference level. Fills are
    # opaque, so they are drawn first and the reference fill is mixed with the grey
    fill.fns <- colour.blend("gray50",0.2)
    for (l in which(fried.ns)){
      graphics::rect(l-0.5,0.5,l+0.5,k+0.5,col=fill.fns,border=NA)
    }
    if (length(ref.col) > 0){
      fill.ref <- colour.blend(col.best,0.08,bg=if (fried.ns[ref.col]) fill.fns else graphics::par("bg"))
      graphics::rect(ref.col-0.5,0.5,ref.col+0.5,k+0.5,col=fill.ref,border=NA)
    }
    # Guides
    graphics::abline(h=1:k,col="gray85",lty=3)
    if (L > 1){
      graphics::abline(v=(1:(L-1))+0.5,col="gray70")
    }
    if (length(ref.col) > 0){
      graphics::rect(ref.col-0.5,0.5,ref.col+0.5,k+0.5,col=NA,border=col.best,lty=2,lwd=1.5)
    }

    y <- k:1
    for (l in 1:L){
      bcol <- if (fried.ns[l]) "gray50" else col.ns
      if (plottype == "groups"){
        # Bars for groups with more than one member
        g <- nemenyi.groups(means, cd[l])
        g <- g[g[,1] != g[,2],,drop=FALSE]
        m <- nrow(g)
        if (m > 0){
          xs <- l - 0.2 + (1:m)*(0.55/(m+1))
          for (i in 1:m){
            graphics::lines(rep(xs[i],2), y[g[i,]], col=bcol, lwd=3)
            graphics::points(rep(xs[i],2), y[g[i,]], col=bcol, pch=20)
          }
        }
        chg <- changed$groups[,l]
        graphics::points(rep(l-0.32,k), y, pch=20, col=ifelse(chg,col.best,fg),
                         cex=ifelse(chg,1.6,1))
      } else {
        chg <- changed$letters[,l]
        graphics::text(rep(l,k), y, grp[,l], col=ifelse(chg,col.best,fg),
                       font=ifelse(chg,2,1), cex=ax$cex.axis, family=ax$family)
      }
    }

    nemenyi.axis(1,lay,ax,at=1:L,labels=lab,las=2)
    nemenyi.axis(2,lay,ax,at=y,labels=lbl.rank,las=2)
    graphics::box(which="plot")
    nemenyi.titles(xlab,ylab,main,lay,ax)

  }

  # Pairwise p-value matrix
  if (plottype == "pvalue"){

    # Shading on log10(p), from 0.001 (none) to 1 (full)
    pmin <- 0.001
    shade <- function(p){
      t <- (log10(pmax(p,pmin)) - log10(pmin))/(-log10(pmin))
      colour.blend(col.ns, 0.05+0.85*t)
    }
    key.at <- c(0.001,0.01,0.05,0.1,1)
    key.lab <- c("<0.001","0.01","0.05","0.1","1")

    # Margins: labels on both axes and colour key with its labels on the right
    if (is.na(line["group"])){
      line["group"] <- 0.8
    }
    lay <- nemenyi.layout(mar.user,line,ax,xlab,ylab,main,xtxt=lbl,ytxt=lbl.rank)
    key.pos <- lay$line[["group"]]
    lay$mar[4] <- max(lay$mar[4], key.pos + 1.3 + nemenyi.textlines(key.lab,ax) + 0.3)
    graphics::par(mar=lay$mar)
    if (is.null(args[["xlim"]])){
      args$xlim <- c(0.5,k+0.5)
    }
    do.call(graphics::plot,args)

    cex.cell <- min(ax$cex.axis, 0.85/graphics::strwidth("<0.001"), 0.6/graphics::strheight("0.00"))
    # Cells are not clipped and axes have no line, so all borders have the same width
    for (i in 1:k){
      yy <- k - i + 1
      for (j in 1:i){
        if (i == j){
          graphics::rect(j-0.5,yy-0.5,j+0.5,yy+0.5,col=fg,border=fg,xpd=NA)
        } else {
          graphics::rect(j-0.5,yy-0.5,j+0.5,yy+0.5,col=shade(pval[i,j]),border=fg,xpd=NA)
          if (cex.cell >= 0.4){
            ptxt <- if (pval[i,j] < pmin) "<0.001" else sprintf("%1.3f",pval[i,j])
            graphics::text(j,yy,ptxt,cex=cex.cell,col=fg,family=ax$family)
          }
        }
      }
    }

    # Colour key in the right margin
    usr <- graphics::par("usr")
    ux <- function(lines){usr[2] + lines*graphics::par("csi")*graphics::par("mex")/(graphics::par("pin")[1]/diff(usr[1:2]))}
    ky <- function(p){0.5 + k*(log10(pmax(p,pmin)) - log10(pmin))/(-log10(pmin))}
    pseq <- 10^seq(log10(pmin),0,length.out=101)
    for (s in 1:100){
      graphics::rect(ux(key.pos),ky(pseq[s]),ux(key.pos+0.8),ky(pseq[s+1]),col=shade(pseq[s+1]),border=NA,xpd=NA)
    }
    graphics::rect(ux(key.pos),0.5,ux(key.pos+0.8),k+0.5,border=fg,xpd=NA)
    graphics::segments(ux(key.pos+0.8),ky(key.at),ux(key.pos+1),ky(key.at),col=fg,xpd=NA)
    graphics::text(ux(key.pos+1.2),ky(key.at),key.lab,adj=c(0,0.5),xpd=NA,cex=ax$cex.axis,
                   font=ax$font.axis,col=ax$col.axis,family=ax$family)
    # Reference significance level
    if (!is.null(ref)) graphics::segments(ux(key.pos-0.1),ky(1-ref),ux(key.pos+0.9),ky(1-ref),col=col.best,lwd=2,xpd=NA)

    nemenyi.axis(1,lay,ax,at=1:k,labels=lbl,las=2,lwd=0,lwd.ticks=graphics::par("lwd"))
    nemenyi.axis(2,lay,ax,at=k:1,labels=lbl.rank,las=2,lwd=0,lwd.ticks=graphics::par("lwd"))
    nemenyi.titles(xlab,ylab,main,lay,ax)

  }

  return(structure(list("means"=means,"pval"=pval,"pairs"=pairs,"fpval"=res$fpval,
                        "levels"=levels,"range"=rbind("from"=lo,"to"=hi),"cd"=cd,
                        "group"=grp,"changed"=changed,"k"=k,"n"=n),class="nemenyi.sens"))

}

# Confidence levels as percentages, with enough digits so that levels below 100%
# are not shown as 100% and, if distinct=TRUE, different levels have different labels
nemenyi.pct <- function(x, distinct=TRUE){
  ok <- function(s, v){all(v >= 1 | as.numeric(s) < 100)}
  if (distinct){
    u <- unique(x)
    for (d in 1:6){
      s <- formatC(100*u, format="f", digits=d)
      if (!anyDuplicated(s) && ok(s, u)) break
    }
    return(formatC(100*x, format="f", digits=d))
  }
  # Otherwise about two significant digits in the distance from 100%
  sapply(x, function(v){
    d <- if ((v <= 0) || (v >= 1)) 1 else min(6, max(1, ceiling(-log10(1-v))))
    formatC(100*v, format="f", digits=d)
  })
}

#' @export
#' @method summary nemenyi.sens
summary.nemenyi.sens <- function(object,...){
  print(object)
}

#' @export
#' @method print nemenyi.sens
print.nemenyi.sens <- function(x,...){

  writeLines("Friedman and Nemenyi Tests: sensitivity to the confidence level")
  writeLines(paste0("Number of observations is ", x$n, " and number of methods is ", x$k))
  writeLines(paste0("Friedman test p-value: ", format(round(x$fpval,4),nsmall=4),
                    " (significant for confidence levels below ", nemenyi.pct(1-x$fpval,FALSE), "%)"))
  writeLines("")
  writeLines("Pairwise Nemenyi p-values (pairs are not significantly different above the confidence level shown):")
  out <- x$pairs
  out[["Difference"]] <- format(round(out[["Difference"]],4),nsmall=4)
  out[["p-value"]] <- format(round(out[["p-value"]],4),nsmall=4)
  out[["Same above"]] <- paste0(nemenyi.pct(x$pairs[["Same above"]],FALSE),"%")
  print(out, row.names=FALSE)
  writeLines("")
  writeLines("Groups across confidence levels (%):")
  out <- data.frame("Mean rank"=format(round(x$means,4),nsmall=4), x$group, check.names=FALSE)
  print(out)

}
