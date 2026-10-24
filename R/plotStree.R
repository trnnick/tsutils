#' Plot hierarchy as a tree
#'
#' Plots a hierarchy from its summing matrix \code{S} as the familiar tree of nodes and edges. Each node is a row of \code{S}, levels are arranged from the top to the bottom of the hierarchy, and each node is connected to the nodes directly below it. For large hierarchies \code{\link{plotSblock}} is usually cleaner.
#'
#' @param S summing matrix (\eqn{n \times m}). Non-zero entries indicate that the bottom level series in that column is aggregated by the node in that row. This can describe cross-sectional, grouped or temporal hierarchies.
#' @param labels \code{TRUE} to use the row names of \code{S} as node labels (no labels if \code{S} has no row names), \code{FALSE} for no labels, or one label per row of \code{S}. Labels can be a character vector, an expression vector, or a list mixing strings, quoted calls and expressions.
#' @param parse if \code{TRUE} and \code{labels=TRUE}, row names are parsed as \link[grDevices]{plotmath}, e.g. \code{"H[1]"} is plotted as H with subscript 1. Names that are not valid plotmath are plotted as text, with a warning.
#' @param skip how to plot edges that skip levels:
#' \itemize{
#' \item \code{"empty"}: the edge goes directly from the parent to the child.
#' \item \code{"expand"}: the child gets a vertical stem that rises through every level above it that covers none of its bottom level series, as in \code{\link{plotSblock}}. The edge joins the top of the stem.
#' }
#' With \code{edges="elbow"} this sets which end has the long vertical segment: the parent (\code{"empty"}) or the child (\code{"expand"}).
#' @param levels optional vector with the level of each row of \code{S}. Levels are ordered top to bottom in order of appearance. If \code{NULL}, levels are found automatically, as in \code{\link{plotSblock}}.
#' @param layout horizontal position of the nodes:
#' \itemize{
#' \item \code{"span"}: each node is centred over the bottom level series it aggregates, so that the plot lines up with the columns of \code{S} and with \code{\link{plotSblock}}. Nodes closer than one column apart are pushed apart.
#' \item \code{"even"}: nodes of each level are spread evenly, keeping their order.
#' }
#' @param edges edge style: \code{"straight"} or \code{"elbow"}. Elbows give an organisation chart look that suits hierarchies where each node has a single parent. When nodes have several parents (grouped and temporal hierarchies) straight edges are clearer.
#' @param shape node shape: \code{"box"} or \code{"circle"}.
#' @param fill fill colour of the nodes. Either a single colour or one per row of \code{S}.
#' @param border border colour of the nodes.
#' @param edge.col edge colour. Either a single colour, or one per level, in which case the edge takes the colour of the level of its parent. If \code{NULL}, shades of grey from dark (top) to light (bottom).
#' @param edge.lwd edge line width.
#' @param hgap horizontal space between nodes, as a fraction of the distance between neighbouring nodes in a level.
#' @param vgap vertical space between levels, as a fraction of the height of a level.
#' @param node.size size of nodes without labels, in inches. Nodes with labels are sized to fit them. In both cases nodes are not larger than \code{hgap} and \code{vgap} allow.
#' @param cex maximum label size.
#' @param min.cex labels that would need to be smaller than this to fit in their node are not plotted.
#' @param level.labels \code{FALSE} for none, \code{TRUE} for \code{"L1", "L2", ...}, or one label per level, starting from the top, as a character or expression vector. These are plotted in the left margin.
#' @param level.font font of the level labels, as in \code{par("font")}: 1 plain, 2 bold, 3 italic, 4 bold italic.
#' @param main optional plot title.
#'
#' @details
#' A node is connected to another node below it if it aggregates all the bottom level series of that node, and there is no other node in between that does so as well. Therefore, each node is connected only to its nearest parents. In temporal hierarchies this gives the usual lattice, e.g., a monthly node is connected to both its 2-monthly and its quarterly node, but not to the annual node. In grouped hierarchies each node has one parent per grouping. Edges are shaded by the level of their parent, which helps to follow the different groupings.
#'
#' @return Invisibly, a list with:
#' \itemize{
#' \item \code{nodes}: data frame with one row per node and columns \code{node} (row name of \code{S}), \code{row} (row index in \code{S}), \code{level}, \code{x}, \code{y} (centre of the node), \code{width} and \code{height}. Positions and sizes are in plot coordinates, with one unit per column of \code{S} and per level.
#' \item \code{edges}: data frame with one row per edge and columns \code{from}, \code{to} (row names of \code{S}), \code{from.row} and \code{to.row} (row indices in \code{S}).
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{plotSblock}}, \code{\link{plotS}}, \code{\link{plotThief}}.
#'
#' @keywords ts
#'
#' @examples
#' # Helper to construct S from a list of aggregates
#' makeS <- function(aggs,bottom){
#'   A <- t(sapply(aggs,function(a){as.numeric(bottom %in% a)}))
#'   S <- rbind(A,diag(length(bottom)))
#'   dimnames(S) <- list(c(names(aggs),bottom),bottom)
#'   return(S)
#' }
#'
#' # Simple hierarchy
#' S1 <- makeS(list(A=c("Q1","Q2","Q3","Q4"),H1=c("Q1","Q2"),H2=c("Q3","Q4")),
#'             c("Q1","Q2","Q3","Q4"))
#' plotStree(S1)
#' plotStree(S1,edges="elbow")
#'
#' # Row names parsed as plotmath, circles and level labels
#' rownames(S1) <- c("A","H[1]","H[2]","Q[1]","Q[2]","Q[3]","Q[4]")
#' plotStree(S1,parse=TRUE,shape="circle",
#'           level.labels=expression(bold(S)[4],bold(S)[2],bold(S)[1]))
#'
#' # Q5 skips a level and connects directly to A
#' S2 <- makeS(list(A=c("Q1","Q2","Q3","Q4","Q5"),H1=c("Q1","Q2"),H2=c("Q3","Q4")),
#'             c("Q1","Q2","Q3","Q4","Q5"))
#' plotStree(S2,skip="empty")
#' tr <- plotStree(S2,skip="expand")
#' tr$edges
#'
#' # Temporal hierarchy for monthly data
#' plotStree(Sthief(12),labels=FALSE,shape="circle",level.labels=paste("k =",c(12,6,4,3,2,1)))
#'
#' # Grouped hierarchy: region x product, with colours for each grouping
#' bottom <- c("R1P1","R1P2","R2P1","R2P2")
#' S3 <- makeS(list(Total=bottom,R1=c("R1P1","R1P2"),R2=c("R2P1","R2P2"),
#'                  P1=c("R1P1","R2P1"),P2=c("R1P2","R2P2")),bottom)
#' cmp.fill <- RColorBrewer::brewer.pal(5,"Pastel1")
#' cmp.edge <- RColorBrewer::brewer.pal(5,"Set1")
#' plotStree(S3,fill=c("grey88",cmp.fill[2],cmp.fill[2],cmp.fill[5],cmp.fill[5],rep("grey88",4)),
#'           edge.col=c("grey20",cmp.edge[2],cmp.edge[5]))
#'
#' # Larger hierarchy without labels
#' b <- paste0("b",1:60)
#' S4 <- makeS(c(list(Total=b),
#'               setNames(lapply(0:4,function(g){b[g*12+1:12]}),paste0("G",1:5)),
#'               setNames(lapply(0:19,function(g){b[g*3+1:3]}),paste0("SG",1:20))),b)
#' plotStree(S4,labels=FALSE,shape="circle",fill=RColorBrewer::brewer.pal(3,"Pastel1")[2])
#' plotStree(S4,labels=FALSE,edges="elbow")
#'
#' @export plotStree

plotStree <- function(S,labels=TRUE,parse=c(FALSE,TRUE),skip=c("empty","expand"),
                      levels=NULL,layout=c("span","even"),edges=c("straight","elbow"),
                      shape=c("box","circle"),fill="grey88",border="grey30",
                      edge.col=NULL,edge.lwd=1.5,hgap=0.2,vgap=0.5,node.size=0.2,
                      cex=1,min.cex=0.35,level.labels=FALSE,level.font=1,main=NULL){

  # Defaults
  parse <- parse[1]
  skip <- match.arg(skip,c("empty","expand"))
  layout <- match.arg(layout,c("span","even"))
  edges <- match.arg(edges,c("straight","elbow"))
  shape <- match.arg(shape,c("box","circle"))

  # Get structure, labels and levels
  prep <- plotS.prep(S,labels,parse,levels)
  M <- prep$M
  n <- prep$n
  m <- prep$m
  L <- prep$L
  lev <- prep$lev
  lab.list <- prep$lab.list
  fill <- rep_len(fill,n)
  use <- setdiff(1:n,prep$empty)      # Rows to plot

  # Find edges
  # A node is a parent of another if it is in a higher level and aggregates all its bottom series.
  # Keep only the nearest parents, i.e., drop edges that have a node in between.
  Mn <- M*1
  R <- (Mn %*% t(Mn) == matrix(rowSums(Mn),n,n,byrow=TRUE)) & outer(lev,lev,"<")
  R[prep$empty,] <- FALSE
  R[,prep$empty] <- FALSE
  R <- R & ((R*1) %*% (R*1) == 0)
  ed <- which(R,arr.ind=TRUE)
  ed <- ed[order(lev[ed[,1]],ed[,1],lev[ed[,2]],ed[,2]),,drop=FALSE]

  # Horizontal location of nodes
  cen <- rep(NA,n)
  for (i in use){
    cen[i] <- mean(which(M[i,])) - 0.5
  }
  x <- rep(NA,n)
  spacing <- rep(m,L)                 # Smallest distance between nodes in each level
  for (l in 1:L){
    idx <- use[lev[use]==l]
    if (length(idx) == 0){
      next
    }
    idx <- idx[order(cen[idx])]
    K <- length(idx)
    if (layout == "even"){
      xx <- ((1:K)-0.5)*m/K
    } else {
      # Centre over the bottom series, but push apart nodes that are too close (grouped hierarchies)
      d <- min(1,m/K)
      xx <- cen[idx]
      if (K > 1){
        for (j in 2:K){
          xx[j] <- max(xx[j],xx[j-1]+d)
        }
        xx[K] <- min(xx[K],m-d/2)
        for (j in (K-1):1){
          xx[j] <- min(xx[j],xx[j+1]-d)
        }
      }
    }
    x[idx] <- xx
    if (K > 1){
      spacing[l] <- min(diff(xx))
    }
  }
  y <- L - lev + 0.5

  # Start plotting
  op <- plotS.setup(L,m,level.labels,level.font,main)
  on.exit(par(op))

  # Node sizes
  # Work in inches, so that circles are circles, and then convert to plot units
  xpi <- diff(par("usr")[1:2])/par("pin")[1]
  ypi <- diff(par("usr")[3:4])/par("pin")[2]
  lab.w <- lab.h <- rep(0,n)
  for (i in use){
    if (!is.null(lab.list[[i]])){
      lab.w[i] <- strwidth(lab.list[[i]],units="inches",cex=cex)
      lab.h[i] <- strheight(lab.list[[i]],units="inches",cex=cex)
    }
  }
  pad <- 0.5*strheight("M",units="inches",cex=cex)
  w <- h <- rep(NA,n)
  for (l in 1:L){
    idx <- use[lev[use]==l]
    if (length(idx) == 0){
      next
    }
    # Space available and space needed for the labels. All nodes in a level get the same size.
    w.av <- spacing[l]*(1-hgap)/xpi
    h.av <- (1-vgap)/ypi
    has.lab <- !sapply(lab.list[idx],is.null)
    if (!any(has.lab)){
      w.need <- h.need <- node.size
    } else if (shape == "circle"){
      w.need <- h.need <- max(sqrt(lab.w[idx]^2+lab.h[idx]^2)) + pad
    } else {
      h.need <- max(lab.h[idx]) + 2*pad
      w.need <- max(max(lab.w[idx])+2*pad,h.need)
    }
    if (shape == "circle" || !any(has.lab)){
      w[idx] <- h[idx] <- min(w.av,h.av,w.need)
    } else {
      w[idx] <- min(w.av,w.need)
      h[idx] <- min(h.av,h.need)
    }
  }
  hw <- w*xpi/2                       # Half width and height in plot units
  hh <- h*ypi/2

  # Plot edges
  if (is.null(edge.col)){
    edge.col <- gray(seq(0.15,0.6,length.out=max(L-1,1)))
  }
  edge.col <- rep_len(edge.col,L)
  # Edges join boxes at their border and circles at their centre (circles are plotted on top)
  if (shape == "box"){
    y.bot <- y - hh
    y.top <- y + hh
  } else {
    y.bot <- y.top <- y
  }
  sg <- NULL                          # Segments: x0, y0, x1, y1, level of parent
  for (e in seq_len(nrow(ed))){
    i <- ed[e,1]
    j <- ed[e,2]
    if (edges == "straight"){
      # With skip="expand" the child gets a stem through the levels above that do not cover it
      top <- lev[j]
      if (skip == "expand"){
        cols <- which(M[j,])
        k <- lev[j] - 1
        while (k > lev[i] && !any(prep$coverage[k,cols])){
          top <- k
          k <- k - 1
        }
      }
      if (top < lev[j]){
        ys <- L - top + 0.5 + (1-vgap)/2
        sg <- rbind(sg,c(x[j],y.top[j],x[j],ys,lev[i]),c(x[i],y.bot[i],x[j],ys,lev[i]))
      } else {
        sg <- rbind(sg,c(x[i],y.bot[i],x[j],y.top[j],lev[i]))
      }
    } else {
      # Elbows: horizontal bar above the child (skip="empty") or below the parent (skip="expand")
      # Bars of different levels in the same gap are spread out so that they do not merge
      if (skip == "empty"){
        yb <- L - lev[j] + 1
        key <- lev[i]
        keys <- sort(unique(lev[ed[lev[ed[,2]]==lev[j],1]]))
      } else {
        yb <- L - lev[i]
        key <- -lev[j]
        keys <- sort(unique(-lev[ed[lev[ed[,1]]==lev[i],2]]))
      }
      if (length(keys) > 1){
        yb <- yb + seq(0.35,-0.35,length.out=length(keys))[keys==key]*vgap
      }
      sg <- rbind(sg,c(x[i],y.bot[i],x[i],yb,lev[i]),c(x[i],yb,x[j],yb,lev[i]),
                  c(x[j],yb,x[j],y.top[j],lev[i]))
    }
  }
  if (!is.null(sg)){
    segments(sg[,1],sg[,2],sg[,3],sg[,4],col=edge.col[sg[,5]],lwd=edge.lwd,
             lend=ifelse(edges=="elbow",2,1))
  }

  # Plot nodes
  if (shape == "box"){
    rect(x[use]-hw[use],y[use]-hh[use],x[use]+hw[use],y[use]+hh[use],col=fill[use],border=border)
  } else {
    th <- seq(0,2*pi,length.out=61)
    for (i in use){
      polygon(x[i]+hw[i]*cos(th),y[i]+hh[i]*sin(th),col=fill[i],border=border)
    }
  }

  # Add labels, shrinking them to fit
  # For circles use the square inside the circle
  shrink <- ifelse(shape=="circle",sqrt(0.5),1)
  for (i in use){
    lab <- lab.list[[i]]
    if (is.null(lab)){
      next
    }
    fit <- min(0.9*shrink*2*hw[i]/strwidth(lab,cex=1),0.8*shrink*2*hh[i]/strheight(lab,cex=1))
    cex.i <- min(cex,fit)
    if (cex.i >= min.cex){
      text(x[i],y[i],lab,cex=cex.i)
    }
  }

  # Prepare output
  nodes <- data.frame(node=prep$node.names[use],row=use,level=lev[use],x=x[use],y=y[use],
                      width=2*hw[use],height=2*hh[use],stringsAsFactors=FALSE)
  links <- data.frame(from=prep$node.names[ed[,1]],to=prep$node.names[ed[,2]],
                      from.row=ed[,1],to.row=ed[,2],stringsAsFactors=FALSE)
  rownames(links) <- NULL

  return(invisible(list(nodes=nodes,edges=links)))

}
