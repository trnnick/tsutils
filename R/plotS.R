#' Plot hierarchy
#'
#' Plots a hierarchy from its summing matrix \code{S}, either as a tree or as blocks. This is a convenience function that calls \code{\link{plotStree}} or \code{\link{plotSblock}}. These can also be used directly, and their help files describe all the options available.
#'
#' @param S summing matrix (\eqn{n \times m}). Non-zero entries indicate that the bottom level series in that column is aggregated by the node in that row. This can describe cross-sectional, grouped or temporal hierarchies.
#' @param view type of plot:
#' \itemize{
#' \item \code{"tree"}: nodes connected with edges, see \code{\link{plotStree}}.
#' \item \code{"block"}: each node is a block spanning the bottom level series it aggregates, see \code{\link{plotSblock}}. This is cleaner for large hierarchies.
#' }
#' @param labels \code{TRUE} to use the row names of \code{S} as node labels (no labels if \code{S} has no row names), \code{FALSE} for no labels, or one label per row of \code{S}. Labels can be a character vector, an expression vector, or a list mixing strings, quoted calls and expressions.
#' @param parse if \code{TRUE} and \code{labels=TRUE}, row names are parsed as \link[grDevices]{plotmath}, e.g. \code{"H[1]"} is plotted as H with subscript 1.
#' @param skip how to plot nodes that skip a level: \code{"empty"} or \code{"expand"}. See \code{\link{plotStree}} and \code{\link{plotSblock}} for details.
#' @param levels optional vector with the level of each row of \code{S}. If \code{NULL}, levels are found automatically.
#' @param fill fill colour of the nodes. Either a single colour or one per row of \code{S}.
#' @param border border colour of the nodes.
#' @param level.labels \code{FALSE} for none, \code{TRUE} for \code{"L1", "L2", ...}, or one label per level, starting from the top. These are plotted in the left margin.
#' @param level.font font of the level labels, as in \code{par("font")}: 1 plain, 2 bold, 3 italic, 4 bold italic.
#' @param main optional plot title.
#' @param ... additional arguments passed to \code{\link{plotStree}} or \code{\link{plotSblock}}, depending on \code{view}.
#'
#' @return Invisibly, the output of \code{\link{plotStree}} or \code{\link{plotSblock}}.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{plotStree}}, \code{\link{plotSblock}}, \code{\link{plotThief}}.
#'
#' @keywords ts
#'
#' @examples
#' # Small hierarchy: A = H1 + H2, H1 = Q1 + Q2, H2 = Q3 + Q4
#' S <- rbind(c(1,1,1,1),c(1,1,0,0),c(0,0,1,1),diag(4))
#' rownames(S) <- c("A","H1","H2","Q1","Q2","Q3","Q4")
#' plotS(S)
#' plotS(S,view="block")
#'
#' # Pass additional arguments to plotStree
#' plotS(S,edges="elbow",shape="circle",level.labels=TRUE)
#'
#' @export plotS

plotS <- function(S,view=c("tree","block"),labels=TRUE,parse=c(FALSE,TRUE),
                  skip=c("empty","expand"),levels=NULL,fill="grey88",border="grey30",
                  level.labels=FALSE,level.font=1,main=NULL,...){

  # Defaults
  view <- match.arg(view,c("tree","block"))
  parse <- parse[1]
  skip <- match.arg(skip,c("empty","expand"))

  # Pick plotting function
  if (view == "tree"){
    fun <- plotStree
    fun.name <- "plotStree"
  } else {
    fun <- plotSblock
    fun.name <- "plotSblock"
  }

  # Check that additional arguments are meant for the selected view
  args <- list(...)
  if (length(args) > 0){
    if (is.null(names(args)) || any(names(args) == "")){
      stop("Additional arguments must be named.")
    }
    bad <- setdiff(names(args),names(formals(fun)))
    if (length(bad) > 0){
      stop(paste0("Not an argument of ",fun.name," (view=\"",view,"\"): ",
                  paste(bad,collapse=", "),"."))
    }
  }

  # Plot
  # quote=TRUE so that expressions in labels are passed as they are
  args <- c(list(S=S,labels=labels,parse=parse,skip=skip,levels=levels,fill=fill,
                 border=border,level.labels=level.labels,level.font=level.font,
                 main=main),args)
  out <- do.call(fun,args,quote=TRUE)

  return(invisible(out))

}

# Internal function shared by plotStree and plotSblock.
# Gets the summing matrix structure, the labels and the levels.
plotS.prep <- function(S,labels,parse,levels){

  M <- as.matrix(S) != 0
  n <- nrow(M)
  m <- ncol(M)
  node.names <- rownames(S)
  if (is.null(node.names)){
    node.names <- paste0("row",1:n)
  }

  # Labels
  # Each row gets either a string, an expression or NULL (no label)
  if (is.null(labels)){
    labels <- FALSE
  }
  lab.list <- vector("list",n)
  if (isTRUE(labels)){
    if (!is.null(rownames(S))){
      for (i in 1:n){
        if (parse == TRUE){
          lab.list[[i]] <- tryCatch(parse(text=node.names[i],keep.source=FALSE),
                                    error=function(e){
                                      warning(paste0("Row name '",node.names[i],"' is not valid plotmath. Plotted as text."))
                                      node.names[i]
                                    })
        } else {
          lab.list[[i]] <- node.names[i]
        }
      }
    }
  } else if (!isFALSE(labels)){
    if (length(labels) != n){
      stop(paste0("Argument labels must be TRUE, FALSE, or have one entry per row of S (",n,")."))
    }
    for (i in 1:n){
      if (is.character(labels[[i]])){
        lab.list[[i]] <- labels[[i]]
      } else {
        lab.list[[i]] <- as.expression(labels[[i]])
      }
    }
  }

  # Rows without any bottom level series are not plotted
  empty <- which(rowSums(M) == 0)
  if (length(empty) > 0){
    warning(paste0("Rows with no bottom level series are not plotted: ",
                   paste(node.names[empty],collapse=", "),"."))
  }

  # Assign rows to levels
  # Scan S from the top and start a new level whenever a row overlaps with the current level
  if (is.null(levels)){
    lev <- rep(0,n)
    current <- 1
    covered <- rep(FALSE,m)       # Bottom series used by the current level
    for (i in 1:n){
      if (any(covered & M[i,])){
        current <- current + 1
        covered[] <- FALSE
      }
      lev[i] <- current
      covered <- covered | M[i,]
    }
  } else {
    if (length(levels) != n){
      stop("Argument levels must have one entry per row of S.")
    }
    lev <- as.integer(factor(levels,levels=unique(levels)))
  }
  L <- max(lev)

  # Which bottom level series are covered by each level
  coverage <- matrix(FALSE,L,m)
  for (l in 1:L){
    coverage[l,] <- colSums(M[lev==l,,drop=FALSE]) > 0
  }

  return(list(M=M,n=n,m=m,node.names=node.names,lab.list=lab.list,empty=empty,
              lev=lev,L=L,coverage=coverage))

}

# Internal function shared by plotStree and plotSblock.
# Opens the plot and adds the level labels. Returns the old par settings.
plotS.setup <- function(L,m,level.labels,level.font,main){

  # Level labels
  if (isFALSE(level.labels) || is.null(level.labels)){
    ln <- NULL
  } else if (isTRUE(level.labels)){
    ln <- paste0("L",1:L)
  } else {
    ln <- level.labels
    if (length(ln) != L){
      stop(paste0("Argument level.labels must have one entry per level (",L," levels found)."))
    }
  }

  # Left margin from the width of the level labels (in lines)
  left <- 0.5
  if (!is.null(ln)){
    if (dev.cur() == 1){
      dev.new()
    }
    left <- 0.8 + max(strwidth(ln,units="inches",font=level.font))/par("csi")
  }
  op <- par(mar=c(0.5,left,ifelse(is.null(main),0.5,2.5),0.5))

  # Start plot
  # Level l is plotted in the band [L-l, L-l+1] and each bottom series takes one unit horizontally
  plot.new()
  plot.window(xlim=c(0,m),ylim=c(0,L),xaxs="i",yaxs="i")
  if (!is.null(main)){
    title(main)
  }
  if (!is.null(ln)){
    # mtext ignores par("cex") unless told, so match strwidth above
    mtext(ln,side=2,at=L-(1:L)+0.5,las=1,line=0.3,adj=1,cex=par("cex"),font=level.font)
  }

  return(op)

}
