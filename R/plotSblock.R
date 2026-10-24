#' Plot hierarchy as blocks
#'
#' Plots a hierarchy from its summing matrix \code{S}, with each node (row of \code{S}) shown as the block of bottom level series it aggregates. Levels are arranged from the top to the bottom of the hierarchy. The horizontal axis has one unit per bottom level series (columns of \code{S}). This view remains clean for large hierarchies. For the more traditional view use \code{\link{plotStree}}.
#'
#' @param S summing matrix (\eqn{n \times m}). Non-zero entries indicate that the bottom level series in that column is aggregated by the node in that row. This can describe cross-sectional, grouped or temporal hierarchies.
#' @param labels \code{TRUE} to use the row names of \code{S} as block labels (no labels if \code{S} has no row names), \code{FALSE} for no labels, or one label per row of \code{S}. Labels can be a character vector, an expression vector, or a list mixing strings, quoted calls and expressions.
#' @param parse if \code{TRUE} and \code{labels=TRUE}, row names are parsed as \link[grDevices]{plotmath}, e.g. \code{"H[1]"} is plotted as H with subscript 1. Names that are not valid plotmath are plotted as text, with a warning.
#' @param skip how to treat the space left by a node that skips a level:
#' \itemize{
#' \item \code{"empty"}: leave it empty.
#' \item \code{"expand"}: extend the block upwards through every level above that covers none of its columns. This stops at the first level that covers any of them.
#' }
#' @param levels optional vector with the level of each row of \code{S}. Levels are ordered top to bottom in order of appearance. If \code{NULL}, levels are found automatically. See details.
#' @param fill fill colour of the blocks. Either a single colour or one per row of \code{S}.
#' @param border border colour of the blocks.
#' @param hgap horizontal space between blocks, as a fraction of a column width.
#' @param vgap vertical space between levels, as a fraction of the height of a level.
#' @param cex maximum label size.
#' @param min.cex labels that would need to be smaller than this to fit in their block are not plotted.
#' @param level.labels \code{FALSE} for none, \code{TRUE} for \code{"L1", "L2", ...}, or one label per level, starting from the top, as a character or expression vector. These are plotted in the left margin.
#' @param level.font font of the level labels, as in \code{par("font")}: 1 plain, 2 bold, 3 italic, 4 bold italic.
#' @param main optional plot title.
#'
#' @details
#' By default, rows are assigned to levels by scanning \code{S} from the top and starting a new level whenever a row overlaps with a row already in the current level. Therefore, rows within a level aggregate different bottom level series. This works with the usual ordering of \code{S} for cross-sectional, grouped and temporal hierarchies. For unusual row orders, provide the level of each row with \code{levels}.
#'
#' In grouped hierarchies some nodes aggregate bottom level series that are not next to each other (e.g., a product node in a region x product grouping). These are plotted as several blocks, all with the same label and colour.
#'
#' @return Invisibly, a data frame with one row per block and columns:
#' \itemize{
#' \item \code{node}: row name of \code{S}.
#' \item \code{row}: row index in \code{S}.
#' \item \code{level}: level of the node.
#' \item \code{top}: highest level reached, after any expansion.
#' \item \code{first}, \code{last}: first and last column of \code{S} covered by the block.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' For the design and explanation of this plot see \href{https://arxiv.org/abs/2609.23223}{arXiv:2609.23223}.
#'
#' @seealso \code{\link{plotStree}}, \code{\link{plotS}}, \code{\link{plotThief}}.
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
#' plotSblock(S1)
#'
#' # Row names parsed as plotmath and level labels
#' rownames(S1) <- c("A","H[1]","H[2]","Q[1]","Q[2]","Q[3]","Q[4]")
#' plotSblock(S1,parse=TRUE,level.labels=expression(bold(S)[4],bold(S)[2],bold(S)[1]))
#'
#' # Labels mixing text and expressions
#' plotSblock(S1,labels=list("Total",quote(H[1]),quote(H[2]),"Q1","Q2","Q3",quote(Q[4]^"*")))
#'
#' # Q5 skips a level and connects directly to A
#' S2 <- makeS(list(A=c("Q1","Q2","Q3","Q4","Q5"),H1=c("Q1","Q2"),H2=c("Q3","Q4")),
#'             c("Q1","Q2","Q3","Q4","Q5"))
#' plotSblock(S2,skip="empty")
#' plotSblock(S2,skip="expand")
#'
#' # Temporal hierarchy for monthly data
#' plotSblock(Sthief(12),labels=FALSE,level.labels=paste("k =",c(12,6,4,3,2,1)))
#'
#' # Grouped hierarchy: region x product
#' bottom <- c("R1P1","R1P2","R2P1","R2P2")
#' S3 <- makeS(list(Total=bottom,R1=c("R1P1","R1P2"),R2=c("R2P1","R2P2"),
#'                  P1=c("R1P1","R2P1"),P2=c("R1P2","R2P2")),bottom)
#' cmp <- RColorBrewer::brewer.pal(5,"Pastel1")
#' plotSblock(S3,fill=c("grey88",cmp[2],cmp[2],cmp[5],cmp[5],rep("grey88",4)))
#'
#' # Larger hierarchy without labels
#' b <- paste0("b",1:60)
#' S4 <- makeS(c(list(Total=b),
#'               setNames(lapply(0:4,function(g){b[g*12+1:12]}),paste0("G",1:5)),
#'               setNames(lapply(0:19,function(g){b[g*3+1:3]}),paste0("SG",1:20))),b)
#' plotSblock(S4,labels=FALSE)
#'
#' @export plotSblock

plotSblock <- function(S,labels=TRUE,parse=c(FALSE,TRUE),skip=c("empty","expand"),
                       levels=NULL,fill="grey88",border="grey30",hgap=0.08,vgap=0.35,
                       cex=1,min.cex=0.35,level.labels=FALSE,level.font=1,main=NULL){

  # Defaults
  parse <- parse[1]
  skip <- match.arg(skip,c("empty","expand"))

  # Get structure, labels and levels
  prep <- plotS.prep(S,labels,parse,levels)
  M <- prep$M
  n <- prep$n
  m <- prep$m
  L <- prep$L
  lev <- prep$lev
  fill <- rep_len(fill,n)
  for (l in 1:L){
    if (any(colSums(M[lev==l,,drop=FALSE]) > 1)){
      warning(paste0("Level ",l," contains overlapping nodes. Blocks will overlap."))
    }
  }

  # Find blocks
  # Each block is a run of adjacent bottom series of a row
  blocks <- vector("list",n)
  for (i in setdiff(1:n,prep$empty)){
    r <- rle(M[i,])
    last <- cumsum(r$lengths)
    first <- last - r$lengths + 1
    blocks[[i]] <- data.frame(node=prep$node.names[i],row=i,level=lev[i],top=lev[i],
                              first=first[r$values],last=last[r$values],
                              stringsAsFactors=FALSE)
  }
  blocks <- do.call(rbind,blocks)
  rownames(blocks) <- NULL

  # Expand blocks upwards into empty space if requested
  if (skip == "expand"){
    for (b in 1:nrow(blocks)){
      cols <- blocks$first[b]:blocks$last[b]
      k <- blocks$level[b] - 1
      while (k >= 1 && !any(prep$coverage[k,cols])){
        blocks$top[b] <- k
        k <- k - 1
      }
    }
  }

  # Start plotting
  op <- plotS.setup(L,m,level.labels,level.font,main)
  on.exit(par(op))

  # Plot blocks
  x0 <- blocks$first - 1 + hgap/2
  x1 <- blocks$last - hgap/2
  y0 <- L - blocks$level + vgap/2
  y1 <- L - blocks$top + 1 - vgap/2
  rect(x0,y0,x1,y1,col=fill[blocks$row],border=border)

  # Add labels, shrinking them to fit
  for (b in 1:nrow(blocks)){
    lab <- prep$lab.list[[blocks$row[b]]]
    if (is.null(lab)){
      next
    }
    fit <- min(0.85*(x1[b]-x0[b])/strwidth(lab,cex=1),0.7*(y1[b]-y0[b])/strheight(lab,cex=1))
    cex.b <- min(cex,fit)
    if (cex.b >= min.cex){
      text((x0[b]+x1[b])/2,(y0[b]+y1[b])/2,lab,cex=cex.b)
    }
  }

  return(invisible(blocks))

}
