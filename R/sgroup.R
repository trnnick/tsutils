#' Grouped hierarchy S matrix from dimensions
#'
#' Construct the summing matrix S of a grouped hierarchy from the labels of the bottom level series on one or more dimensions, e.g., product type, region, and organic. Redundant aggregate nodes are removed.
#'
#' @param x data frame or matrix of labels, with one row per bottom level series and one column per dimension. Column names are used as dimension names.
#' @param names optional names of the bottom level series. If \code{NULL}, the row names of \code{x} are used, or else the labels of each series joined with \code{sep}.
#' @param sep separator used to construct the names of the nodes.
#' @param cross which combinations of dimensions to include, in addition to the total and each dimension on its own:
#' \itemize{
#' \item \code{"all"}: all combinations of two or more dimensions.
#' \item \code{"none"}: no combinations.
#' \item a vector of specific combinations, each with the dimension names separated by \code{":"}, e.g., \code{c("product:region","region:organic")}. A list of character vectors, e.g., \code{list(c("product","region"))}, is also accepted.
#' }
#' @param na how to treat missing labels:
#' \itemize{
#' \item \code{"exclude"}: a series with a missing label on a dimension is not part of the nodes of that dimension, but it is still part of the total and of the nodes of the other dimensions.
#' \item \code{"group"}: missing labels form their own group, named \code{"NA"}.
#' }
#' @param total name of the top level node.
#' @param sign optional vector of signs (\code{1} or \code{-1}), one per series, with which they enter aggregate nodes, e.g., \code{-1} for imports when net exports are needed. See details.
#' @param sort.levels if \code{TRUE}, then reorder the levels so that levels with fewer nodes are placed higher in the S matrix. Levels with the same number of nodes keep their original order. The bottom level is always last.
#' @param outplot if \code{TRUE}, then plot the hierarchy using \code{\link{plotS}}, with the names of the levels in the left margin.
#' @param ... additional arguments passed to \code{\link{plotS}}, when \code{outplot=TRUE}.
#'
#' @details
#' For the total, each dimension on its own, and each combination of dimensions requested in \code{cross}, the series are grouped by their labels, giving one aggregate node per combination of values. A combination of dimensions forms one level of the hierarchy. Levels are ordered from the total, to single dimensions, to combinations of two dimensions, and so on, with the bottom level series last.
#'
#' With \code{sign}, series enter aggregate nodes with their sign. Nodes whose series are all negative are gross totals and are kept positive, e.g., total imports, while nodes that mix signs are net values, e.g., net exports. Bottom level series are always positive.
#'
#' Aggregate nodes are redundant and removed when they:
#' \itemize{
#' \item aggregate a single series, as they repeat a bottom level series;
#' \item aggregate the same series as a node higher up, e.g., a dimension with a single value repeats the total, or a combination of nested dimensions (country within region) repeats the finer dimension. The first, higher, node is kept.
#' }
#' Levels left without any nodes disappear. If \code{sort.levels=TRUE}, the levels are reordered by their number of nodes, after removing the redundant nodes. Levels are named by their dimensions, e.g., \code{"product x region"}, and the bottom level is named \code{"Bottom"}. Nodes are named by their values joined with \code{sep}. If the same name arises from different dimensions, the dimension name is added, e.g., \code{"organic:yes"}.
#'
#' @return S matrix (\eqn{n \times m}), with row and column names. The attribute \code{"levels"} gives the level (combination of dimensions) of each row, which can be passed to \code{\link{plotS}}. The attribute \code{"removed"} gives the names of the removed redundant nodes.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{plotS}}, \code{\link{Saugment}}, \code{\link{Sthief}}.
#'
#' @keywords ts
#'
#' @examples
#' # Six products in five regions, organic or not
#' x <- expand.grid(organic=c("no","yes"),
#'                  region=c("denmark","sweden","norway","iceland","finland"),
#'                  product=c("milk","yoghurt","cheese","quark","whey","casein"),
#'                  stringsAsFactors=FALSE)[,3:1]
#' S <- Sgroup(x,outplot=TRUE)
#' dim(S)
#' table(attr(S,"levels"))
#'
#' # Only the total and each dimension on its own
#' S <- Sgroup(x,cross="none",outplot=TRUE)
#'
#' # Add only the product by region combination
#' S <- Sgroup(x,cross="product:region",outplot=TRUE)
#'
#' # Levels with fewer nodes higher up: organic, region, product, ...
#' S <- Sgroup(x,cross="product:region",sort.levels=TRUE,outplot=TRUE)
#' unique(attr(S,"levels"))
#'
#' # Block view, passing arguments to plotS
#' S <- Sgroup(x,cross="product:region",outplot=TRUE,view="block")
#'
#' # Redundant dimensions: market has a single value and union is nested in region
#' x$market <- "nordics"
#' x$union <- ifelse(x$region %in% c("norway","iceland"),"non-eu","eu")
#' S <- Sgroup(x)
#' head(attr(S,"removed"))
#'
#' # Missing labels
#' x <- x[,1:3]
#' x$region[1:4] <- NA
#' S <- Sgroup(x,cross="none")
#' S[c("Total","denmark"),1:6]
#'
#' # Net values: exports and imports of goods and services
#' x <- data.frame(flow=c("exports","exports","imports","imports"),
#'                 type=c("goods","services","goods","services"),
#'                 row.names=c("X.goods","X.services","M.goods","M.services"))
#' S <- Sgroup(x,cross="none",sign=ifelse(x$flow=="imports",-1,1),total="NetExports")
#' S
#' # NetExports, goods and services are net values, exports and imports are gross totals
#'
#' @export Sgroup

Sgroup <- function(x,names=NULL,sep="_",cross=c("all","none"),na=c("exclude","group"),
                   total="Total",sign=NULL,sort.levels=c(FALSE,TRUE),
                   outplot=c(FALSE,TRUE),...){

  # Defaults
  sort.levels <- sort.levels[1]
  outplot <- outplot[1]
  if (is.character(cross) && identical(cross,c("all","none"))){
    cross <- "all"
  }
  na <- match.arg(na,c("exclude","group"))
  x <- as.data.frame(x,stringsAsFactors=FALSE)
  m <- nrow(x)
  p <- ncol(x)
  if (m < 2){
    stop("Argument x must describe at least two series.")
  }
  if (is.null(colnames(x))){
    colnames(x) <- paste0("D",1:p)
  }
  dims <- colnames(x)
  for (j in 1:p){
    x[,j] <- as.character(x[,j])
    if (na == "group"){
      x[is.na(x[,j]),j] <- "NA"
    }
  }

  # Names of the bottom level series
  if (is.null(names)){
    if (!is.null(rownames(x)) && !all(rownames(x) == as.character(1:m))){
      names <- rownames(x)
    } else {
      names <- do.call(paste,c(x,sep=sep))
      if (any(duplicated(names)) || any(is.na(x))){
        names <- paste0("s",1:m)
      }
    }
  }
  if (length(names) != m || any(duplicated(names))){
    stop("Argument names must have one unique name per series.")
  }

  # Combinations of dimensions: the total, each dimension, and the requested combinations
  subs <- c(list(integer(0)),as.list(1:p))
  if (is.character(cross) && length(cross) == 1 && cross %in% c("all","none")){
    if (cross == "all" && p > 1){
      for (k in 2:p){
        subs <- c(subs,combn(p,k,simplify=FALSE))
      }
    }
  } else {
    if (is.character(cross)){
      cross <- strsplit(cross,":",fixed=TRUE)
    }
    if (!is.list(cross)){
      stop("Argument cross must be \"all\", \"none\", or a vector of combinations of dimensions.")
    }
    for (cmb in cross){
      cmb <- trimws(cmb)
      if (any(!(cmb %in% dims))){
        stop(paste0("Unknown dimension in cross: ",paste(setdiff(cmb,dims),collapse=", "),
                    ". Dimensions are: ",paste(dims,collapse=", "),"."))
      }
      if (length(unique(cmb)) < 2){
        stop("Each combination in cross must have at least two different dimensions.")
      }
      subs <- c(subs,list(sort(match(unique(cmb),dims))))
    }
    # Drop repeated combinations and order by number of dimensions
    subs <- subs[!duplicated(sapply(subs,paste,collapse=","))]
    subs <- subs[order(sapply(subs,length))]
  }

  # Aggregate nodes of each combination
  A <- list()
  A.names <- NULL
  A.full <- NULL
  A.level <- NULL
  for (s in subs){
    if (length(s) == 0){
      A[[length(A)+1]] <- rep(1,m)
      A.names <- c(A.names,total)
      A.full <- c(A.full,total)
      A.level <- c(A.level,total)
      next
    }
    key <- do.call(paste,c(x[,s,drop=FALSE],sep=sep))
    key[apply(is.na(x[,s,drop=FALSE]),1,any)] <- NA
    for (v in unique(key[!is.na(key)])){
      A[[length(A)+1]] <- as.numeric(!is.na(key) & key == v)
      A.names <- c(A.names,v)
      A.full <- c(A.full,paste0(paste(dims[s],collapse=":"),":",v))
      A.level <- c(A.level,paste(dims[s],collapse=" x "))
    }
  }
  A <- do.call(rbind,A)

  # Signs: series enter with their sign, unless all are negative (gross total)
  if (!is.null(sign)){
    if (length(sign) != m || any(is.na(sign)) || any(!(sign %in% c(-1,1)))){
      stop("Argument sign must have one entry per series, equal to -1 or 1.")
    }
    A <- A * matrix(sign,nrow(A),m,byrow=TRUE)
    neg <- apply(A,1,function(a){all(a[a != 0] < 0)})
    A[neg,] <- -A[neg,]
  }

  # Remove redundant nodes: single series, or same series as a node higher up
  rem <- rowSums(A != 0) <= 1 | duplicated(A)
  removed <- A.names[rem]
  A <- A[!rem,,drop=FALSE]
  A.names <- A.names[!rem]
  A.full <- A.full[!rem]
  A.level <- A.level[!rem]

  # Reorder levels by number of nodes, keeping the original order for ties
  if (sort.levels == TRUE){
    lev.names <- unique(A.level)
    lev.size <- table(factor(A.level,levels=lev.names))
    lev.ord <- lev.names[order(as.numeric(lev.size),1:length(lev.names))]
    idx <- order(match(A.level,lev.ord),1:length(A.level))
    A <- A[idx,,drop=FALSE]
    A.names <- A.names[idx]
    A.full <- A.full[idx]
    A.level <- A.level[idx]
  }

  # Unique names
  idx <- A.names %in% c(A.names[duplicated(A.names)],names)
  A.names[idx] <- A.full[idx]

  # S matrix
  S <- rbind(A,diag(m))
  dimnames(S) <- list(c(A.names,names),names)
  attr(S,"levels") <- c(A.level,rep("Bottom",m))
  attr(S,"removed") <- removed

  # Plot hierarchy, with the names of the levels
  if (outplot == TRUE){
    args <- list(...)
    if (is.null(args$levels)){
      args$levels <- attr(S,"levels")
    }
    if (is.null(args$level.labels)){
      args$level.labels <- unique(args$levels)
    }
    do.call(plotS,c(list(S=S),args),quote=TRUE)
  }

  return(S)

}
