#' Augment hierarchy with discrepancy series
#'
#' Augments a hierarchy whose observed series are not coherent with additional bottom level series that capture the observed discrepancies. The augmented hierarchy is coherent, so that standard forecast reconciliation methods can be applied to it without modification.
#'
#' @param S summing matrix (\eqn{n \times m}). Non-zero entries indicate that the bottom level series in that column is aggregated by the node in that row. Rows of \code{S} correspond to columns of \code{y}.
#' @param y matrix or \code{ts} matrix (\eqn{T \times n}) of observed series. If both \code{rownames(S)} and \code{colnames(y)} exist, the columns of \code{y} are matched to the rows of \code{S} by name. Otherwise the columns of \code{y} must follow the rows of \code{S}. Missing values are not allowed.
#' @param tol absolute tolerance that separates genuine discrepancies from rounding. Set this to the rounding precision of the data, e.g., \code{0.01} for data published with three decimals. The default is only suitable for data without rounding.
#' @param nodes optional vector of names or row indices of aggregate nodes that must receive their own discrepancy series. See details.
#' @param outplot if \code{TRUE}, then plot the augmented hierarchy using \code{\link{plotS}}, highlighting the discrepancy series.
#' @param ... additional arguments passed to \code{\link{plotS}}, when \code{outplot=TRUE}.
#'
#' @details
#' Bottom level series are identified as the rows of \code{S} that are unit vectors. All other rows are aggregate nodes. For aggregate node \eqn{i} the discrepancy is measured against the bottom level, \eqn{\varepsilon_{i,t} = y_{i,t} - s_i' b_t}, where \eqn{s_i} is the row of \code{S} and \eqn{b_t} the bottom level series.
#'
#' Aggregate nodes are processed from the smallest to the largest, in terms of the number of bottom level series they aggregate, and are classified as:
#' \itemize{
#' \item \code{"coherent"}: the discrepancy is within \code{tol} of zero.
#' \item \code{"implied"}: the discrepancy is, within \code{tol}, a linear combination of the discrepancies of nodes that already have their own discrepancy series. For example, net exports defined as exports minus imports. The coefficients of the combination are stored in the corresponding row of \code{D}. Coefficients within \code{1e-3} of an integer are rounded to that integer.
#' \item \code{"own"}: the node receives its own discrepancy series.
#' }
#' When discrepancies are linearly dependent, which nodes receive their own discrepancy series is not unique. Use \code{nodes} to choose them.
#'
#' @return A list with:
#' \itemize{
#' \item \code{S}: augmented summing matrix (\eqn{(n+k) \times (m+k)}). Rows are the rows of the input \code{S}, in the same order, followed by the \eqn{k} discrepancy series. Columns are the bottom level series followed by the discrepancy series.
#' \item \code{y}: input series with the discrepancy series appended. This is a \code{ts} object if \code{y} is one.
#' \item \code{D}: matrix (\eqn{n_a \times k}) that maps the discrepancy series to the aggregate nodes.
#' \item \code{d}: the discrepancy series (\eqn{T \times k}).
#' \item \code{eps}: discrepancies of all aggregate nodes (\eqn{T \times n_a}).
#' \item \code{status}: classification of each aggregate node.
#' \item \code{error}: maximum absolute incoherence of the augmented series, i.e., \code{max(abs(y - B \%*\% t(S)))}, with \code{B} the augmented bottom level series. This should be within \code{tol}.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{plotS}}, \code{\link{Sthief}}.
#'
#' @keywords ts
#'
#' @examples
#' # Small hierarchy: Total = A + B, A = A1 + A2, B = B1 + B2
#' S <- rbind(c(1,1,1,1),c(1,1,0,0),c(0,0,1,1),diag(4))
#' rownames(S) <- c("Total","A","B","A1","A2","B1","B2")
#' set.seed(1)
#' b <- matrix(rnorm(40,10),10,4)
#' y <- b %*% t(S)
#' colnames(y) <- rownames(S)
#'
#' # Add discrepancies at Total and A
#' y[,"Total"] <- y[,"Total"] + rnorm(10)
#' y[,"A"] <- y[,"A"] + rnorm(10)
#' aug <- Saugment(S,ts(y,frequency=4),outplot=TRUE)
#' aug$status
#' aug$S
#'
#' # Block view, passing arguments to plotS
#' aug <- Saugment(S,y,outplot=TRUE,view="block")
#'
#' # Total is published as A + B, while A and B differ from their bottom level
#' # series. The discrepancy of Total is then the sum of those of A and B.
#' y <- b %*% t(S)
#' colnames(y) <- rownames(S)
#' y[,"A"] <- y[,"A"] + rnorm(10)
#' y[,"B"] <- y[,"B"] + rnorm(10)
#' y[,"Total"] <- y[,"A"] + y[,"B"]
#'
#' # By default A and B get their own discrepancy series and Total is implied
#' Saugment(S,y)$status
#'
#' # Use nodes to give Total its own discrepancy series instead. B is now implied
#' # as the difference between the discrepancies of Total and A
#' aug <- Saugment(S,y,nodes="Total")
#' aug$status
#' aug$D
#'
#' @export Saugment

Saugment <- function(S,y,tol=sqrt(.Machine$double.eps)*max(abs(y)),nodes=NULL,
                     outplot=c(FALSE,TRUE),...){

  # Defaults
  outplot <- outplot[1]
  S <- as.matrix(S)
  n <- nrow(S)
  m <- ncol(S)
  if (is.null(rownames(S))){
    rownames(S) <- paste0("y",1:n)
  }
  if (is.null(colnames(S))){
    colnames(S) <- paste0("b",1:m)
  }

  # Match y to the rows of S
  if (any(class(y) == "ts")){
    y.tsp <- tsp(y)
  } else {
    y.tsp <- NULL
  }
  y <- unclass(as.matrix(y))
  attr(y,"tsp") <- NULL
  if (!is.null(colnames(y)) && all(rownames(S) %in% colnames(y))){
    y <- y[,rownames(S),drop=FALSE]
  } else if (ncol(y) == n){
    colnames(y) <- rownames(S)
  } else {
    stop("Columns of y must match the rows of S, either by name or by position.")
  }
  if (any(is.na(y))){
    stop("Argument y contains missing values. Trim the series to a common sample first.")
  }

  # Bottom level series are the unit rows of S, the rest are aggregate nodes
  isUnit <- apply(S,1,function(s){sum(s != 0) == 1 && sum(s) == 1})
  idxB <- sapply(1:m,function(j){which(isUnit & S[,j] == 1)[1]})
  if (any(is.na(idxB))){
    stop("Every column of S needs a bottom level series, i.e., a unit row.")
  }
  idxA <- setdiff(1:n,idxB)
  if (length(idxA) == 0){
    stop("S has no aggregate nodes.")
  }

  # Discrepancies against the bottom level
  B <- y[,idxB,drop=FALSE]
  eps <- y[,idxA,drop=FALSE] - B %*% t(S[idxA,,drop=FALSE])
  colnames(eps) <- rownames(S)[idxA]

  # Classify aggregate nodes, from the smallest to the largest
  if (is.character(nodes)){
    if (any(!(nodes %in% rownames(S)[idxA]))){
      stop("Argument nodes must refer to aggregate nodes of S.")
    }
    nodes <- match(nodes,rownames(S))
  }
  idxOrd <- idxA[order(rowSums(S[idxA,,drop=FALSE] != 0))]
  if (!is.null(nodes)){
    idxOrd <- c(intersect(idxOrd,nodes),setdiff(idxOrd,nodes))
  }
  status <- rep("",length(idxA))
  names(status) <- rownames(S)[idxA]
  idxOwn <- NULL
  coefs <- list()
  for (i in idxOrd){
    nam <- rownames(S)[i]
    e <- eps[,nam]
    # Forced to have own discrepancy series
    if (!is.null(nodes) && (i %in% nodes)){
      status[nam] <- "own"
      idxOwn <- c(idxOwn,i)
      next
    }
    # Coherent
    if (max(abs(e)) <= tol){
      status[nam] <- "coherent"
      next
    }
    # Implied by the discrepancies of other nodes
    if (length(idxOwn) > 0){
      X <- eps[,rownames(S)[idxOwn],drop=FALSE]
      cf <- qr.coef(qr(X),e)
      cf[is.na(cf)] <- 0
      cf <- ifelse(abs(cf - round(cf)) < 1e-3,round(cf),cf)
      names(cf) <- paste0("d.",rownames(S)[idxOwn])
      if (max(abs(e - X %*% cf)) <= tol){
        status[nam] <- "implied"
        coefs[[nam]] <- cf
        next
      }
    }
    status[nam] <- "own"
    idxOwn <- c(idxOwn,i)
  }

  # Coherent hierarchy, nothing to augment
  k <- length(idxOwn)
  if (k == 0){
    message("No discrepancies found. The hierarchy is coherent within tol.")
    if (!is.null(y.tsp)){
      y <- ts(y,start=y.tsp[1],frequency=y.tsp[3])
    }
    return(list(S=S,y=y,D=NULL,d=NULL,eps=eps,status=status,error=0))
  }

  # Discrepancy series and their mapping to the aggregate nodes
  d.names <- paste0("d.",rownames(S)[idxOwn])
  d <- eps[,rownames(S)[idxOwn],drop=FALSE]
  colnames(d) <- d.names
  D <- matrix(0,n,k,dimnames=list(rownames(S),d.names))
  for (j in 1:k){
    D[idxOwn[j],j] <- 1
  }
  for (nam in names(coefs)){
    D[nam,names(coefs[[nam]])] <- coefs[[nam]]
  }

  # Augmented S and series
  S.aug <- rbind(cbind(S,D),cbind(matrix(0,k,m),diag(k)))
  dimnames(S.aug) <- list(c(rownames(S),d.names),c(colnames(S),d.names))
  y.aug <- cbind(y,d)
  err <- max(abs(y.aug - cbind(B,d) %*% t(S.aug)))
  if (!is.null(y.tsp)){
    y.aug <- ts(y.aug,start=y.tsp[1],frequency=y.tsp[3])
    d <- ts(d,start=y.tsp[1],frequency=y.tsp[3])
  }

  # Plot augmented hierarchy, highlighting the discrepancy series
  if (outplot == TRUE){
    args <- list(...)
    if (is.null(args$fill)){
      args$fill <- c(rep("grey88",n),rep(brewer.pal(3,"Pastel1")[2],k))
    }
    do.call(plotS,c(list(S=S.aug),args),quote=TRUE)
  }

  return(list(S=S.aug,y=y.aug,D=D[idxA,,drop=FALSE],d=d,eps=eps,status=status,error=err))

}
