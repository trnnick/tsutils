#' Mean aggregation S matrix
#'
#' Converts a summing matrix \code{S}, where aggregate nodes are sums of bottom level series, to one where aggregate nodes are means of the bottom level series they aggregate. The means can use equal weights or a vector of bottom level weights. Negative entries, e.g., imports in net exports, are allowed.
#'
#' @param S summing matrix (\eqn{n \times m}). Negative entries indicate series that enter a node with a negative sign, e.g., from \code{\link{Sgroup}} with argument \code{sign}.
#' @param w optional vector of non-negative bottom level weights, one per column of \code{S}. If named, it is matched to the column names of \code{S}. If \code{NULL}, all bottom level series have equal weight, giving simple means.
#'
#' @details
#' Each aggregate row \eqn{i} of \code{S} is replaced by \eqn{s_{ij} w_j / \sum_k |s_{ik}| w_k}, so that every aggregate node is a convex combination of the bottom level series it aggregates, entering with their sign. For example, net exports with equal weights become \eqn{(X_g + X_s - M_g - M_s)/4}. For non-negative \code{S} the normalisation is the usual weighted sum. The weights do not need to sum to one, as they are normalised within each node. Bottom level rows (unit rows) are not changed, so a bottom level series with zero weight is excluded from the means but remains in the hierarchy. Entries of \code{S} other than one, e.g., from \code{\link{Sgroup}} or temporal hierarchies, are kept as relative weights.
#'
#' The result is still a linear aggregation, \eqn{y = S b}, so it can be used in place of the summing matrix for forecast reconciliation, e.g., for series that are rates, prices, or stocks rather than flows.
#'
#' If structural scaling is used with the resulting \code{S}, it should be derived from its assumptions rather than computed from \code{S}: with uncorrelated bottom level errors of equal variance, the error variance of node \eqn{i} is proportional to \eqn{\sum_j s_{ij}^2 w_j^2 / (\sum_j |s_{ij}| w_j)^2}, which is \eqn{1/k} for an equal weight mean over \eqn{k} series. Using \code{diag(rowSums(S))}, as for sums, gives the identity matrix, which corresponds to OLS rather than structural scaling.
#'
#' @return S matrix with weighted mean aggregation, with the same dimensions, names and attributes as \code{S}.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{Sgroup}}, \code{\link{Sthief}}, \code{\link{Saugment}}.
#'
#' @keywords ts
#'
#' @examples
#' # Small hierarchy: A = H1 + H2, H1 = Q1 + Q2, H2 = Q3 + Q4
#' S <- rbind(c(1,1,1,1),c(1,1,0,0),c(0,0,1,1),diag(4))
#' rownames(S) <- c("A","H1","H2","Q1","Q2","Q3","Q4")
#'
#' # Equal weight means
#' Smean(S)
#'
#' # Weighted means, e.g., by sales volume of each bottom level series
#' Smean(S,w=c(10,30,20,40))
#'
#' # Temporal hierarchy for quarterly data where the annual value is the average
#' Smean(Sthief(4))
#'
#' # Net exports: exports and imports of goods and services
#' S <- rbind(NetExports=c(1,1,-1,-1),Exports=c(1,1,0,0),Imports=c(0,0,1,1),diag(4))
#' colnames(S) <- rownames(S)[4:7] <- c("X.goods","X.services","M.goods","M.services")
#' Smean(S)
#'
#' # Structural scaling for the mean hierarchy, derived from its assumptions
#' Sm <- Smean(S)
#' W <- diag(rowSums(Sm^2)/rowSums(abs(Sm))^2)
#' diag(W)
#'
#' @export Smean

Smean <- function(S,w=NULL){

  # Defaults
  S.in <- S
  S <- as.matrix(S)
  m <- ncol(S)

  # Bottom level weights
  if (is.null(w)){
    w <- rep(1,m)
  } else {
    if (length(w) != m){
      stop(paste0("Argument w must have one weight per column of S (",m,")."))
    }
    if (!is.null(names(w)) && !is.null(colnames(S))){
      if (!all(colnames(S) %in% names(w))){
        stop("Names of w must match the column names of S.")
      }
      w <- w[colnames(S)]
    }
    if (any(is.na(w)) || any(w < 0)){
      stop("Argument w must contain non-negative weights.")
    }
  }

  # Weighted S, normalised within each aggregate row
  isUnit <- apply(S,1,function(s){sum(s != 0) == 1 && sum(s) == 1})
  Sw <- S * matrix(w,nrow(S),m,byrow=TRUE)
  tot <- rowSums(abs(Sw))
  bad <- which(!isUnit & tot == 0)
  if (length(bad) > 0){
    if (!is.null(rownames(S))){
      bad <- rownames(S)[bad]
    }
    stop(paste0("Aggregate nodes with zero total weight: ",paste(head(bad,10),collapse=", "),
                ifelse(length(bad) > 10,", ...",""),". Check w."))
  }
  Sw[!isUnit,] <- Sw[!isUnit,,drop=FALSE]/tot[!isUnit]
  Sw[isUnit,] <- S[isUnit,,drop=FALSE]

  # Keep names and other attributes of S, e.g., levels from Sgroup
  if (is.matrix(S.in)){
    S.in[] <- Sw
    Sw <- S.in
  }

  return(Sw)

}
