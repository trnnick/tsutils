#' ABC-XYZ visualisation
#'
#' Jointly visualise ABC and XYZ analyses.
#'
#' @param imp an obkect of class \code{abc} that is the output of function \code{\link{abc}}.
#' @param frc an obkect of class \code{abc} that is the output of function \code{\link{xyz}}.
#' @param outplot if \code{TRUE}, then provide a visualisation of the analyses.
#' @param error vector of forecast errors for each series that will be distributed in each class, presented as an average.
#' @param ... additional arguments passed to the plot.
#'
#' @return A list containing:
#' \itemize{
#' \item{\code{class}}{: a matrix containing the number of time series in each class.}
#' \item{\code{error}}{: a matrix containing the averaged error for each class, if the argument \code{error} was used.}
#' }
#'
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' Ord K., Fildes R., Kourentzes N. (2017) \href{https://kourentzes.com/forecasting/2017/10/16/new-forecasting-book-principles-of-business-forecasting-2e/}{Principles of Business Forecasting, 2e}. \emph{Wessex Press Publishing Co.}, p.515-518.
#'
#' @seealso \code{\link{abc}}, \code{\link{xyz}}.
#'
#' @keywords ts
#'
#' @examples
#' x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
#' abcxyz(abc(x),xyz(x,type="cv"))
#'
#' @export abcxyz

abcxyz <- function(imp,frc,outplot=c(TRUE,FALSE),error=NULL,...){

  outplot <- outplot[1]

  k.abc <- length(imp$conc)
  nam.abc <- rownames(imp$conc)
  k.xyz <- length(frc$conc)
  nam.xyz <- rownames(frc$conc)

  n <- length(imp$value)
  matrix.abcxyz <- matrix.error <- matrix(0,k.abc,k.xyz,dimnames=list(nam.abc,nam.xyz))

  for (i in 1:k.abc){
    for (j in 1:k.xyz){
      matrix.abcxyz[i,j] <- sum(imp$class == nam.abc[i] & frc$class == nam.xyz[j])
      if (!is.null(error)){
        matrix.error[i,j] <- mean(error[imp$class == nam.abc[i] & frc$class == nam.xyz[j]])
      }
    }
  }

  if (outplot==TRUE){

    # Get colours
    cmp <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlGn")[2:10])(k.abc+k.xyz-1)
    # Create colour map matrix
    cmp.loc <- array(NA, c(1,1)*(k.abc+k.xyz-1))
    cmp.loc[,1] <- 1:(k.abc+k.xyz-1)
    for (i in 2:(k.abc+k.xyz-1)){
      cmp.loc[,i] <- c(tail(cmp.loc[,1], -(i-1)), head(cmp.loc[,1], (i-1)))
    }
    # cmp.loc <- cmp.loc[(k.abc):(k.abc+k.xyz-1),1:k.abc]
    cmp.loc <- cmp.loc[1:k.xyz,1:k.abc]

    # Allow user to override plot defaults
    args <- list(...)
    if (!("main" %in% names(args))){
      args$main <- "ABC-XYZ analyses"
    }
    if (!("xlab" %in% names(args))){
      args$xlab <- "Importance"
    }
    if (!("ylab" %in% names(args))){
      args$ylab <- "Forecastability"
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    if (!("yaxs" %in% names(args))){
      args$yaxs <- "i"
    }
    # Remaining defaults
    args$x <- args$y <- NA
    args$xlim <- c(0,k.abc)
    args$ylim <- c(0,k.xyz)
    args$xaxt <- args$yaxt <- "n"

    x.abc <- c(0,1:k.abc)
    y.xyz <- c(0,1:k.xyz)

    do.call(plot,args)

    for (i in 1:k.abc){
      for (j in 1:k.xyz){
        polygon(c(x.abc[i],x.abc[i+1],x.abc[i+1],x.abc[i]),
                c(y.xyz[j],y.xyz[j],y.xyz[j+1],y.xyz[j+1]),col=cmp[cmp.loc[j,i]])
        if (is.null(error)){
          text((x.abc[i]+x.abc[i+1])/2,(y.xyz[j]+y.xyz[j+1])/2,
               paste0(nam.abc[i],nam.xyz[j],":",matrix.abcxyz[i,j]))
        } else {
          text((x.abc[i]+x.abc[i+1])/2,(y.xyz[j]+y.xyz[j+1])/2,
               paste0(nam.abc[i],nam.xyz[j],":",matrix.abcxyz[i,j],"\n",
                 round(matrix.error[i,j],2)))
        }
      }
    }
    axis(side = 1, at = c(0.5, k.abc-0.5), labels = c("High","Low"))
    axis(side = 2, at = c(0.5, k.xyz-0.5), labels = c("Low","High"))
  }

  if (is.null(error)){
    matrix.error <- NULL
  }
  return(list("class"=matrix.abcxyz,"error"=matrix.error))

}
