#' XYZ analysis
#'
#' Perform XYZ analysis on a set of time series.
#'
#' @param x this can either be an array, where each column is a series, or a vector of values. If \code{x} is a vector of values forecastability is not calculated and the input is used as such.
#' @param m seasonal length for time series. Required when type is \code{"naive"} or \code{"ets"}.
#' @param prc a vector of percentages indicating how many items are included in each class. By default this is \code{c(0.2,0.3,0.5)}, but any set of percentage values can be used as long as \code{0<=prc[i]<=1} and \code{sum(prc)==1}.
#' @param type the type of forecastability calculation. This can be:
#' \itemize{
#' \item{\code{"naive"}}{: fit naive and seasonal naive and calculate forecastability using RMSE/mean level.}
#' \item{\code{"ets"}}{: fit ets and calculate and calculate forecastability using RMSE/mean level.}
#' \item{\code{"cv"}}{: use coefficient of variation as a proxy of forecastability.}
#' }
#'
#'
#' @return Return object of class \code{abc} and contains:
#' \itemize{
#' \item{\code{value}}{: a vector containing the forecastability value of each series.}
#' \item{\code{class}}{: a vector containing the class membership of each series.}
#' \item{\code{rank}}{: a vector containing the rank of each series, with 1 being the lowest forecastability series.}
#' \item{\code{conc}}{: the forecastability concentration of each class, as percentage of total value.}
#' \item{\code{model}}{: fitted model for each series.}
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @seealso \code{\link{abc}}, \code{\link{plot.abc}}, \code{\link{abcxyz}}.
#'
#' @references
#' Ord K., Fildes R., Kourentzes N. (2017) \href{http://kourentzes.com/forecasting/2017/10/16/new-forecasting-book-principles-of-business-forecasting-2e/}{Principles of Business Forecasting, 2e}. \emph{Wessex Press Publishing Co.}, p.515-518.
#'
#' @keywords ts
#'
#' @examples
#' x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
#' z <- xyz(x,m=12)
#' print(z)
#' plot(z)
#'
#' @export xyz

xyz <- function(x,m=NULL,prc=c(0.2,0.3,0.5),type=c("naive","ets","cv")){

    type <- type[1]

    n <- dim(x)[2]              # Number of series total
    if (is.null(n)){            # If x is not an array, then get length
        n <- length(x)
    }

    if (sum(dim(x)==1)>0 | any(class(x)=="numeric")){
        x.mean <- x
        x.model <- NULL
    } else {
        x.mean <- array(0,c(1,n))
        x.model <- array(NA,c(1,n))
        for (i in 1:n){
            if (type=="cv"){
                x.mean[i] <- sd(x[,i],na.rm=TRUE)/mean(x[,i],na.rm=TRUE)
                x.model[i] <- "cv"
            } else {
                if (is.null(m)){
                    stop("Seasonality m not provided")
                }
                x.temp <- x[,i]
                x.temp <- x.temp[!is.na(x.temp)]
                if (type=="ets"){
                    fit <- ets(ts(x.temp,frequency=m))
                    x.mean[i] <- sqrt(fit$mse)/mean(x.temp)
                    x.model[i] <- fit$method
                } else {
                    l <- length(x.temp)
                    n1.fit <- c(NA,x.temp[1:(l-1)])
                    nm.fit <- c(rep(NA,m),x.temp[1:(l-m)])
                    n.err <- c(mean((x.temp[(m+1):l] - n1.fit[(m+1):l])^2),
                               mean((x.temp[(m+1):l] - nm.fit[(m+1):l])^2))
                    cmp.err <- (l-m)*log(n.err) + c(2,2*m)
                    cmp <- cmp.err == min(cmp.err)
                    x.mean[i] <- sqrt(n.err[cmp])/mean(x.temp)
                    x.model[i] <- c("Naive","Seasonal Naive")[cmp]
                }
            }
        }
    }
    # Find rank and percentage contribution of each series
    x.rank <- order(x.mean, decreasing=TRUE)
    x.sort <- x.mean[x.rank]
    x.sort <- (x.sort/sum(x.sort))*100

    k <- length(prc)            # Number of classes
    p <- array(0,c(k,1))        # Number of series in each class
    x.ind <- array(k,c(1,n))    # Indicator for class of each series
    x.class <- array(NA,c(1,n)) # Class of each series
    nam.xyz <- LETTERS[26:(26-k+1)]
    x.imp <- array(0,c(k,1),dimnames=list(nam.xyz,"Errors"))    # Percentage importance of each class

    # Calculate classes
    for (i in 1:(k)){
        p[i] <- round(n*prc[i])
        if (i==1){
            x.ind[x.rank[1:p[i]]] <- i
            x.imp[i] <- sum(x.sort[1:sum(p[1:i])])
        } else if (i!=k) {
            x.ind[x.rank[(sum(p[1:(i-1)])+1):sum(p[1:i])]] <- i
            x.imp[i] <- sum(x.sort[1:sum(p[1:i])]) - sum(x.imp[1:(i-1)])
        } else {
            p[i] <- n - sum(p[1:(i-1)])
            x.imp[i] <- sum(x.sort[1:sum(p[1:i])]) - sum(x.imp[1:(i-1)])
        }
        x.class[x.ind==i] <- nam.xyz[i]
    }

    names(prc) <- nam.xyz

    return(structure(list("type"="XYZ","prc"=prc,"value"=x.mean,"class"=x.class,"rank"=x.rank,"conc"=x.imp,"model"=x.model),class="abc"))

}
