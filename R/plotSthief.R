#' Plot temporal hierarchy
#'
#' Plots the temporal hierarchy for a given time series of seasonal periodicity.
#'
#' @param y input time series (a \code{ts} object) or an integer.
#' @param labels if \code{TRUE} labels will be added for the temporal aggregation levels if the seasonal period is 4 (quarters), 7 (days in a week), 12 (months), 24 (hours), 48 (half-hours), 52 (weeks) or 364 (days).
#' @param ... additional arguments passed to the plotting function.
#'
#' @return Produces a plot of the temporal hierarchy.
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' Athanasopoulos, G., Hyndman, R. J., Kourentzes, N., & Petropoulos, F. (2017). \href{http://kourentzes.com/forecasting/2017/02/27/forecasting-with-temporal-hierarchies-3/}{Forecasting with temporal hierarchies}. European Journal of Operational Research, 262(1), 60-74.
#'
#' @keywords ts
#'
#' @examples
#' plotSthief(AirPassengers)
#'
#' @export plotSthief

plotSthief <- function(y, labels=c(TRUE,FALSE), ...){

  labels <- labels[1]

  # Get S matrix
  # This also checks the nature of the input
  S <- Sthief(y)
  n <- dim(S)[1]
  m <- dim(S)[2] # Frequency

  # Check ellipsis for plotting arguments
  args <- list(...)
  args.nms <- names(args)

  # Solf defaults for plot
  if (!("bty" %in% args.nms)){
    args$bty <- "n"
  }

  if (!("xlab" %in% args.nms)){
    args$xlab <- ""
  }

  if (!("ylab" %in% args.nms)){
    args$ylab <- ""
  }

  # Hard defaults for plot
  args$x <- args$y <- NA
  args$xaxt <- args$yaxt <- "n"
  args$xlim <- args$ylim <- c(0,1)

  # Create plot map
  # I need to find all nodes that are made up by lower levels
  Cons <- vector("list",n)
  bottom <- which(apply(S,1,sum) == 1)
  for (i in n:1){
    zerocon <- S[i,] == 0
    # Find relevant nodes
    # First remove those that are non-zero when this is zero
    use <- apply((S == 0)[,zerocon]==TRUE,1,function(x){all(x==TRUE)})
    # Then make sure they have at least one non-zero
    use <- use & apply((S == 1)[,!zerocon,drop=FALSE]==TRUE,1,function(x){any(x==TRUE)})
    use[i] <- FALSE       # Exclude self
    # Make sure that the number of elements can actually make up this level (i.e., multiples of lower make up upper)
    loc <- which(use)
    loc <- loc[sum(S[i,]) %% apply(S[use,,drop=FALSE],1,sum) == 0]
    # Check if these nodes are already used up and are not bottom level
    loc <- setdiff(loc,setdiff(unique(unlist(Cons[n:(i+1)])),bottom))
    # Check if bottom level and have been used up
    locB <- unique(unlist(Cons[setdiff(loc,bottom)]))
    if (setequal(intersect(loc,bottom),locB)){
      loc <- setdiff(loc,locB)
    }
    # That's it!
    Cons[[i]] <- loc
  }
  if (n > (length(bottom)+1)){
    Cons[[1]] <- setdiff(Cons[[1]],bottom)
  }

  # Initialiase plotting
  Sn <- apply(S,1,sum)        # Get number of elements per level
  Kn <- unique(Sn)            # Find unique composition sizes
  N <- length(Kn)             # Which leads to size of hierarchy
  Ln <- sapply(Kn,function(x){sum(Sn==x)})
  # Find locations
  mrg <- 0.2                  # Margins
  vloc <- seq(mrg/2,1-mrg/2,length.out=N)
  hloc <- vector("list",N)
  for (i in 1:N){
    hloc[[i]] <- seq(mrg/2,1-mrg/2,length.out=Ln[i]+2)
  }

  # Start plotting
  do.call(plot,args)

  # Plot connection lines
  cmp <- colorRampPalette(RColorBrewer::brewer.pal(9,"Greys")[9:3])(N-1)
  for (i in (N-1):1){
    # Find nodes that belong to this level
    loc <- which(Sn == Kn[N-i])
    for (j in 1:length(loc)){
      nd <- Cons[[loc[j]]]
      for (l in 1:length(nd)){
        # Find layer
        locV <- which(Sn[nd[l]]==Kn)
        locH <- which(nd[l] == which(Sn == Sn[nd[l]]))
        lines(c(hloc[[N-i]][j+1],hloc[[locV]][locH+1]),c(vloc[i+1],vloc[N-locV+1]),lwd=2,col=cmp[i])
      }
    }
  }

  # Plot nodes
  rd <- max(0.015,1/((max(Ln)+2)*2.5))
  rd <- min(rd,0.04)
  for (i in 1:N){
    for (j in 1:Ln[i]){
      plotrix::draw.circle(hloc[[i]][j+1],vloc[N+1-i],rd,col="#377EB8")
    }
  }

  # Get labels
  if (labels==TRUE && any(m == c(4, 7, 12, 24, 48, 52, 364))){
    switch(as.character(m),
           "4" = {lbls <- c("Quarter","Semi-annual","Annual")},
           "7" = {lbls <- c("Daily","Weekly")},
           "12" = {lbls <- c("Monthly","2-monthly","Quarterly","4-monthly","Half-yearly","Yearly")},
           "24" = {lbls <- c("Hourly","2-hourly","3-hourly","4-hourly","Quarter-daily","8-hourly","Half-daily","Daily")},
           "48" = {lbls <- c("Half-hourly","Hourly","1.5-hourly","2-hourly","3-hourly","4-hourly","Quarter-daily","8-hourly","Half-daily","Daily")},
           "52" = {lbls <- c("Weekly","2-weekly","4-weekly","Quarterly","Half-yearly","Yearly")},
           "364" = {lbls <- c("Daily","2-daily","4-daily","Weekly","13-daily","2-weekly","26-daily","52-daily","Quarterly","Half-yearly","Yearly")})

    for (i in 1:N){
      text(mrg/2,vloc[i],lbls[i],xpd=TRUE,pos=2)
    }
  }

}

