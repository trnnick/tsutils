#' Plot temporal hierarchy
#'
#' Plots a temporal hierarchy, either as a tree or as blocks, and names its levels descriptively, e.g., "Monthly", "2-monthly", "Quarterly", "Yearly", instead of by their temporal aggregation level \eqn{k}. The plots are produced by \code{\link{plotS}}, which calls \code{\link{plotStree}} or \code{\link{plotSblock}}.
#'
#' @param y input time series (a \code{ts} object), an integer, or a temporal hierarchy summing matrix. For a time series or an integer the summing matrix is constructed with \code{\link{Sthief}}, using its frequency or the integer as the seasonal periodicity.
#' @param freq frequency of the bottom level series, see details. If \code{NULL}, this is inferred from the size of the top level.
#' @param k optional vector of temporal aggregation levels to include, passed to \code{\link{Sthief}}. Ignored if \code{y} is a summing matrix.
#' @param view type of plot: \code{"tree"} (see \code{\link{plotStree}}) or \code{"block"} (see \code{\link{plotSblock}}). The block view is cleaner for large hierarchies, such as daily data within a year.
#' @param labels \code{FALSE} for no labels in the nodes, \code{TRUE} to number the nodes within each level, or one label per row of the summing matrix.
#' @param level.labels \code{TRUE} for the descriptive names of the levels, \code{FALSE} for none, or one label per level, starting from the top.
#' @param show.k if \code{TRUE}, the temporal aggregation level is added to the level names, e.g., "Quarterly (k = 3)".
#' @param shape node shape for the tree view: \code{"circle"} or \code{"box"}. Ignored for the block view.
#' @param node.size size of nodes without labels for the tree view, in inches. This is larger than the default of \code{\link{plotStree}}. Nodes are not made larger than the space available in the plot. Ignored for the block view.
#' @param fill fill colour of the nodes. Either a single colour or one per row of the summing matrix.
#' @param border border colour of the nodes.
#' @param ... additional arguments passed to \code{\link{plotStree}} or \code{\link{plotSblock}}, depending on \code{view}. For example, \code{edges}, \code{skip}, \code{level.font} or \code{main}.
#'
#' @details
#' Rows of the summing matrix are grouped into levels by the number of bottom level periods they aggregate (\eqn{k}), with the largest \eqn{k} at the top. Therefore, the order of the rows does not matter.
#'
#' The level names depend on the frequency of the bottom level. If \code{freq=NULL}, this is inferred from the number of bottom level periods in the top level, assuming that this is the usual seasonal cycle:
#' \itemize{
#' \item 4, 12, 13, 52: quarterly, monthly, 4-weekly and weekly data, with yearly top level.
#' \item 260, 364, 365: business daily and daily data, with yearly top level.
#' \item 5, 7, 168, 336: business daily, daily, hourly and half-hourly data, with weekly top level.
#' \item 24, 48, 96, 144, 288, 1440: hourly, half-hourly, 15-minute, 10-minute, 5-minute and minutely data, with daily top level.
#' }
#' For other cases, or when the hierarchy does not end in the usual cycle (e.g., 24 months with a 2-yearly top level, which would otherwise be taken as hourly data), use \code{freq}. This can be \code{"secondly"}, \code{"minutely"}, \code{"quarter-hourly"}, \code{"half-hourly"}, \code{"hourly"}, \code{"daily"}, \code{"business-daily"}, \code{"weekly"}, \code{"fortnightly"}, \code{"half-monthly"}, \code{"monthly"}, \code{"quarterly"}, \code{"half-yearly"} or \code{"yearly"}, optionally with a multiple, e.g., \code{"15-minute"}, \code{"2-hourly"} or \code{"4-weekly"}. Case, spaces, hyphens and underscores are ignored and nouns such as \code{"minute"} or \code{"month"} work as well.
#'
#' Each level is named using the largest time unit that divides it exactly, e.g., 120 minutes is "2-hourly" and 14 days is "2-weekly". The following conventions are used: 7 days or 5 business days are a week; 13, 26 and 52 weeks are "Quarterly", "Half-yearly" and "Yearly"; and 365 days are a year. Levels that do not fit any of these keep the unit of the bottom level, e.g., "26-daily".
#'
#' @return Invisibly, a list with the output of \code{\link{plotStree}} (\code{nodes} and \code{edges}) or of \code{\link{plotSblock}} (\code{blocks}), with row indices referring to the rows of the summing matrix as given, and:
#' \itemize{
#' \item \code{levels}: data frame with the \code{level}, its temporal aggregation level \code{k} and its \code{name}.
#' \item \code{freq}: name of the frequency of the bottom level.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
#'
#' @references
#' Athanasopoulos, G., Hyndman, R. J., Kourentzes, N., & Petropoulos, F. (2017). \href{https://kourentzes.com/forecasting/2017/02/27/forecasting-with-temporal-hierarchies-3/}{Forecasting with temporal hierarchies}. European Journal of Operational Research, 262(1), 60-74.
#'
#' @seealso \code{\link{Sthief}}, \code{\link{plotStree}}, \code{\link{plotSblock}}, \code{\link{plotS}}.
#'
#' @keywords ts
#'
#' @examples
#' # Full temporal hierarchies
#' plotThief(AirPassengers)
#' plotThief(4,show.k=TRUE)
#' plotThief(24)
#' plotThief(7,labels=TRUE)
#'
#' # Plain level labels and smaller nodes
#' plotThief(12,level.font=1,node.size=0.25)
#'
#' # Same structure, but half-monthly data
#' plotThief(24,freq="half-monthly")
#'
#' # Hourly data aggregated to days and weeks, but not years
#' plotThief(168,k=c(168,24,12))
#'
#' # Two years of monthly data
#' plotThief(24,k=c(24,12,6),freq="monthly")
#'
#' # 15-minute data, with numbered boxes and elbow edges
#' plotThief(4,freq="15-minute",labels=TRUE,shape="box",edges="elbow")
#'
#' # Block view, better for large hierarchies
#' plotThief(12,view="block")
#' plotThief(364,view="block",main="Daily data")
#'
#' @export plotThief

plotThief <- function(y,freq=NULL,k=NULL,view=c("tree","block"),labels=FALSE,
                      level.labels=TRUE,show.k=c(FALSE,TRUE),shape=c("circle","box"),
                      node.size=0.5,
                      fill=RColorBrewer::brewer.pal(3,"Pastel1")[2],
                      border=RColorBrewer::brewer.pal(3,"Set1")[2],...){

  # Defaults
  view <- match.arg(view,c("tree","block"))
  show.k <- show.k[1]
  shape <- match.arg(shape,c("circle","box"))
  args <- list(...)
  if (any(c("levels","parse") %in% names(args))){
    stop("Arguments levels and parse are set by plotThief.")
  }

  # Get S matrix
  # Time series or integers go through Sthief, which also checks the input
  if (any(class(y) == "ts") || is.null(dim(y))){
    S <- Sthief(y,k=k)
  } else {
    S <- as.matrix(y)
    if (!is.null(k)){
      warning("Argument k is ignored when y is a summing matrix.")
    }
  }
  n <- nrow(S)
  M <- S != 0
  kn <- rowSums(M)                    # Temporal aggregation level of each row
  if (any(kn == 0)){
    stop(paste0("Rows of the summing matrix with no bottom level series: ",
                paste(which(kn == 0),collapse=", "),"."))
  }
  node.names <- rownames(S)
  if (is.null(node.names)){
    node.names <- paste0("row",1:n)
  }

  # Order rows by temporal aggregation level, and within each level by position
  first <- apply(M,1,function(x){which(x)[1]})
  ord <- order(-kn,first)
  kn.ord <- kn[ord]
  K <- unique(kn.ord)
  for (i in K){
    if (any(colSums(M[kn==i,,drop=FALSE]) > 1)){
      warning(paste0("Nodes aggregating ",i," periods overlap. This is not a standard temporal hierarchy."))
    }
  }

  # Find frequency of the bottom level and name the levels
  if (is.null(freq)){
    freq <- thief.freq.infer(max(kn))
  }
  base <- thief.freq.parse(freq)
  freq.name <- thief.period.name(base$unit,base$n)
  lev.names <- sapply(K,function(x){thief.period.name(base$unit,base$n*x)})
  if (isTRUE(level.labels)){
    ll <- lev.names
    if (show.k == TRUE){
      ll <- paste0(ll," (k = ",K,")")
    }
  } else if (isFALSE(level.labels)){
    ll <- FALSE
  } else {
    if (length(level.labels) != length(K)){
      stop(paste0("Argument level.labels must have one entry per level (",length(K)," levels found)."))
    }
    ll <- level.labels
  }

  # Node labels and colours, reordered as the rows
  if (isTRUE(labels)){
    lab <- as.character(ave(1:n,kn.ord,FUN=seq_along))
  } else if (isFALSE(labels) || is.null(labels)){
    lab <- FALSE
  } else {
    if (length(labels) != n){
      stop(paste0("Argument labels must be TRUE, FALSE, or have one entry per row of the summing matrix (",n,")."))
    }
    lab <- labels[ord]
  }
  fill <- rep_len(fill,n)[ord]
  S.ord <- S[ord,,drop=FALSE]
  rownames(S.ord) <- node.names[ord]

  # Plot
  plot.args <- c(list(S=S.ord,view=view,labels=lab,parse=FALSE,levels=kn.ord,fill=fill,
                      border=border,level.labels=ll),args)
  if (view == "tree"){
    plot.args$shape <- shape
    plot.args$node.size <- node.size
  }
  out <- do.call(plotS,plot.args,quote=TRUE)

  # Map rows back to the original order of the summing matrix
  if (view == "tree"){
    out$nodes$row <- ord[out$nodes$row]
    out$edges$from.row <- ord[out$edges$from.row]
    out$edges$to.row <- ord[out$edges$to.row]
  } else {
    out$row <- ord[out$row]
    out <- list(blocks=out)
  }
  out$levels <- data.frame(level=1:length(K),k=K,name=lev.names,stringsAsFactors=FALSE)
  out$freq <- freq.name

  return(invisible(out))

}

# Internal function for plotThief.
# Frequency of the bottom level implied by the number of periods in the top level.
thief.freq.infer <- function(k){

  tab <- c("4"="quarterly","5"="business-daily","7"="daily","12"="monthly","13"="4-weekly",
           "24"="hourly","48"="half-hourly","52"="weekly","96"="15-minute","144"="10-minute",
           "168"="hourly","260"="business-daily","288"="5-minute","336"="half-hourly",
           "364"="daily","365"="daily","1440"="minutely")
  out <- tab[as.character(k)]
  if (is.na(out)){
    stop(paste0("Cannot infer the frequency of the data from a top level of ",k,
                " periods. Use argument freq, e.g., freq=\"monthly\"."))
  }

  return(unname(out))

}

# Internal function for plotThief.
# Splits a frequency, such as "monthly", "15-minute" or "half-hourly", into a time unit and a multiple of it.
thief.freq.parse <- function(freq){

  if (!is.character(freq) || length(freq) != 1){
    stop("Argument freq must be a single string, e.g., \"monthly\".")
  }
  f <- tolower(gsub("[[:space:]_-]","",freq))

  # Get the multiple, if any
  mlt <- 1
  num <- regmatches(f,regexec("^([0-9]+)(.+)$",f))[[1]]
  if (length(num) == 3){
    mlt <- as.integer(num[2])
    f <- num[3]
  }

  # Units and some special cases
  units <- list(second=c("second","seconds","secondly","sec","secs","s"),
                minute=c("minute","minutes","minutely","min","mins"),
                hour=c("hour","hours","hourly","hr","hrs","h"),
                day=c("day","days","daily","d"),
                bday=c("businessday","businessdays","businessdaily","bday","bdaily",
                       "weekday","weekdays","workingday","workday"),
                week=c("week","weeks","weekly","wk","w"),
                halfmonth=c("halfmonth","halfmonthly","semimonthly"),
                month=c("month","months","monthly","mon"),
                year=c("year","years","yearly","annual","annually","y"))
  special <- list(halfhourly=c("minute",30),quarterhourly=c("minute",15),
                  fortnightly=c("week",2),fortnight=c("week",2),
                  quarterly=c("month",3),quarter=c("month",3),q=c("month",3),
                  halfyearly=c("month",6),halfyear=c("month",6),semiannual=c("month",6),
                  semiannually=c("month",6),biannual=c("month",6))
  if (f %in% names(special)){
    return(list(unit=special[[f]][1],n=as.numeric(special[[f]][2])*mlt))
  }
  unit <- names(units)[sapply(units,function(x){f %in% x})]
  if (length(unit) == 0){
    stop(paste0("Unknown frequency \"",freq,"\". Use for example \"monthly\", \"hourly\", \"15-minute\" or \"weekly\"."))
  }
  # Years are counted in months, so that they work with quarters and half-years
  if (unit == "year"){
    return(list(unit="month",n=12*mlt))
  }

  return(list(unit=unit,n=mlt))

}

# Internal function for plotThief.
# Descriptive name of a period of n time units, using the largest unit that divides it exactly.
thief.period.name <- function(unit,n){

  # "hourly" or "2-hourly", etc.
  ly <- function(n,word){
    if (n == 1){
      paste0(toupper(substring(word,1,1)),substring(word,2))
    } else {
      paste0(n,"-",word)
    }
  }

  out <- switch(unit,
                second = if (n %% 60 == 0){thief.period.name("minute",n/60)
                } else if (n == 1){"Secondly"
                } else {paste0(n,"-second")},
                minute = if (n %% 60 == 0){thief.period.name("hour",n/60)
                } else if (n == 30){"Half-hourly"
                } else if (n == 15){"Quarter-hourly"
                } else if (n == 1){"Minutely"
                } else {paste0(n,"-minute")},
                hour = if (n %% 24 == 0){thief.period.name("day",n/24)
                } else {ly(n,"hourly")},
                day = if (n %% 7 == 0){thief.period.name("week",n/7)
                } else if (n %% 365 == 0){thief.period.name("year",n/365)
                } else {ly(n,"daily")},
                bday = if (n %% 5 == 0){thief.period.name("week",n/5)
                } else if (n == 1){"Business-daily"
                } else {paste0(n,"-business-day")},
                week = if (n %% 52 == 0){thief.period.name("year",n/52)
                } else if (n == 13){"Quarterly"
                } else if (n == 26){"Half-yearly"
                } else {ly(n,"weekly")},
                halfmonth = if (n %% 2 == 0){thief.period.name("month",n/2)
                } else if (n == 1){"Half-monthly"
                } else {paste0(n/2,"-monthly")},
                month = if (n %% 12 == 0){thief.period.name("year",n/12)
                } else if (n == 3){"Quarterly"
                } else if (n == 6){"Half-yearly"
                } else {ly(n,"monthly")},
                year = ly(n,"yearly"),
                stop(paste0("Unknown time unit: ",unit,".")))

  return(out)

}
