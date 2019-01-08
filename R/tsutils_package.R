#' tsutils: Time Series Exploration, Modelling and Forecasting
#'
#' The tsutils package provides functions to support various aspects of time series and forecasting modelling. In particular this package includes: (i) tests and visualisations that can help the modeller explore time series components and perform decomposition; (ii) modelling shortcuts, such as functions to construct lagmatrices and seasonal dummy variables of various forms; (iii) an implementation of the Theta method; (iv) tools to facilitate the design of the forecasting process, such as ABC-XYZ analyses; and (v) "quality of life" tools, such as treating time series for trailing and leading values.
#'
#' @section Time series exploration:
#' \itemize{
#' \item{\code{\link{cmav}}}{: centred moving average.}
#' \item{\code{\link{coxstuart}}}{: Cox-Stuart test for location/dispersion.}
#' \item{\code{\link{decomp}}}{: classical time series decomposition.}
#' \item{\code{\link{seasplot}}}{: construct seasonal plots.}
#' \item{\code{\link{trendtest}}}{: test a time series for trend.}
#' }
#'
#' @section Time series modelling:
#' \itemize{
#' \item{\code{\link{getOptK}}}{: optimal temporal aggregation level for AR(1), MA(1), ARMA(1,1).}
#' \item{\code{\link{lagmatrix}}}{: create leads/lags of variable.}
#' \item{\code{\link{residout}}}{: construct control chart of residuals.}
#' \item{\code{\link{seasdummy}}}{: create seasonal dummies.}
#' \item{\code{\link{theta}}}{: Theta method.}
#' }
#'
#' @section Forecasting process modelling:
#' \itemize{
#' \item{\code{\link{abc}}}{: ABC analysis.}
#' \item{\code{\link{xyz}}}{: XYZ analysis.}
#' \item{\code{\link{abcxyz}}}{: ABC-XYZ analyses visualisation.}
#' }
#'
#' @section Quality of life:
#' \itemize{
#' \item{\code{\link{geomean}}}{: geometric mean.}
#' \item{\code{\link{lambdaseq}}}{: generate sequence of lambda for LASSO regression.}
#' \item{\code{\link{leadtrail}}}{: remove leading/training zeros/NAs.}
#' \item{\code{\link{wins}}}{: winsorisation, including vectorised versions \code{colWins} and \code{rowWins}.}
#' }
#'
#' @section Time series data:
#' \itemize{
#' \item{\code{\link{referrals}}}{: A&E monthly referrals.}
#' }
#'
#' @docType package
#' @keywords package
#'
#' @name tsutils
#'
#' @importFrom graphics plot lines par points polygon axis image abline box text boxplot legend mtext
#' @importFrom stats frequency friedman.test qtukey na.exclude sd ts pbinom median density quantile start cor.test deltat end optim time tsp
#' @importFrom forecast ets forecast.ets Arima forecast
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette gray rgb rainbow
#' @importFrom utils tail head
#' @importFrom MAPA tsaggr
#'
NULL
