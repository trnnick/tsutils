% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seasdummy.R
\name{seasdummy}
\alias{seasdummy}
\title{Create seasonal dummy variables.}
\usage{
seasdummy(n, m = NULL, y = NULL, type = c("bin", "trg"), full = c(FALSE, TRUE))
}
\arguments{
\item{n}{number of observations to create.}

\item{m}{seasonal periodicity. If \code{NULL} it will take the information from the provided time series (\code{y} argument). See notes.}

\item{y}{this is an optional time series input that can be used to get seasonal periodicity (\code{m}) and the start point.}

\item{type}{type of seasonal dummies to create.
\itemize{
\item \code{"bin}: binary dummies.
\item \code{"trg}: trigonometric dummies. See notes.
}}

\item{full}{If full is \code{TRUE}, then keeps the m-th dummy that is co-linear to the rest. See notes.}
}
\value{
An array with seasonal dummies, where rows correspond observations and columns to dummy variables.
}
\description{
Create binary or trigonometric seasonal dummies.
}
\note{
If the seasonal periodicity is fractional then the the type will be overriden to trigonometric and only two seasonal dummies with be produced. One cosine and one sine.
}
\examples{
seasdummy(24,12)

}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
