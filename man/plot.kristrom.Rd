\name{plot.kristrom}
\alias{plot.kristrom}

\title{
Plotting \code{kristrom} objects
}

\description{
Plotting method for objects of class \code{"kristrom"}. The empirical survival curve 
is plotted.
}

\usage{
\method{plot}{kristrom}(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
    lwd = NULL, lty = NULL, \dots)
}

\arguments{
    \item{x}{an object of class \code{"kristrom"}.}
    \item{main}{the main title of the plot. If unspecified, no main title is displayed.}
    \item{sub}{the sub-title of the plot. If unspecified, no sub-title is displayed.}
    \item{xlab}{the x label of the plot. If missing, \code{xlab = "Bid"} is used. 
        Setting \code{xlab = ""} displays no x label.}
    \item{ylab}{the y label of the plot. If missing, \code{ylab = "Survival Probability"} is 
        used. Setting \code{ylab = ""} displays no y label.}
    \item{lwd}{the line width for the plot. If missing, \code{lwd = 3} is used.}
    \item{lty}{the line type for the plot. If missing, \code{lty = 1} is used.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\seealso{
    \code{\link{kristrom}}, \code{\link{summary.kristrom}}
}

\keyword{DCchoice}
\keyword{single-bounded}
\keyword{nonparametric}
\keyword{methods}
