\name{plot.turnbull}
\alias{plot.turnbull}

\title{
Plotting \code{turnbull} objects
}

\description{
Plotting method for objects of class \code{"turnbull"}. The empirical survival curve and confidence interval (if computed) are plotted.
}

\usage{
\method{plot}{turnbull}(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
        lwd = NULL, lty = NULL, plotCI = FALSE, ltyCI = 5, \dots)
}

\arguments{
    \item{x}{an object of class \code{"turnbull"}.}
    \item{main}{the main title of the plot. If unspecified, no main title is displayed.}
    \item{sub}{the sub-title of the plot. If unspecified, no sub-title is displayed.}
    \item{xlab}{the x label of the plot. If missing, \code{xlab = "Bid"} is used. 
        Setting \code{xlab = ""} displays no x label.}
    \item{ylab}{the y label of the plot. If missing, \code{ylab = "Survival Probability"} is 
        used. Setting \code{ylab = ""} displays no y label.}
    \item{lwd}{the line width for the plot. If missing, \code{lwd = 3} is used.}
    \item{lty}{the line type for the plot. If missing, \code{lty = 1} is used.}
    \item{plotCI}{logical. If \code{TRUE} and \code{x} contains the estimates of the confidence intervals, these are plotted along with the survival function.}
    \item{ltyCI}{a graphical parameter defining the line type of the confidence interval. By default, \code{ltyCI = 5} (dashed line).}
    \item{\dots}{optional arguments. Currently not in use.}
}

\seealso{
    \code{\link{turnbull.sb}}, \code{\link{turnbull.db}}, \code{\link{summary.turnbull}}
}

\keyword{DCchoice}
\keyword{double-bounded}
\keyword{single-bounded}
\keyword{nonparametric}
\keyword{methods}
