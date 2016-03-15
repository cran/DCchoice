\name{plot.sbchoice}
\alias{plot.sbchoice}

\title{
Plotting sbchoice objects
}

\description{
Plotting method for objects of class \code{"sbchoice"}.
}

\usage{
\method{plot}{sbchoice}(x, type = NULL, main = NULL, sub = NULL, xlab = NULL, 
ylab = NULL, lwd = NULL, lty = NULL, xlim = NULL, ylim = NULL, bid = NULL, ...)
}

\arguments{
  \item{x}{
an object of class \code{"kristrom"}.
}
  \item{type}{
type of plot.
}
  \item{main}{
the main title of the plot. If unspecified, no main title is displayed.
}
  \item{sub}{
the sub-title of the plot. If unspecified, no sub-title is displayed.
}
  \item{xlab}{
the x label of the plot. If missing, \code{xlab = "Bid"} is used. 
Setting \code{xlab = ""} displays no x label.
}
  \item{ylab}{
the y label of the plot. If missing, \code{ylab = "Probability of selecting yes"} is 
used. Setting \code{ylab = ""} displays no y label.
}
  \item{lwd}{
the line width for the plot. If missing, \code{lwd = 3} is used.
}
  \item{lty}{
the line type for the plot. If missing, \code{lty = 1} is used.
}
  \item{xlim}{
the x limits of the plot.
}
  \item{ylim}{
the y limits of the plot.
}
  \item{bid}{
the bid limits should be drawn. If missing, the minimum and maximum values of the bid
variable(s) in the original dataset is used.
}
  \item{\dots}{
optional arguments.
}
}

\details{
The function \code{plot()} draws choice probabilities of yes according to 
the range of bid (covariates are set on average).
}

\seealso{
\code{\link{sbchoice}}
}

\examples{
## See Examples in sbchoice.
}

\keyword{DCchoice}
\keyword{single-bounded}
\keyword{methods}
