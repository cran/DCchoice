\name{summary.kristrom}
\alias{summary.kristrom}
\alias{print.summary.kristrom}

\title{
Summarizing the Kristr\"{o}m's nonparametric estimation of WTP
}

\description{
Summary method for objects of class \code{"kristrom"}.
}

\usage{
\method{summary}{kristrom}(object, digits = max(3, getOption("digits") - 1), \dots)

\method{print}{summary.kristrom}(x, digits = max(3, getOption("digits") - 1), ...)
}

\arguments{
    \item{object}{an object of class \code{"kristrom"}.}
    \item{x}{an object of class \code{"kristrom"}.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\seealso{
    \code{\link{plot.kristrom}}, \code{\link{kristrom}}, 
    \code{\link[Ecdat]{NaturalPark}}, \code{\link{sbchoice}} 
}


\keyword{DCchoice}
\keyword{single-bounded}
\keyword{nonparametric}
\keyword{methods}
