\name{summary.turnbull}
\alias{summary.turnbull}
\alias{print.summary.turnbull}

\title{
Summarizing the Kaplan-Meier-Turnbull nonparametric estimation of WTP
}

\description{
Summary method for objects of class \code{"turnbull"}.
}

\usage{
\method{summary}{turnbull}(object, printCI = TRUE, \dots)

\method{print}{summary.turnbull}(x, digits = max(3, getOption("digits") - 1), \dots)
}

\arguments{
    \item{object}{an object of class \code{"turnbull"}.}
    \item{printCI}{an argument controlling whether the summary of confidence interval estimates are printed. The CIs are not summarized and printed unless \code{conf.int = TRUE} in \code{turnbull.sb} or \code{turnbull.db}. }
    \item{x}{an object of class \code{"turnbull"}.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\keyword{DCchoice}
\keyword{double-bounded}
\keyword{single-bounded}
\keyword{nonparametric}
\keyword{methods}
