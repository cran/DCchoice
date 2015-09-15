\name{summary.dbchoice}
\alias{summary.dbchoice}
\alias{print.summary.dbchoice}

\title{
Summarizing \code{dbchoice} estimation
}

\description{
Summary method for objects of class \code{"dbchoice"}.
}

\usage{
\method{summary}{dbchoice}(object, \dots)

\method{print}{summary.dbchoice}(x, digits = max(3, getOption("digits") - 1), \dots)
}

\arguments{
    \item{object}{an object of class \code{"dbchoice"}.}
    \item{x}{an object of class \code{"summary.dbchoice"}.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\details{
  The function \code{summary.dbchoice()} computes and returns a list of summary
  statistics of the fitted model in \code{object} of the \code{"dbchoice"} class. 

  Some of the values are printed up to certain decimal places. Actual values of 
  individual components are displayed separately, for instance, by \code{summary(object)$coefficients}. 
  See the \bold{Value} section for a list of components.
  
}


\value{
  In addition to those available in the \code{object} of the \code{"dbchoice"} class, 
  the following list components are added.
    \item{medianWTP}{the estimated median WTP.}
    \item{meanWTP}{the estimated mean WTP.}
    \item{trunc.meanWTP}{the estimated mean WTP truncated at the maximum bid.}
    \item{adj.trunc.meanWTP}{the truncated mean WTP with the adjustment of 
    Boyle \emph{et~al.}(1988).}
    \item{LR.test}{a vector containing the likelihood ratio test statistic, the degrees 
        of freedom and the associated p-value. The null is that all the parameters on the explanatory variables other than constant and the bid variable are jointly zero.}
    \item{coef}{a table of estimates including their s.e., z-value, and p-value.}
    \item{AIC}{information criterion (AIC and BIC).}
}

\references{
Boyle KJ, Welsh MP, Bishop RC (1988).
\dQuote{Validation of Empirical Measures of Welfare Change: Comment.}
\emph{Land Economics}, \bold{64}(1), 94--98.
}

\seealso{
\code{\link{dbchoice}}
}

\keyword{DCchoice}
\keyword{double-bounded}
\keyword{methods}
