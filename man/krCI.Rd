\name{krCI}
\alias{krCI}
\alias{print.krCI}


\title{Calculating confidence intervals for WTP using a parametric simulation}
\description{
    This function calculates confidence intervals for WTP using the method of Krinsky 
    and Robb (1986, 1990).
}
\usage{
krCI(obj, nsim = 1000, CI = 0.95)

\method{print}{krCI}(x, \dots)

}

\arguments{
    \item{obj}{an S3 class object \code{"dbchoice"} or \code{"sbchoice"}.}
    \item{nsim}{the number of draws of the parameters.}
    \item{CI}{a percentile of the confidence intervals to be estimated.}
    \item{x}{an object of class \code{"krCI"}.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\details{
    In the method of Krinsky and Robb (1986, 1990), a set of parameters is drawn \code{nsim} times 
    from a multivariate normal distribution with a vector of the parameter estimates as a mean 
    and the variance-covariance matrix of the parameter estimates. Then, various WTPs are computed for 
    each draw of simulated parameters. As a result, we are able to build an empirical distribution 
    of the WTPs concerned, and hence the confidence intervals. For each WTP, and when \code{nsim = 1000}, 
    the lower and the upper bound of the 95\% confidence interval (\code{CI = 0.95}) correspond to 
    the 26th and the 975th sorted estimates, respectively.

    Confidence intervals based on the bootstrap method are calculated by \code{\link{bootCI}}. 
    
    Hole (2007) conducted simulation experiments to compare the performance of the method of Krinsky 
    and Robb (1986, 1990) with the bootstrap one. 
}

\value{
  The function \code{krCI()} returns an object of S3 class \code{"krCI"}. 
  An object of \code{"krCI"} is a list with the following components.
  
  \item{out}{the output table with simulated confidence intervals as well as 
  the four type of WTP estimates (mean, truncated mean, truncated mean with adjustment 
  and median) from the ML estimation.}
  \item{mWTP}{a vector of simulated mean WTP. When |beta| < 1, this item is set 
  to \code{-999}.}
  \item{tr.mWTP}{a vector of simulated mean WTP truncated at the maximum bid.}
  \item{adj.tr.mWTP}{a vector of simulated mean WTP truncated at the maximum bid with the adjustment.}
  \item{medWTP}{a vector of simulated median WTP.}

  When the parameter estimate on the bid does not satisfy the condition for the existence 
  of the finite mean WTP (|beta|>1), the values of the lower and the upper bound of the 
  confidence interval are coerced to set to \code{-999}.

  The table contains the confidence intervals for the four types (mean, truncated mean, 
  truncated mean with adjustment and median) of WTP estimate from the ML estimation. 
  The adjustment for the truncated mean WTP is implemented by the method of Boyle 
  \emph{et~al.}(1988). 

  The generic function \code{print()} is available for the object of class 
  \code{"krCI"} and displays the table of simulated confidence intervals. 

}

\references{
Boyle KJ, Welsh MP, Bishop RC (1988).
\dQuote{Validation of Empirical Measures of Welfare Change: Comment.}
\emph{Land Economics}, \bold{64}(1), 94--98.

Hole AR (2007).
\dQuote{A Comparison of Approaches to Estimating Confidence Intervals for 
    Willingness to Pay Measure.}
\emph{Health Economics}, \bold{16}, 827--840.

Krinsky I, Robb AL (1986).
\dQuote{On Approximating the Statistical Properties of Elasticities.} 
\emph{The Review of Economics and Statistics}, \bold{68}, 715--719.

Krinsky I, Robb AL (1990).
\dQuote{On Approximating the Statistical Properties of Elasticities:
 A Correction.}
\emph{The Review of Economics and Statistics}, \bold{72}, 189--190.
}

\seealso{
\code{\link{bootCI}}, \code{\link{dbchoice}}, \code{\link{sbchoice}}
}

\examples{
## See Examples in dbchoice and sbchoice.
}

\keyword{DCchoice}
\keyword{double-bounded}
\keyword{single-bounded}
