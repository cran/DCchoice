\name{bootCI}
\alias{bootCI}
\alias{print.bootCI}


\title{Calculating confidence intervals for WTP using a nonparametric bootstrap method}
\description{
    This function calculates confidence intervals for WTP using the nonparametric bootstrap method.
}

\usage{
bootCI(obj, nboot = 1000, CI = 0.95)

\method{print}{bootCI}(x, \dots)

}

\arguments{
  \item{obj}{an S3 class object \code{"dbchoice"} or \code{"sbchoice"}.}
  \item{nboot}{the number of bootstrap resampling.}
  \item{CI}{a percentile of the confidence intervals to be estimated.}
  \item{x}{an object of class \code{"bootCI"}.}
  \item{\dots}{optional arguments. Currently not in use.}
}

\details{
    The bootstrap method resamples the data at our hands and repeatedly estimates the model 
    with the bootstrapped data to formulate an empirical distribution of the associated WTP. 
    This is a clear contrast with the method of Krinsky and Robb (1986, 1990) where the parameters are 
    directly drawn from the multivariate normal distribution.
    
    The upper and the lower bound of the interval is determined similarly to the case of 
    the function \code{\link{krCI}}.

    Hole (2007) conducted simulation experiments to compare the performance of the method of Krinsky 
    and Robb (1986, 1990) with the bootstrap one. 
}

\section{Warning}{
    It is time consuming (usually takes several minutes) to implement this function.
}

\value{
  The function \code{bootCI()} returns an object of S3 class \code{"bootCI"}. 
  An object of \code{"bootCI"} is a list with the following components.
  
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

  The generic function \code{print()} is available for the object of class 
  \code{"bootCI"} and displays the table of simulated confidence intervals. 

  The table contains the confidence intervals for the four types (mean, truncated mean, 
  truncated mean with adjustment and median) of WTP from the ML estimation. 
  The adjustment for the truncated mean WTP is implemented by the method of Boyle 
  \emph{et~al.}(1988). 

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
\code{\link{krCI}}, \code{\link{dbchoice}}, \code{\link{sbchoice}}
}

\examples{
## See Examples in dbchoice and sbchoice.
}


\keyword{DCchoice}\keyword{double-bounded}
\keyword{single-bounded}

