\name{summary.sbchoice}
\alias{summary.sbchoice}
\alias{print.summary.sbchoice}

\title{
Summarizing \code{sbchoice} estimation
}

\description{
Summary method for objects of class \code{sbchoice}.
}

\usage{
\method{summary}{sbchoice}(object, \dots)

\method{print}{summary.sbchoice}(x, digits = max(3, getOption("digits") - 1), \dots)

}

\arguments{
    \item{object}{an object of class \code{"sbchoice"}.}
    \item{x}{an object of class \code{"summary.sbchoice"}.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\details{
  The function \code{summary.sbchoice()} computes and returns a list of summary
  statistics of the fitted model in \code{object} of the \code{"sbchoice"} class. 

  Some of the values are printed up to certain decimal places. Actual values of 
  individual components are displayed separately, for instance, by \code{summary(object)$coefficients}. 
  See the \bold{Value} section for a list of components.
  
  Since the model for the single-bounded dichotomous choice CV data is estimated by \code{\link[stats]{glm}}, 
  an object of class\cr
  \code{"summary.sbchoice"} is constructed based on a \code{"summary.glm"} 
  class object. The summary of the \code{"summary.glm"} class object is available by 
  \code{summary(object)$glm.summary}. Other components computed by 
  \code{\link[stats]{summary.glm}} are also accessible. See \code{\link[stats]{summary.glm}} 
  for details.
}

\value{
  In addition to those available in the \code{object} of the \code{"sbchoice"} class, 
  the following list components are added.
    \item{glm.summary}{a summary of the \code{glm} estimation. It is an object of class 
        \code{"summary.glm"}.}
    \item{glm.null.summary}{a summary of the \code{glm} estimation of the null model 
        (i.e., only with the intercept). It is an object of class \code{"summary.glm"}.}
    \item{loglik}{the value of the log-likelihood of the model.}
    \item{loglik.null}{the value of the log-likelihood of the null model.}
    \item{psdR2}{McFadden's pseudo-R2 measure.}
    \item{adjpsdR2}{McFadden's pseudo-R2 measure adjusted for the degrees of freedom.}
    \item{medianWTP}{the estimated median WTP.}
    \item{meanWTP}{the estimated mean WTP.}
    \item{trunc.meanWTP}{the estimated mean WTP truncated at the maximum bid.}
    \item{adj.trunc.meanWTP}{the truncated mean WTP with the adjustment of 
    Boyle \emph{et~al.}(1988).}
    \item{LR.test}{a vector containing the likelihood ratio test statistic, the degrees 
        of freedom and the associated p-value.}
    \item{AIC}{information criterion (AIC and BIC).}
}

\references{
Boyle KJ, Welsh MP, Bishop RC (1988).
\dQuote{Validation of Empirical Measures of Welfare Change: Comment.}
\emph{Land Economics}, \bold{64}(1), 94--98.
}

\seealso{
\code{\link{sbchoice}}, \code{\link[stats]{glm}}, \code{\link[stats]{summary.glm}}
}

\keyword{DCchoice}
\keyword{single-bounded}
\keyword{methods}
