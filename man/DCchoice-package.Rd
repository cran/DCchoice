\name{DCchoice-package}
\alias{DCchoice-package}
\docType{package}

\title{DCchoice: a package for analyzing dichotomous choice contingent valuation data}
\description{
The package provides functions for analyzing single- and double-bounded dichotomous 
choice contingent valuation (CV) data
}

\details{
This package provides functions for analyzing single- and double-bounded dichotomous 
choice contingent valuation (CV) data.

In the single-bounded dichotomous choice CV that was first proposed by Bishop and Heberlein (1979) 
respondents are requested to answer a question like the following:

\strong{If the environmental policy burdens you with USD \eqn{X} annually, would you agree or disagree to it?}\cr

This format, respondents only states "yes (I agree)" or "no (I disagree)," meaning whether their 
willingness to pay (WTP) is greater or lower than the bid (USD \eqn{X}) they are offered.

The double-bounded dichotomous choice CV was proposed by Hanemann (1985) 
and Carson (1985) to improve the efficiency of single-bounded dichotomous choice CV. 
In the CV format, respondents are requested to answer the second (follow-up) question 
just after they answer the single-bounded dichotomous choice CV style question 
(the first/initial question). An example of double-bounded dichotomous choice CV 
questions is as follows (\eqn{Xl < X < Xh}):

First question\cr
\strong{If the environmental policy burdens you with USD \eqn{X} annually, would you agree or disagree to it?}

Second question for the respondents who agree to the policy in the first question\cr
\strong{If the amount of tax is USD \eqn{Xh}, would you agree or disagree to it?}

Second question for the respondents who disagree to the policy in the first question\cr
\strong{If the amount of tax is USD \eqn{Xl}, would you agree or disagree to it?}

In the double-bounded dichotomous choice CV question, there are four possible response 
outcomes: (yes, yes); (yes, no); (no, yes); and (no, no). If the respondent \eqn{i}'s answer 
is (yes, yes), the analyst can tell \eqn{WTPi > Xh} (\eqn{WTPi} is the WTP of the 
respondent \eqn{i}). Similarly, (yes, no) means \eqn{X < WTPi < Xh}, (no, yes) means 
\eqn{Xl < WTPi < X}, and (no, no) means \eqn{WTPi <Xl}.

There are two ways of estimating WTP from the single- and double-bounded dichotomous choice CV: 
parametric and nonparametric approaches. In this package, the functions \code{\link{sbchoice}} 
and \code{\link{dbchoice}}, which are based on the utility difference approach (Hanemann, 1984), 
are developed for the parametric approach to single- and double-bounded dichotomous choice data, 
respectively. 

Confidence intervals for the estimates of WTPs are constructed by two methods. 
These are the Krinsky and Robb (1986, 1990)'s method and the bootstrap one. 
The former is implemented by  \code{\link{krCI}} while the latter by \code{\link{bootCI}}. 

Both of the methods rely on simulation techniques with different settings. Usually, 
a bootstrap method takes much longer time than the Krinsky and Robb's method does. 
It has been pointed out that each method has both advantages and disadvantages, see, 
for instance, the discussions in Hole (2007) and the references therein. 

Functions for nonparametric approaches are also included in the package. \code{\link{kristrom}} 
(Kristrom, 1990) and \code{\link{turnbull.sb}} (Carson and Steinberg, 1990) are designed for 
analyses for single-bounded dichotomous choice data whereas \code{\link{turnbull.db}} 
(Carson and Hanemann, 2005) for double-bounded ones.

}


\note{
\pkg{DCchoice} depends on \pkg{interval} (Fay and Shaw, 2010) and \pkg{Ecdat} (Croissant, 2011) 
packages. The former further depends on several packages, among which \pkg{Icens} 
(Gentleman and Vandal, 2011) package may not be installed by the GUI Package Installer 
accessible from the menu bar (Windows and Mac OS) 
because it is available only from \href{http://bioconductor.org/}{Bioconductor}. 

To install \pkg{DCchoice} and other dependent packages simultaneously, use the 
\code{\link{install.packages}} function instead. See the \bold{Exapmles} section for the 
code. 

}



\references{
Bishop RC, Heberlein TA (1979).
\dQuote{Measuring Values of Extra-Market Goods: Are Indirect
  Measures Biased?}
\emph{American Journal of Agricultural Economics}, \bold{61}(5),
  926--930.

Carson RT (1985).
\dQuote{Three Essays on Contingent Valuation}.
Dissertation, University of California Berkeley.

Carson RT, Hanemann WM (2005).
\dQuote{Contingent Valuation.}
in KG M\"{a}ler, JR Vincent (eds.), \emph{Handbook of Environmental Economics}. 
Elsevier, New York.

Carson RT, Steinberg D (1990).
\dQuote{Experimental Design for Discrete Choice Voter Preference
  Surveys.}
in \emph{1989 Proceeding of the Survey Methodology Section of the
  American Statistical Association}, 821--822.

Croissant Y (2011).
\emph{\pkg{Ecdat}: Data Sets for Econometrics,} 
\R package version 0.1-6.1, 
\url{http://CRAN.R-project.org/package=Ecdat}.

Fay MP, Shaw PA (2010).
\dQuote{Exact and Asymptotic Weighted Logrank Tests for Interval Censored Data: 
    The \pkg{interval} \R Package}, 
\emph{Journal of Statistical Software}, 
\bold{36}(2), 1-34. \url{http://www.jstatsoft.org/v36/i02/}.

Gentleman R, Vandal A (2011).
\emph{Icens: NPMLE for Censored and Truncated Data}. 
\R package version 1.24.0, 
\url{http://CRAN.R-project.org/package=Icens}.

Hanemann, WM (1984). 
\dQuote{Welfare Evaluations in Contingent Valuation Experiments with Discrete Responses}, 
\emph{American Journal of Agricultural Economics},
\bold{66}(2), 332--341.

Hanemann WM (1985).
\dQuote{Some Issues in Continuous- and Discrete-Response Contingent
  Valuation Studies.}
\emph{Northeastern Journal of Agricultural Economics}, \bold{14},
  5--13.

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

Kristr\"{o}m B (1990). 
\dQuote{A Non-Parametric Approach to the Estimation of Welfare 
Measures in Discrete Response Valuation Studies.}
\emph{Land Economics}, \bold{66}(2), 135--139.

}

\examples{
## Installation of DCchoice along with dependent packages is carried out 
## by the following lines of code:
\dontrun{
install.packages("DCchoice", 
  repos = c("http://R-Forge.R-project.org", 
    "http://www.bioconductor.org/packages/release/bioc", 
    "http://cran.rstudio.com/"), 
 dep = TRUE)}
## For Mac and Unix/Linux users, add the option argument 'type="source"'. 
}



\keyword{DCchoice}
\keyword{package}

