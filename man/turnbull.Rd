\name{turnbull}
\alias{turnbull.db}
\alias{turnbull.sb}
\alias{turnbull.oohb}
\alias{print.turnbull}

\title{The Kaplan-Meier-Turnbull nonparametric approach to analyze dichotomous choice contingent valuation data}
\description{
This function analyzes single-, one-and-one-half-, or double-bounded dichotomous choice contingent valuation (CV) 
data on the basis of the Kaplan-Meier-Turnbull method.
}

\usage{
## for the single-bounded data
turnbull.sb(formula, data, subset, conf.int = FALSE, B = 200, 
            conf.level = 0.95, timeMessage = FALSE, \dots)

## for the one-and-one-half-bound data
turnbull.oohb(formula, data, subset, conf.int = FALSE, B = 200, 
              conf.level = 0.95, timeMessage = FALSE, \dots)


## for the double-bounded data
turnbull.db(formula, data, subset, conf.int = FALSE, B = 200, 
            conf.level = 0.95, timeMessage = FALSE, \dots)

\method{print}{turnbull}(x, digits = max(3, getOption("digits") - 1), \dots)
}

\arguments{
    \item{formula}{a formula specifying the model structure.}
    \item{data}{a data frame containing the variables in the model formula.}
    \item{subset}{an optional vector specifying a subset of observations.}
    \item{x}{an object of class \code{"turnbull"}.}
    \item{conf.int}{logical. If \code{TRUE}, the confidence intervals are computed.}
    \item{B}{number of bootstrap resampling for confidence interval estimation. Defaulted to 200.}
    \item{conf.level}{a confidence coefficient. Defaulted to 0.95.}
    \item{timeMessage}{logical. If \code{TRUE}, the estimated time for bootstrapping is displayed.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\details{
The functions \code{turnbull.sb()}, \code{turnbull.oohb()}, and \code{turnbull.db()} analyze single-,  
one-and-one-half-, and double-bounded dichotomous choice contingent valuation (CV) data, respectively, on 
the basis of the modified Kaplan-Meier-Turnbull method (Carson and Steinberg, 1990).

\strong{For single-bounded dichotomous choice data}\cr
Most of the details of \code{turnbull.sb()} is the same as those of \code{turnbull.db()}. 
See the subsequent section for details. 

A difference between \code{turnbull.sb()} and \code{turnbull.db()} appears in the 
definition of the model formula. In \code{turnbull.sb()}, the first part contains 
only one response variable (e.g., \code{R1}) and the other part contains only one 
bid variable (e.g., \code{BD1}) because respondents are requested to answer a CV 
question in the single-bounded dichotomous choice CV. A typical structure of the 
formula is given by

\code{R1 ~ BD1}

\code{\link{kristrom}} is an alternative nonparametric method for analyzing single-bounded 
dichotomous choice data. A parametric analysis can be done by \code{\link{sbchoice}}.

\strong{For one-and-one-half-bound dichotomous choice data}\cr
The details of \code{turnbull.oohb()} is the same as those of \code{turnbull.db()}.
See the the subsequent section for details.

A difference between \code{turnbull.oohb()} and \code{turnbull.db()} appears in the
definition of variables. See \code{\link{oohbchoice}} and \code{\link{oohbsyn}} 
for the details of creating response and bid variables.  

\strong{For double-bounded dichotomous choice data}\cr
A generic call to \code{turnbull.db()} is given by

\code{turnbull.db(formula, data)}
  
The argument \code{formula} defines the response variables and bid variables. The argument 
\code{data} is set as a data frame containing the variables in the model. 

A typical structure of the formula for \code{turnbull.db()} is defined as follows:

\code{R1 + R2 ~ BD1 + BD2}

The formula consists of two parts. The first part, the left-hand side of the 
tilde sign (\code{~}), must contain the response variables for the suggested 
prices in the first and the second stage of CV questions. In the above example, 
\code{R1} denotes a binary or two-level factor response variable for a bid in 
the first stage and \code{R2} for a bid in the second stage. Each of \code{R1} 
and \code{R2} contains \code{"Yes"} or \code{"No"} to the bid or \code{1} for 
\code{"Yes"} and \code{0} for \code{"No"}. 
The two variables are connected with the arithmetic operator (\code{+}). The other 
part, which starts after the tilde sign, must contain bid variables (\code{BD1}, 
\code{BD2}) containing suggested prices in the first and second stage of double-bounded 
dichotomous choice CV question. The two variables are also connected with the 
arithmetic operator (\code{+}).

A parametric approach for analyzing double-bounded dichotomous choice data can be 
carried out by \code{\link{dbchoice}}. 

The structure of data set which assigned to the argument \code{data} is the same as that 
of \code{dbchoice()}. See \code{\link{dbchoice}} for details in the data set structure.

The confidence intervals are computed by the bootstrap method in \code{\link[interval]{icfit}} of the \pkg{interval} (Fay and Shaw, 2010) package. The arguments \code{conf.int}, \code{B}, \code{conf.level} and \code{timeMessage} are passed to \code{icfit()}. The bootstrap can be time consuming, so that it is in general not advisable to increase the number of bootstrap resamplint \code{B}. See the help of \code{\link[interval]{icfit}} for futher detail.

The generic function \code{print()} is available for a fitted model object 
of class \code{"turnbull"} and displays simply estimated probabilities of 
the distribution: this is the same as the element \code{pf} in the function 
\code{\link[interval]{icfit}} of the \pkg{interval} (Fay and Shaw, 2010) package.

The extractor function \code{summary()} is also available for a fitted model 
object of class \code{"turnbull"} and the generic function \code{print.summary()} 
displays the survival probability and three types of WTP estimates (a Kaplan-Meier mean, 
a Spearman-Karber mean, and median WTPs). In estimating the two types of mean WTP, 
the area under the survivor function is truncated at the maximum bid because 
there seems no unified approach to determine an ending point of bids at which 
the acceptance probability is zero. Therefore, we leave the decision of how 
the area is treated to the user. In practice, once the ending point, 
$bid_{end}$, is found, it is straightforward to compute the triangular 
area by $0.5(bid_{end} - bid_{max})P_{max}$ where $bid_{max}$ is 
the maximum bid and $P_{max}$ is the acceptance probability for 
$bid_{max}$, both of which are reported in the summarized output.


Furthermore, the generic function \code{plot()} is available for a fitted model object 
of class \code{"turnbull"} and displays the empirical survival function. See 
\code{\link{plot.turnbull}} for details.
}

\value{
Both \code{turnbull.db()} and \code{turnbull.sb()} return an object of S3 class 
\code{"turnbull"} that is a list with the following components.
\item{left}{a vector of left endpoints of censored interval. The vector is 
    internally assigned to the argument \code{L} in \code{icfit()} of the 
    \pkg{interval} package.}
\item{right}{a vector of right endpoints of censored interval. The vector is 
    internally assigned to the argument \code{R} in \code{icfit()} of the 
    \pkg{interval} package.}
\item{turnbull}{a list of components returned from \code{icfit()}.}
\item{unq.dib}{a vector of unique bids including Inf.}
\item{CI}{estimates for confidence intervals from \code{icfit()}.}
}



\references{
Croissant Y (2011).
\emph{\pkg{Ecdat}: Data Sets for Econometrics,} 
\R package version 0.1-6.1, 
\url{http://CRAN.R-project.org/package=Ecdat}.

Fay MP, Shaw PA (2010).
\dQuote{Exact and Asymptotic Weighted Logrank Tests for Interval Censored Data: 
    The \pkg{interval} \R Package}, 
\emph{Journal of Statistical Software}, 
\bold{36}(2), 1-34. \url{http://www.jstatsoft.org/v36/i02/}.

Krinsky I, Robb AL (1986).
\dQuote{On Approximating the Statistical Properties of Elasticities.}
\emph{The Review of Economics and Statistics}, \bold{68}, 715--719.

Krinsky I, Robb AL (1990).
\dQuote{On Approximating the Statistical Properties of Elasticities: 
  A Correction.}
\emph{The Review of Economics and Statistics}, \bold{72}, 189--190.

Turnbull BW (1976).
\dQuote{The Empirical Distribution Function with Arbitrarily 
Grouped, Censored and Truncated Data.}
\emph{Journal of the Royal Statistical Society. Series B},
  \bold{38}(3), 290--295.

}

\seealso{
\code{\link{summary.turnbull}}, \code{\link{plot.turnbull}}, \code{\link{kristrom}}, 
\code{\link{sbchoice}}, \code{\link{dbchoice}}, \code{\link[Ecdat]{NaturalPark}}, \code{\link[stats]{glm}}, 
\code{\link[interval]{icfit}}
}

\examples{
## Examples are based on turnbull.db(). turnbull.sb() is similarly 
## implemented. A difference appears in the definition of the model 
## formula. See "Details" section of the help.

## A data set used here is NaturalPark in the package Ecdat (Croissant, 
## 2011): double-bounded dichotomous choice CV style question for 
## measuring willingness to pay for the preservation of the Alentejo 
## Natural Park. The data frame contains seven variables: bid1 (bid in 
## the initial question), bidh (higher bid in the follow-up question), 
## bidl (lower bid in the follow-up question), answers (response 
## outcomes in a factor format with four levels of "nn", "ny", "yn", "yy"), 
## respondents' characteristic variables such as age, sex and income (see 
## NaturalPark for details).
data(NaturalPark, package = "Ecdat")

## The variable answers are converted into a format that is suitable for 
## the function turnbull.db() as follows:
NaturalPark$R1 <- ifelse(substr(NaturalPark$answers, 1, 1) == "y", 1, 0)
NaturalPark$R2 <- ifelse(substr(NaturalPark$answers, 2, 2) == "y", 1, 0)

## The variables bidh and bidl are integrated into one variable (bid2) 
## as follows:
NaturalPark$bid2 <- ifelse(NaturalPark$R1 == 1, NaturalPark$bidh, NaturalPark$bidl)

## The formula for turnbull.sb and turnbull.db are defined respectively as follows:
fmts <- R1 ~ bid1
fmtd <- R1 + R2 ~ bid1 + bid2

## The function turnbull.db() with the formula fmtd and the data frame 
## NaturalPark is executed as follows:
NPts <- turnbull.sb(fmts, data = NaturalPark)
NPts
NPtss <- summary(NPts)
NPtss
plot(NPts)


## The function turnbull.db() with the formula fmtd and the data frame 
## NaturalPark is executed as follows:
NPtd <- turnbull.db(fmtd, data = NaturalPark)
NPtd
NPtds <- summary(NPtd)
NPtds
plot(NPtd)
}

\keyword{DCchoice}
\keyword{double-bounded}
\keyword{single-bounded}
\keyword{nonparametric}
\keyword{survival}
