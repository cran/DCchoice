\name{sbchoice}
\alias{sbchoice}
\alias{print.sbchoice}

\title{Parametric approach to analyze single-bounded dichotomous choice contingent valuation data}
\description{
This function analyzes single-bounded dichotomous choice contingent valuation (CV) data 
on the basis of the utility difference approach.
}

\usage{
sbchoice(formula, data, dist = "log-logistic", \dots)

\method{print}{sbchoice}(x, digits = max(3, getOption("digits") - 1), \dots)
}

\arguments{
    \item{formula}{an S3 class object \code{"formula"} and specifies the model structure.}
    \item{data}{a data frame containing the variables in the model formula.}
    \item{dist}{a character string setting the error distribution in the model, which takes one 
        of \code{"logistic"}, \code{"normal"}, \code{"log-logistic"}, \code{"log-normal"} or 
        \code{"weibull"}.}
    \item{x}{an object of class \code{"sbchoice"}.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\details{
The function \code{sbchoice()} implements an analysis of single-bounded dichotomous choice 
contingent valuation (CV) data on the basis of the utility difference approach 
(Hanemann, 1984).

The extractor function \code{summary()} is available for a \code{"sbchoice"} class object. 
See\cr
\code{\link{summary.sbchoice}} for details.

There are two functions available for computing the confidence intervals for the estimates of WTPs. 
\code{\link{krCI}} implements simulations to construct empirical distributions of the WTP while 
\code{\link{bootCI}} carries out nonparametric bootstrapping. 

Most of the details of \code{sbchoice()} is the same as those of \code{dbchoice()}, a double-bounded 
analogue of \code{sbchoice}. See the section \strong{Details} in \code{\link{dbchoice}}. Differences between 
the two functions are as follows:
\enumerate{
    \item In the model formula, the first part contains only one response variable (e.g., \code{R1}) 
    and the third part contains only one bid variable (e.g., \code{BD1}) because respondents are 
    requested to answer a CV question in the single-bounded dichotomous choice CV. The following 
    is a typical structure of the formula:\cr
      \code{R1 ~ (the names of the covariates) | BD1}
    \item The function \code{sbchoice()} analyzes the responses to single-bounded dichotomous choice CV 
    questions internally using the function \code{glm()} with the argument \cr
    \code{family = binomial(link = "logit")}\cr 
    or \code{family = binomial(link = "probit")}.\cr
    When \code{dist = "weibull"}, optimization is carries out using the \code{optim()} function 
    with a hard-coded log-likelihood function.
    \item Outputs from \code{sbchoice()} are slightly different from those from \code{dbchoice()} 
    because the analysis in \code{sbchoice()} internally depends on the function \code{glm()} 
    for the (log-) normal or (log-) logistic distributions. (see the \bold{Value} section).
}

Nonparametric analysis of single-bounded dichotomous choice data can be done by \code{\link{turnbull.sb}} 
or by \code{\link{kristrom}}. 

}

\value{
This function returns an object of S3 class \code{"sbchoice"} that is a list with the following components.
\item{coefficients}{a named vector of estimated coefficients.}
\item{call}{the matched call.}
\item{formula}{the formula supplied.}
\item{glm.out}{a list of components returned from \code{glm()} with the data set and the formula. In case 
      of the Weibull distribution, a list of components from the \code{optim()} is saved.}
\item{glm.null}{a list of components returned from \code{glm()} with the data set and a formula 
    containing only constant (null model). In case of the Weibull distribution, a list of 
    components from the \code{optim()} is saved.}
\item{distribution}{a character string showing the error distribution used.}
\item{nobs}{a number of observations.}
\item{covariates}{a named matrix of the covariates used in the model.}
\item{bid}{a named matrix of the bids used in the model.}
\item{yn}{a named matrix of the responses to the CV question used in the model.}
\item{data.name}{the data matrix.}
}


\references{
Bateman IJ, Carson RT, Day B, Hanemann M, Hanley N, Hett T, Jones-Lee M, Loomes
  G, Mourato S, \"{O}zdemiro\={g}lu E, Pearce DW, Sugden R, Swanson J (eds.) (2002).
\emph{Economic Valuation with Stated Preference Techniques: A Manual.} 
Edward Elger, Cheltenham, UK.

Boyle KJ, Welsh MP, Bishop RC (1988).
\dQuote{Validation of Empirical Measures of Welfare Change: Comment.}
\emph{Land Economics}, \bold{64}(1), 94--98.

Carson RT, Hanemann WM (2005).
\dQuote{Contingent Valuation.}
in KG M\"{a}ler, JR Vincent (eds.), \emph{Handbook of Environmental Economics}. 
Elsevier, New York.

Croissant Y (2011).
\emph{\pkg{Ecdat}: Data Sets for Econometrics,} 
\R package version 0.1-6.1, 
\url{http://CRAN.R-project.org/package=Ecdat}.

Hanemann, WM (1984). 
\dQuote{Welfare Evaluations in Contingent Valuation Experiments with Discrete Responses},
\emph{American Journal of Agricultural Economics},
\bold{66}(2), 332--341.

Hanemann M, Kanninen B (1999).
\dQuote{The Statistical Analysis of Discrete-Response {CV} Data.}, 
in IJ Bateman, KG Willis (eds.), 
\emph{Valuing Environmental Preferences: Theory and Practice of the Contingent 
Valuation Methods in the {US}, {EU}, and Developing Countries}, 
302--441. 
Oxford University Press, New York.

}

\seealso{
\code{\link{summary.sbchoice}}, 
\code{\link{krCI}}, \code{\link{bootCI}}, 
\code{\link[Ecdat]{NaturalPark}}, 
\code{\link{turnbull.sb}}, \code{\link{kristrom}}
\code{\link[stats]{glm}}, \code{\link[stats]{formula}}
\code{\link{dbchoice}}
}


\examples{
## Examples for sbchoice() are also based on a data set NaturalPark 
## in the package Ecdat (Croissant 2011): so see the section Examples 
## in the dbchoice() for details.
data(NaturalPark, package = "Ecdat")

## The variable answers are converted into a format that is suitable for 
## the function sbchoice() as follows:
NaturalPark$R1 <- ifelse(substr(NaturalPark$answers, 1, 1) == "y", 1, 0)
NaturalPark$R2 <- ifelse(substr(NaturalPark$answers, 2, 2) == "y", 1, 0)


## We assume that the error distribution in the model is a log-logistic; 
## therefore, the bid variables bid1 is converted into LBD1 as follows:
NaturalPark$LBD1 <- log(NaturalPark$bid1)

## The utility difference function is assumed to contain covariates 
## (sex, age, and income) as well as the bid variable (LBD1) as follows 
## (R2 is not used because of single-bounded dichotomous choice CV format):
fmsb <- R1 ~ sex + age + income | LBD1

\dontrun{
## The formula may be alternatively defined as
fmsb <- R1 ~ sex + age + income | log(bid1)
}

## The function sbchoice() with the function fmsb and the data frame NP 
## is executed as follows:
NPsb <- sbchoice(fmsb, data = NaturalPark)
NPsb
NPsbs <- summary(NPsb)
NPsbs

\dontrun{
## Generic functions such as summary() and coefficients() work for glm.out
summary(NPsb$glm.out)
coefficients(NPsb$glm.out)

## The confidence intervals for these WTPs are calculated using the 
## function krCI() or bootCI() as follows:
krCI(NPsb)
bootCI(NPsb)
}
}

\keyword{DCchoice}
\keyword{nonlinear}
