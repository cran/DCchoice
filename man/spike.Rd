\name{spike-models}
\alias{spike-models}
\alias{spike}
\alias{spike models}
\alias{sbspike}
\alias{oohbspike}
\alias{dbspike}
\alias{spikeCoef}
\alias{print.spike}
\alias{summary.spike}
\alias{print.summary.spike}
\alias{vcov.spike}
\alias{logLik.spike}
\alias{plot.spike}


\title{Parametric approach to analyze dichotomous choice contingent valuation data on the basis of a simple spike model}

\description{These functions implement a simple spike model analysis of single-, one-and-one-half-, and double-bounded dichotomous choice contingent valuation data using the maximum likelihood method.}

\usage{
## for the single-bounded data
sbspike(formula, data, subset, na.action = na.omit, par = NULL, \dots)

## for the one-and-one-half-bounded data
oohbspike(formula, data, subset, na.action = na.omit, par = NULL, \dots)

## for the double-bounded data
dbspike(formula, data, subset, na.action = na.omit, par = NULL, \dots)

\method{print}{spike}(x, digits = max(3, getOption("digits") - 1), \dots)

\method{summary}{spike}(object, \dots)

\method{print}{summary.spike}(x, digits = max(3, getOption("digits") - 1), \dots)

\method{logLik}{spike}(object, \dots)

\method{vcov}{spike}(object, \dots)

\method{plot}{spike}(x, type = "l", main = NULL, sub = NULL, 
  xlab = "Bid", ylab = "Probability", lwd = 3, lty = 1, 
  xlim = c(0, max(x$bid)), ylim = c(0, 1), bid = c(0, max(x$bid)), \dots)
}

\arguments{
    \item{formula}{an object of S3 class \code{'Formula'} specifying the model structure.}
    \item{data}{a data frame containing the variables in the model formula.}
    \item{subset}{an optional vector specifying a subset of observations.}
    \item{na.action}{a function which indicates what should happen when the data contains \code{NA}s.}
    \item{par}{a vector of initial parameters over which the optimization is carried out.}
    \item{x, object}{an object of class \code{'spike'}.}
    \item{digits}{the number of digits to display.}
    \item{type}{type of plot.}
    \item{main}{the main title of the plot. If unspecified, no main title is displayed.}
    \item{sub}{the sub-title of the plot. If unspecified, no sub-title is displayed.}
    \item{xlab}{the x label of the plot. The default is \code{"Bid"}. }
    \item{ylab}{the y label of the plot. The default is \code{"Probability"}.}
    \item{lwd}{the line width for the plot. The default is \code{3}.}
    \item{lty}{the line type for the plot. The default is \code{1}.}
    \item{xlim}{the x limits of the plot. The default is \code{c(0, max(x$bid))}.}
    \item{ylim}{the y limits of the plot. The default is \code{c(0, 1)}.}
    \item{bid}{the bid limits that should be drawn. The default is \code{c(0, max(x$bid))}.}
    \item{\dots}{optional arguments.}
}

\details{

The functions \code{\link{sbspike}}, \code{\link{oohbspike}}, and \code{\link{dbspike}} implement a spike model analysis of single-, one-and-one-half-, and double-bounded dichotomous choice contingent valuation (SB, OOHB, and DB DCCV) data, respectively. A simple spike model assumes a non-zero probability of zero willingness to pay (WTP) for a good/service and a zero probability of negative WTP. These functions are developed according to the original simplest spike model proposed by Kristr\"{o}m (1997) and its follow-up studies (i.e., Yoo and Kwak (2002) for DB DCCV and Kwak et al. (2013) for OOHB DCCV). These functions use a maximum likelihood methods to fit the models with the CV data.

Since the usage of spike model functions \code{\link{sbspike}}, \code{\link{oohbspike}}, and \code{\link{dbspike}} are similar to the un-spike (the ordinary) DCCV model functions \code{\link{sbchoice}}, \code{\link{oohbchoice}}, and \code{\link{dbchoice}}, respectively, this help below explains only the differences in usage between the spike and ordinary model functions. We assume that you understand how to use the ordinary model functions \code{\link{sbchoice}}, \code{\link{oohbchoice}}, and \code{\link{dbchoice}}. If you are unfamiliar with the ordinary model functions, please refer to helps for these at first.

The first difference between the spike and ordinal model functions is that an argument \code{distribution} used in the ordinary model functions is not defined in the spike functions: the functions for spike models assume that the error distribution is only a logistic with a spike, and thus the other distributions (i.e., log-logistic, normal, log-normal, and Weibull) are not available to the spike model functions.

The other difference is about an argument \code{formula}, which is assigned an object of the S3 class \code{'\link[Formula]{Formula}'}. For a model formula for the ordinary model functions, the left-handed side of the tilde (\code{~}) contains only response variable(s) (i.e., the response to SB DCCV question, \code{R1}, for \code{\link{sbchoice}}; the response to the first stage of OOHB/DB DCCV question, \code{R1}, and the second one, \code{R2}, for \code{\link{oohbchoice}} and \code{\link{dbchoice}}), while it contains both the response variable(s) and spike variable for the spike model functions. The spike variable, \code{S}, which must be set in the second part (after the vertical bar [\code{|}]) of the left-handed side of the tilde, takes the value of \code{1} if the respondent has a positive WTP for a good specified in the DCCV question and \code{0} otherwise. See Kristr\"{o}m (1997) for a question to measure whether the respondent has a positive WTP or not. A typical structure of the formula for spike model functions consists of the following four parts:

for \code{sbspike()},
\code{R1 | S ~ <the names of the covariates> | BD1}

for \code{dbspike()},
\code{R1 + R2 | S ~ <the names of the covariates> | BD1 + BD2}

and for \code{oohbspike()},
\code{R1 + R2 | S ~ <the names of the covariates> | BL + BH}

where \code{BD1} and \code{BD2} are variables containing suggested prices in the first and second stages of the SB/DB DCCV question; and \code{BL} and \code{BH} are variables containing suggested lower and higher prices in the OOHB DCCV question.

According to the structure of the formula, a data set (data frame) consists of four parts. An example of the data set for \code{\link{dbspike}} is as follows (\code{sex}, \code{age}, and \code{income} are respondent characteristics and assumed to be covariates):

\tabular{rrrrrrrr}{
  \code{R1} \tab \code{R2} \tab \code{S} \tab \code{sex} \tab \code{age} \tab \code{income} \tab \code{BD1} \tab \code{BD2}\cr
  \code{Yes} \tab \code{Yes} \tab \code{1} \tab \code{Male} \tab \code{20} \tab \code{Low} \tab \code{100} \tab \code{250}\cr
  \code{Yes} \tab \code{No} \tab \code{0} \tab \code{Male} \tab \code{30} \tab \code{Low} \tab \code{500} \tab \code{1000}\cr
  \dots \tab \tab \tab \tab \tab \tab \tab
}

The spike model functions fit the models with DCCV data using the function \code{\link[stats]{optim}} on the basis of the initial coefficients estimated from an un-spike (ordinary) binary logit model analysis of the response to the SB DCCV question, or the first-stage response to the OOHB/DB DCCV question. The binary logit model is estimated internally using the function \code{\link[stats]{glm}} with the argument \code{family = binomial(link = "logit")}.

The spike model functions return an S3 \code{'spike'} class object. Various methods for the S3 \code{"spike"} class object are provided as follows: \code{print()} displays estimated coefficients; \code{summary()} extracts detailed information on the fitted model; \code{summary.print()} displays information extracted by \code{summary()}; \code{logLik()} extracts the value of a log-likelihood function at estimates; \code{vcov()} returns the variance-covariance matrix of the fitted model; and \code{plot()} draws an estimated survival distribution of the WTP according to the fitted model. These S3 methods correspond to those for the ordinary DCCV functions \code{\link{sbchoice}}, \code{\link{oohbchoice}}, and \code{\link{dbchoice}}. Therefore, for details, see helps for the corresponding methods for ordinary DCCV functions. Note that the mean and median WTPs calculated by \code{summary()} for the spike model functions are defined as follows (see Kristr\"{o}m 1997): mean WTP = ln(1 + exp(A))/B if the parameter for a bid variable (B) is positive (A is the constant), and NA otherwise; median WTP = A/B if 1/(1 + exp(-A)) < 0.5, and 0 otherwise. When covariates are included in the fitted model, the constant in the mean and median WTPs is replaced with \bold{x'b}, where \bold{x} is a row vector of covariates at the sample mean including the value of 1 for the constant, and \bold{b} is a column vector of estimates for covariates including the constant. See Yoo and Kwak (2009), Kwak et al. (2013), and Lee et al. (2010) for SB, OOHB, and DB spike models with covariates, respectively. 

The existing functions \code{\link{bootCI}} and \code{\link{krCI}}, which compute the confidence intervals for the estimates of WTPs using non-parametric and parametric bootstrapping approaches respectively, were revised to handle an S3 \code{'spike'} class object. An existing function \code{\link{ct2df}} was also updated to handle a data set in contingency-table format for spike model functions.

Furthermore, a new function \code{\link{spikeCoef}} was developed to estimate a spike for the fitted model as 1/(1 + exp(A)), where A is the constant. This function returns the estimated spike, its standard error, and the corresponding z- and p-values under the null hypothesis where the spike is zero. The standard error is calculated using the delta method. When covariates are included in the fitted model, the constant in the formula is replaced with \bold{x'b} as the mean and median WTP calculations. See the examples below, for details.

}

\value{

These spike model functions return an S3 class object \code{'spike'}, which is a list with the following components.

\item{f.stage}{a list of components returned from the un-spike (ordinary) binary logit mode analysis using the function \code{\link[stats]{glm}} based on the responses to the SB DCCV question, or the first stage of the OOHB/DB DCCV question. If the argument \code{par} is not \code{NULL}, the supplied vector is returned.}

\item{optim.out}{a list of components returned from the function \code{\link[stats]{optim}}.}

\item{coefficients}{a named vector of estimated coefficients.}

\item{call}{the matched call.}

\item{formula}{the formula supplied.}

\item{Hessian}{an estimate of the Hessian. See also Hessian in \code{\link[stats]{optim}}.}

\item{loglik}{a value of the log likelihood at the estimates.}

\item{convergence}{a logical code: \code{TRUE} means a successful convergence.}

\item{niter}{a vector of two integers describing the number of calls to the object function and numerical gradient, respectively. See also counts in \code{\link[stats]{optim}}.}

\item{nobs}{a number of observations.}

\item{covariates}{a named matrix of the covariates used in the model.}

\item{bid}{a named matrix of the bids used in the model.}

\item{yn}{a named matrix of the responses to the SB DCCV question or the first and second stage of the OOHB/DB DCCV question used in the model.}

\item{data.name}{the data matrix.}

\item{terms}{terms.}

\item{contrast}{contrasts used for factors.}

\item{xlevels}{levels used for factors.}

}

\references{

Kristr\"{o}m B. (1997) Spike models in contingent valuation. \emph{American Journal of Agricultural Economics} \bold{79}: 1013--1023.

Yoo S-H, Kwak S-J. (2002) Using a spike model to deal with zero response data from double bounded dichotomous choice contingent valuation surveys. \emph{Applied Economics Letters} \bold{9}: 929--932.

Kwak S-J, Yoo S-H, Kim C-S. (2013) Measuring the willingness to pay for tap water quality improvements: results of a contingent valuation survey in Pusan. \emph{Water} \bold{5}: 1638--1652.

Lee J-S, Yoo S-H, Kwak S-J. (2010) Public's willingness to pay for preventing climate change. \emph{Applied Economics Letters} \bold{17}: 619--622.

Yoo S-H, Kwak S-Y. (2009) Willingness to pay for green electricity in Korea: a contingent valuation study. \emph{Energy Policy} \bold{37}: 5408--5416.

}

\seealso{
\code{\link{sbchoice}}, \code{\link{oohbchoice}}, \code{\link{dbchoice}},
\code{\link{ct2df}}, \code{\link{krCI}}, \code{\link{bootCI}}, 
\code{\link{CarsonSB}}, \code{\link{CarsonDB}}, \code{\link{oohbsyn}},
\code{\link[stats]{glm}}, \code{\link[stats]{optim}}, 
\code{\link[Formula]{Formula}}
}

\examples{
# Example datasets were created by modifying CarsonSB, CarsonDB, and oohbsyn.

# Spike SB Example
sb <- data.frame(
 bid = c(10, 30, 60, 120),
 y   = c(178, 138, 129, 88),
 ny  = c(56, 45, 50, 76),
 nn  = c(30, 84, 76, 93))
SB <- ct2df(sb, bid1 = "bid", type = "single", spike = TRUE)
head(SB)
dim(SB)
SBout <- sbspike(R1 | S ~ 1 | bid1, data = SB)
summary(SBout)
spikeCoef(SBout)
\dontrun{
krCI(SBout)
bootCI(SBout)}
plot(SBout, main = "Spike SB model")

# Spike DB Example
db <- data.frame(
 bidh = c(30, 60, 120),
 bid1 = c(10, 30, 60),
 bidl = c(5, 10, 30),
 yy   = c(119, 69, 54),
 yn   = c(59, 69, 75),
 ny   = c(8, 31, 25),
 nny  = c(47, 61, 35),
 nnn  = c(31, 37, 66))
DB <- ct2df(x = db, type = "double", spike = TRUE)
head(DB)
dim(DB)
DBout <- dbspike(R1 + R2 | S ~ 1 | bid1 + bid2, data = DB)
summary(DBout)
spikeCoef(DBout)
\dontrun{
krCI(DBout)
bootCI(DBout)}
plot(DBout, main = "Spike DB model")

# Spike OOHB Example
oohb <- data.frame(
 bidl = c(2, 4, 6, 8),
 bidh = c(4, 6, 8, 10),
 yy   = c(8, 6, 4, 2),
 yn   = c(1, 3, 1, 1),
 n_y  = c(1, 1, 4, 4),
 n_n  = c(0, 1, 1, 3),
 y    = c(7, 6, 3, 1),
 ny   = c(2, 2, 3, 1),
 nn_y = c(1, 1, 2, 5),
 nn_n = c(0, 0, 2, 3))
OOHB <- ct2df(x = oohb, type = "oohb", spike = TRUE)
head(OOHB)
dim(OOHB)
OOHBout <- oohbspike(R1 + R2 | S ~ 1 | BL + BH, data = OOHB)
summary(OOHBout)
spikeCoef(OOHBout)
\dontrun{
krCI(OOHBout)
bootCI(OOHBout)}
plot(OOHBout, main = "Spike OOHB model")
}

\keyword{DCchoice}
\keyword{spike}
