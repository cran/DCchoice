\name{dbchoice}
\alias{dbchoice}
\alias{print.dbchoice}
\alias{vcov.dbchoice}
\alias{logLik.dbchoice}

\title{Parametric approach to analyze double-bounded dichotomous choice contingent valuation data}
\description{This function analyzes double-bounded dichotomous choice contingent 
            valuation (CV) data on the basis of the utility difference approach.}
\usage{
dbchoice(formula, data, subset, na.action = na.omit, dist = "log-logistic",
         par = NULL, \dots)

\method{print}{dbchoice}(x, digits = max(3, getOption("digits") - 1), \dots)

\method{vcov}{dbchoice}(object, \dots)

\method{logLik}{dbchoice}(object, \dots)

}

\arguments{
    \item{formula}{an object of S3 class \code{"formula"} and specifies the model structure.}
    \item{data}{a data frame containing the variables in the model formula.}
    \item{subset}{an optional vector specifying a subset of observations.}
    \item{na.action}{a function which indicates what should happen when the data contains \code{NA}s.}
    \item{dist}{a character string setting the error distribution in the model, which 
                takes one of \code{"logistic"}, \code{"normal"}, \code{"log-logistic"}, 
                \code{"log-normal"} or \code{"weibull"}.}
    \item{par}{a vector of initial parameters over which the optimization is carried out.}
    \item{x}{an object of class \code{"dbchoice"}.}
    \item{digits}{a number of digits to display.}
    \item{object}{an object of class \code{"dbchoice"}.}
    \item{...}{optional arguments. Currently not in use.}
}

\details{

The function \code{dbchoice()} implements an analysis of double-bounded dichotomous choice 
contingent valuation (CV) data on the basis of the utility difference approach (Hanemann, 1984). 
A generic call to \code{dbchoice()} is given by

\code{dbchoice(formula, data, dist = "log-logistic", ...)}

The extractor function \code{summary()} is available for a \code{"dbchoice"} class object. 
See \code{\link{summary.dbchoice}} for details.

There are two functions available for computing the confidence intervals for the estimates of WTPs. 
\code{\link{krCI}} implements simulations to construct empirical distributions of the WTP while 
\code{\link{bootCI}} carries out nonparametric bootstrapping. 

The argument \code{formula} defines the response variables and covariates. The argument \code{data} 
is mandatory where the data frame containing the variables in the model is specified. The argument \code{dist} 
sets the error distribution. 
Currently, one of \code{"logistic"}, \code{"normal"}, \code{"log-logistic"}, \code{"log-normal"}, 
or \code{"weibull"} is available. 
The default value is \code{dist = "log-logistic"}, so that it may be omitted if the user wants to estimate a model 
with log-logistic error distribution.

The difference between normal and log-normal models or between logistic or 
log-logistic ones is how the bid variable is incorporated into the model to be estimated. 
For the Weibull model, the bid variable must be entered in the natural log. 
Therefore, the user must be careful in defining the model formula that is explained in details 
below.

A typical structure of the formula for \code{dbchoice()} is defined as follows:

\code{R1 + R2 ~ (the names of the covariates) | BD1 + BD2}

The formula is an object of class \code{"formula"} and specifies the model structure. It has to be written in a 
symbolic expression in \R.
The formula consists of three parts. The first part, the left-hand side of the tilde sign (\code{~}), must 
contain the response variables for the suggested prices in the first and the second stage of CV questions. 
In the example below, \code{R1} denotes a binary or two-level factor response variable for a bid in the first 
stage and \code{R2} for a bid in the second stage. Each of \code{R1} and \code{R2} contains \code{"Yes"} or 
\code{"No"} to the bid or \code{1} for \code{"Yes"} and \code{0} for \code{"No"}.

The covariates are defined in the second part in the place of \code{(the names of the covariates)}. Each covariate is 
connected with the arithmetic operator \code{+} and \code{(the names of the covariates)} in the above syntax should be 
replaced with \code{var1 + var2} and the like. The plus sign is nothing to do with addition of the two variables 
in the symbolic expression. When the covariate contains only a constant term, a value of \code{1} is set as the 
covariate (that is, \code{R1 + R2 ~ 1 | BD1 + BD2})

The last part starts after the vertical bar (\code{|}). The names of the two variables (\code{BD1} and \code{BD2}) 
containing suggested prices in the first and second stage of double-bounded dichotomous choice CV question are 
specified in this part. The two variables are also connected with the arithmetic operator (\code{+}).

According to the structure of the formula, a data set (data frame) consists of three parts. An example of the 
data set is as follows (\code{sex}, \code{age}, and \code{income} are respondents characteristics and assumed 
to be covariates):

\tabular{rrrrrrr}{
    \code{R1}  \tab \code{R2}   \tab \code{sex}  \tab \code{age}  \tab \code{income} \tab \code{BD1} \tab \code{BD2}\cr
    Yes \tab Yes  \tab Male \tab 20   \tab Low    \tab 100 \tab 250\cr
    Yes \tab No   \tab Male \tab 30   \tab Low    \tab 500 \tab 1000\cr
    \dots \tab \tab \tab \tab \tab \tab
}

The second bid in the double-bounded dichotomous choice CV question is larger or lower than the first bit 
according to the response to the first stage: if the response to the first stage is \code{"Yes"}, the second bid is 
larger than the first bid; if the response is \code{"No"}, the second bid is lower than the first bid. In the example 
above, \code{BD2} is set as the second bid according to each respondent faced in the second stage. However, 
the followings style of data set is frequently prepared:

\tabular{rrrrrrrr}{
    \code{R1}  \tab \code{R2}  \tab \code{sex}  \tab \code{age} \tab \code{income} \tab \code{BD1} \tab \code{BD2H} \tab \code{BD2L}\cr
    Yes \tab Yes \tab Male \tab 20  \tab Low    \tab 100 \tab 250  \tab 50\cr
    Yes \tab No  \tab Male \tab 30  \tab Low    \tab 500 \tab 1000 \tab 250\cr
    \dots \tab \tab \tab \tab \tab \tab \tab 
}
\code{BD2H} is the second (higher) bid when the respondent answers \code{"Yes"} in the first stage; \code{BD2L} 
is the second (lower) bid when the respondent answers \code{"No"} in the first stage. In this case, the users 
have to convert \code{BD2H} and \code{BD2L} into \code{BD2} (see the section "Examples").

The function \code{dbchoice()} analyzes double-bounded dichotomous choice CV data using the function 
\code{\link[stats]{optim}} on the basis of the initial coefficients that are estimated from a binary 
logit model analysis of the first-stage CV responses (the binary logit model is estimated internally 
by the function \code{\link[stats]{glm}} with the argument \code{family = binomial(link = "logit"))}.

Nonparametric analysis of double-bounded dichotomous choice data can be done by \code{\link{turnbull.db}}. 
A single-bounded analogue of \code{dbchoice} is called \code{\link{sbchoice}}.

}

\value{

This function returns an S3 class object \code{"dbchoice"} that is a list with the following components.

\item{f.stage}{a list of components returned from the function \code{\link[stats]{glm}} based on the 
    responses to the first CV question. The coefficient estimates of the first stage estimation is 
    used as the initial coefficients for full analysis using the function \code{\link[stats]{optim}}.
    If \code{par} is not \code{NULL}, the supplied vector is returned.}
\item{dbchoice}{a list of components returned from the function \code{\link[stats]{optim}}.}
\item{coefficients}{a named vector of estimated coefficients.}
\item{call}{the matched call.}
\item{formula}{the formula supplied.}
\item{Hessian}{an estimate of the Hessian. See also \code{Hessian} in \code{\link[stats]{optim}}.}
\item{distribution}{a character string showing the error distribution used.}
\item{loglik}{a value of the log likelihood at the estimates.}
\item{convergence}{an logical code: \code{TRUE} means a successful convergence.}
\item{niter}{a vector of two integers describing the number of calls to the object function and the 
    numerical gradient, respectively. See also \code{counts} in \code{\link[stats]{optim}}.}
\item{nobs}{a number of observations.}
\item{covariates}{a named matrix of the covariates used in the model.}
\item{bid}{a named matrix of the bids used in the model.}
\item{yn}{a named matrix of the responses to the initial and follow-up CV questions used in the model.}
\item{data.name}{the data matrix.}
\item{terms}{terms}
\item{contrast}{contrasts used for factors}
\item{xlevels}{levels used for factors}
}

\references{
Bateman IJ, Carson RT, Day B, Hanemann M, Hanley N, Hett T, Jones-Lee M, Loomes
  G, Mourato S, \"{O}zdemiro\={g}lu E, Pearce DW, Sugden R, Swanson J (eds.) (2002).
\emph{Economic Valuation with Stated Preference Techniques: A Manual.} 
Edward Elger, Cheltenham, UK.

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

Hanemann WM, Loomis JB, Kanninen BJ (1991).
\dQuote{Statistical Efficiency of Double-Bounded Dichotomous Choice 
Contingent Valuation.} 
\emph{American Journal of Agricultural Economics}, \bold{73}(4), 1255--1263.

}

\seealso{
\code{\link{summary.dbchoice}}, \code{\link{krCI}}, \code{\link{bootCI}}, 
\code{\link{sbchoice}}, \code{\link{turnbull.db}}, \code{\link[Ecdat]{NaturalPark}}, 
\code{\link[stats]{glm}}, \code{\link[stats]{optim}}, \code{\link[stats]{formula}}
}

\examples{
## Examples are based on a data set NaturalPark in the package 
## Ecdat (Croissant 2011): DBDCCV style question for measuring 
## willingness to pay for the preservation of the Alentejo Natural 
## Park. The data set (dataframe) contains seven variables: 
## bid1 (bid in the initial question), bidh (higher bid in the follow-up 
## question), bidl (lower bid in the follow-up question), answers 
## (response outcomes in a factor format with four levels of "nn", 
## "ny", "yn", "yy"), respondents' characteristic variables such 
## as age, sex and income (see NaturalPark for details).
data(NaturalPark, package = "Ecdat")
head(NaturalPark)

## The variable answers are converted into a format that is suitable for the 
## function dbchoice() as follows:
NaturalPark$R1 <- ifelse(substr(NaturalPark$answers, 1, 1) == "y", 1, 0)
NaturalPark$R2 <- ifelse(substr(NaturalPark$answers, 2, 2) == "y", 1, 0)

## We assume that the error distribution in the model is a 
## log-logistic; therefore, the bid variables bid1 is converted 
## into LBD1 as follows:
NaturalPark$LBD1 <- log(NaturalPark$bid1)

## Further, the variables bidh and bidl are integrated into one 
## variable (bid2) and the variable is converted into LBD2 as follows:
NaturalPark$bid2 <- ifelse(NaturalPark$R1 == 1, NaturalPark$bidh, NaturalPark$bidl)
NaturalPark$LBD2 <- log(NaturalPark$bid2)

## The utility difference function is assumed to contain covariates (sex, age, and 
## income) as well as two bid variables (LBD1 and LBD2) as follows:
fmdb <- R1 + R2 ~ sex + age + income | LBD1 + LBD2

\dontrun{
## The formula may be alternatively defined as
fmdb <- R1 + R2 ~ sex + age + income | log(bid1) + log(bid2)
}

## The function dbchoice() with the function fmdb and the dataframe 
## NP is executed as follows:
NPdb <- dbchoice(fmdb, data = NaturalPark)
NPdb
NPdbs <- summary(NPdb)
NPdbs

## The confidence intervals for these WTPs are calculated using the 
## function krCI() or bootCI() as follows:
\dontrun{
krCI(NPdb)
bootCI(NPdb)
}
## The WTP of a female with age = 5 and income = 3 is calculated
## using function krCI() or bootCI() as follows:
\dontrun{
krCI(NPdb, individual = data.frame(sex = "female", age = 5, income = 3))
bootCI(NPdb, individual = data.frame(sex = "female", age = 5, income = 3))
}

## The variable age and income are deleted from the fitted model, 
## and the updated model is fitted as follows:
update(NPdb, .~. - age - income |.)

## The bid design used in this example is created as follows:
bid.design <- unique(NaturalPark[, c(1:3)])
bid.design <- log(bid.design)
colnames(bid.design) <- c("LBD1", "LBDH", "LBDL")
bid.design
## Respondents' utility and probability of choosing Yes-Yes, Yes-No, 
## No-Yes, and No-No under the fitted model and original data are 
## predicted as follows: 
head(predict(NPdb, type = "utility", bid = bid.design))
head(predict(NPdb, type = "probability", bid = bid.design))
## Utility and probability of choosing Yes for a female with age = 5 
## and income = 3 under bid = 10 are predicted as follows:
predict(NPdb, type = "utility",
    newdata = data.frame(sex = "female", age = 5, income = 3, LBD1 = log(10)))
predict(NPdb, type = "probability",
    newdata = data.frame(sex = "female", age = 5, income = 3, LBD1 = log(10)))

## Plot of probabilities of choosing yes is drawn as drawn as follows:
plot(NPdb)
## The range of bid can be limited (e.g., [log(10), log(20)]):
plot(NPdb, bid = c(log(10), log(20)))
}




\keyword{DCchoice}
\keyword{double-bounded}
\keyword{nonlinear}

