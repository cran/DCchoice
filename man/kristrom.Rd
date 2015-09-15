\name{kristrom}
\alias{kristrom}
\alias{print.kristrom}


\title{
The Kristr\"{o}m's nonparametric approach to analyze single-bounded dichotomous choice contingent valuation data
}

\description{
This function analyzes single-bounded dichotomous choice contingent valuation (CV) data on the basis 
of the Kristr\"{o}m's nonparametric method.
}

\usage{
kristrom(formula, data)

\method{print}{kristrom}(x, digits = 4, \dots)
}

\arguments{
    \item{formula}{a formula specifying the model structure.}
    \item{data}{a data frame containing the variables in the model formula.}
    \item{x}{an object of class \code{"kristrom"}.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\details{
The function \code{kristrom()} analyzes single-bounded dichotomous choice contingent 
valuation (CV) data on the basis of Kristr\"{o}m's nonparametric method (Kristr\"{o}m 1990).

The argument \code{formula} defines the response variables and bid variables. 
The argument \code{data} is set as a data frame containing the variables in the model. 

A typical structure of the formula for \code{kristrom()} is defined as follows:

\code{R1 ~ BD1}

The formula consists of two parts. The first part, the left-hand side of the tilde sign 
(\code{~}), must contain the response variable (e.g., \code{R1}) for the suggested prices 
in the CV questions. The response variable contains \code{"Yes"} or \code{"No"} to the bid 
or \code{1} for \code{"Yes"} and \code{0} for \code{"No"}. The other part, which starts 
after the tilde sign, must contain a bid variable (e.g., \code{BD1}) containing suggested 
prices in the CV question.

The structure of data set which assigned to the argument data is the same as that in case 
of \code{dbchoice()}. See \code{\link{dbchoice}} for details in the data set structure.
}

\value{
The function \code{kristrom()} returns an object of S3 class \code{"kristrom"}. 
An object of \code{"kristrom"} is a list with the following components. 

\item{tab.dat}{a matrix describing the number of respondents who answered \code{"yes"} to CV question, 
    total number of respondents, and ratio of respondents who answered \code{"yes"} among the total number of 
    respondents for each value of the suggested bids.}
\item{M}{the number of rows of \code{tab.dat}.}
\item{adj.p}{a vector describing the probability of a yes-answer to the suggested bid, which is the same 
    as the last column of \code{tab.dat}.}
\item{nobs}{the number of observations.}
\item{unq.bid}{a vector of the unique bids.}
\item{estimates}{a matrix of the estimated Kristr\"{o}m's survival probabilities.}

The generic function \code{print()} is available for fitted model object of class \code{"kristrom"} 
and displays the estimated Kristr\"{o}m's survival probabilities. 

The extractor function \code{summary()} is used to display the estimated Kristr\"{o}m's survival 
probabilities as well as three types of WTP estimates (Kaplan-Meier and Spearman-Karber mean, and median estimates). Note that the Spearman-Karber mean estimate is computed 
upto the X-intercept.

A graph of the estimated empirical survival function is depicted by \code{plot()}. 
See \code{\link{plot.kristrom}} for details. 

\code{\link{turnbull.sb}} is an alternative nonparametric method for analyzing single-bounded 
dichotomous choice data. A parametric analysis can be done by \code{\link{sbchoice}}.

}


\references{
Croissant Y (2011).
\emph{\pkg{Ecdat}: Data Sets for Econometrics,} 
\R package version 0.1-6.1, 
\url{http://CRAN.R-project.org/package=Ecdat}.

Kristr\"{o}m B (1990). 
\dQuote{A Non-Parametric Approach to the Estimation of Welfare 
Measures in Discrete Response Valuation Studies.}
\emph{Land Economics}, \bold{66}(2), 135--139.
}

\seealso{
    \code{\link{plot.kristrom}}, \code{\link[Ecdat]{NaturalPark}}, 
    \code{\link{turnbull.sb}}, \code{\link{sbchoice}} 
}

\examples{
## Examples for kristrom() are also based on a data set NaturalPark in the package 
## Ecdat (Croissant 2011): so see the section Examples in the dbchoice() for details.
data(NaturalPark, package = "Ecdat")

## The variable answers are converted into a format that is suitable for the function 
## kristrom() as follows:
NaturalPark$R1 <- ifelse(substr(NaturalPark$answers, 1, 1) == "y", 1, 0)

## The formula is defined as follows:
fmks <- R1 ~ bid1

## The function kristrom() with the function fmks and the data frame NP 
## is executed as follows:
NPks <- kristrom(fmks, data = NaturalPark)
NPks
NPkss <- summary(NPks)
NPkss
plot(NPks)
}

\keyword{DCchoice}
\keyword{single-bounded}
\keyword{nonparametric}
\keyword{survival}
