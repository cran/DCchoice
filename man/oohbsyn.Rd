\name{oohbsyn}
\alias{oohbsyn}
\docType{data}
\title{Synthetic data set for oohbchoice()}
\description{
Dataset created artificially for the examples section of the function \code{oohbhoice()}.
}

\usage{
data(oohbsyn)
}

\format{

A data frame with 80 observations on the following variables.
    \describe{
      \item{id}{a vector of the identification number of the respondent.}
      \item{gender}{a vector containing the gender of the respondent, taking \code{male} or \code{female}.}
      \item{age}{a vector containing the age of the respondent.}
      \item{BL}{a vector of lower bid.}
      \item{BH}{a vector of higher bid.}
      \item{R1}{a vector of response to the first stage CV question, taking a value of \code{1} if the bid is accepted, and \code{0} otherwise.}
      \item{R2}{a vector of response to the second stage CV question, taking a value of \code{1} if the bid is accepted, \code{0} if the bid is not accepted, and \code{-9} if the respondent has no the second stage question.}
    } 
}

\seealso{
\code{\link{oohbchoice}}, \code{\link{turnbull.oohb}}
}

\examples{
## Parametric model
data(oohbsyn)
oohb1 <- oohbchoice(R1 + R2 ~ 1 | log(BL) + log(BH), data = oohbsyn)
summary(oohb1)
oohb2 <- oohbchoice(R1 + R2 ~ age + gender | log(BL) + log(BH), data = oohbsyn)
summary(oohb2)

## Non-parametric model
oohb3 <- turnbull.oohb(R1 + R2 ~ BL + BH, data = oohbsyn)
summary(oohb3)
plot(oohb3)
}


\keyword{DCchoice}
\keyword{datasets}
