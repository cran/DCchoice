\name{dccvResponse}
\alias{dccvResponse}

\title{
Generate artificial responses to dichotomous choice CV questions
}

\description{
This function synthesizes responses to DCCV questions.
}

\usage{
dccvResponse(bid, b, n, seed = NULL)
}

\arguments{
  \item{bid}{
a vector (for SBDC) or matrix (for DBDC) containing a bid design
}
  \item{b}{
a vector containing parameters of independent variables in the model
}
  \item{n}{
an integer value showing the number of respondents in the resultant dataset
}
  \item{seed}{
seed for a random number generator
}
}

\details{
This function synthesizes response to SBDC- or DBDC-CV questions on the basis of the utility difference approach. See Aizaki et al. (2014) for detail on the model.

See the examples below for usage in detail.
}

\value{
This function returns a data frame, in which each row shows a single respondent. It contains the following variables.

For SBDC-CV,
  \item{R1}{
a response to a bid: 1 for "Yes", 0 for "No"
}
  \item{bid1}{
the bid
}

For DBDC-CV,
  \item{R1}{
the response in the first stage, which takes 1 for "Yes", 0 for "No"
}
  \item{R2}{
the response in the second stage, which takes 1 for "Yes", 0 for "No"
}
  \item{bid1}{
the bid in the first stage
}
  \item{bid2}{
the bid in the second stage the respondent faced
}
}

\references{
Aizaki H, Nakatani T, Sato K (2014). Stated Preference Methods Using R. CRC Press, Boca Raton, FL.
}


\seealso{
\code{\link{sbchoice}}, \code{\link{dbchoice}},
}
\examples{
# SBDC-CV
data(CarsonSB)
CarsonSB
CarsonSB[, 1]
SB <- dccvResponse(bid = CarsonSB[, 1], b = c(0.64, -0.011), n = 1000)
head(SB)
sbchoice(R1 ~ 1|bid1, dist = "logistic", data = SB)

# DBDC-CV
data(CarsonDB)
CarsonDB
CarsonDB[, 1:3]
DB <- dccvResponse(bid = data.matrix(CarsonDB[, 1:3]), b = c(0.7, -0.016), n = 1000)
head(DB)
dbchoice(R1 + R2 ~ 1|bid1 + bid2, dist = "logistic", data = DB)
}
\keyword{DCchoice}
\keyword{datasets}
