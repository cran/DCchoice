\name{ct2df}
\alias{ct2df}
\title{
Convert a data frame in contingency-table format into a simple data frame of individual observations
}
\description{
A convinience function converting a data frame object in contingency-table format of bid(s) and responses of dichotomous choice CV into a simple data frame of individual observations. The outcome is suitable for the analysis using functions in the package.
}
\usage{
ct2df(x,
  bid1 = "bid1", bid2h = "bidh", bid2l = "bidl", bidl = "bidl", bidh = "bidh",
  nny = "nny", nnn = "nnn", yy = "yy", yn = "yn", ny = "ny", nn = "nn",
  y = "y", n = "n", n_y = "n_y", n_n = "n_n", nn_y = "nn_y", nn_n = "nn_n", 
  type = "double", spike = FALSE)
}
\arguments{
  \item{x}{
a data frame object in contingency-table format containing bid(s) and responses
}
  \item{bid1}{
a character string showing the bid (for "single") or the bid in the first stage (for "double")
}
  \item{bid2h}{
a character string showing the second (higher) bid when respondents answer "Yes" in the first stage
}
  \item{bid2l}{
a character string showing the second (lower) bid when respondents answer "No" in the first stage
}
  \item{bidh}{
a character string showing the higher bid (for "oohb")
}
  \item{bidl}{
a character string showing the lower bid (for "oohb")
}
  \item{nny}{
a character string showing a number of respondents rejecting the first and the second bids and having a positive willingness-to-pay (WTP) (for a spike "double")
}
  \item{nnn}{
a character string showing a number of respondents rejecting the first and the second bids and having a zero WTP (for a spike "double")
}
  \item{yy}{
a character string showing a number of respondents accepting both the first and the second bids
}
  \item{yn}{
a character string showing a number of respondents accepting only the first bid
}
  \item{ny}{
a character string showing a number of respondents accepting only the second bid (for "double") or rejecting the first bid and having a positive WTP (for a spike "single")
}
  \item{nn}{
a character string showing a number of respondents rejecting the first and the second bids (for "double") or rejecting the first bid and having a zero WTP (for a spike "single")
}
  \item{y}{
a character string showing a number of respondents accepting the bid
}
  \item{n}{
a character string showing a number of respondents rejecting the bid
}
  \item{n_y}{
a character string showing a number of respondents rejecting the bid and having a positive WTP (for a spike "oohb")
}
  \item{n_n}{
a character string showing a number of respondents rejecting the bid and having a zero WTP (for a spike "oohb")
}
  \item{nn_y}{
a character string showing a number of respondents rejecting the first and second bids and having a positive WTP (for a spike "oohb")
}
  \item{nn_n}{
a character string showing a number of respondents rejecting the first and second bid and having a zero WTP (for a spike "oohb")
}
  \item{type}{
a character string setting the elicitation format, which takes one of "single" (single-bounded dichotomous choice format), "oohb" (one-and-one-half-bounded dichotomous choice format), or "double" (double-bounded dichotomous choice format)
}
  \item{spike}{
a logical code: \code{TRUE} for a spike model. The default is \code{FALSE}.
}
}
\details{
The function \code{\link{ct2df}} implements a conversion of a data frame 
containing bid(s) and responses regarding dichotomous choice CV in 
contingency-table format into a data frame suitable for use by the functions 
\code{\link{sbchoice}}, \code{\link{oohbchoice}}, \code{\link{dbchoice}}, \code{\link{sbspike}}, \code{\link{oohbspike}}, \code{\link{dbspike}}, \code{\link{kristrom}}, 
\code{\link{turnbull.sb}}, and \code{\link{turnbull.db}}. 
See \code{\link{CarsonSB}} and \code{\link{CarsonDB}} for dataset 
in contingency-table format.

See the examples below and ones in \code{\link{spike-models}}, for usage in detail.
}
\value{
The function returns a data frame, in which each row shows a single respondent. It contains the following variables.

For "single",
  \item{R1}{
a response to a bid: 1 for "Yes", 0 for "No"
}
  \item{bid1}{
the bid
}

For "double",
  \item{B1}{
a bid in the first stage
}
  \item{B2H}{
a (higher) bid in the second stage when the response is "Yes" in the first stage
}
  \item{B2L}{
a (lower) bid in the second stage when the response is "No" in the first stage
}
  \item{R}{
a combination of responses in the first and second stages, which takes yy for "Yes" and "Yes", yn for "Yes" and "No", ny for "No" and "Yes", or nn for "No" and "No"
}
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

For "oohb",
\item{BH}{a higher bid}
\item{BL}{a lower bid}
\item{R}{a combination of responses in the first and second stages}
\item{R1}{the response in the first stage, which takes \code{1} if the bid is accepted, and \code{0} otherwise}
\item{R2}{the response in the second stage, which takes \code{1} if the bid is accepted, \code{0} if the bid is not accepted, and \code{-9} if the respondent has no second stage}

For spike models,
\item{S}{the response in a simple spike question, which takes \code{1} if the respondent has a positive WTP, and \code{0} otherwise}

}
\references{
Aizaki H, Nakatani T, Sato K (2014). Stated Preference Methods Using R. CRC Press, Boca Raton, FL.
}


\seealso{
\code{\link{sbchoice}}, \code{\link{oohbchoice}}, \code{\link{dbchoice}},
\code{\link{sbspike}}, \code{\link{oohbspike}}, \code{\link{dbspike}},
\code{\link{kristrom}}, \code{\link{turnbull.sb}}, \code{\link{turnbull.db}}
}
\examples{
# Single-bounded dichotomous choice CV format
data(CarsonSB)
CarsonSB
CarsonSB.dat <- ct2df(
  x    = CarsonSB,
  bid1 = "T1",
  y    = "Y",
  n    = "N",
  type = "single")
head(CarsonSB.dat)
summary(turnbull.sb(R1 ~ bid1, data = CarsonSB.dat))

# Double-bounded dichotomous choice CV format
data(CarsonDB)
CarsonDB
CarsonDB.dat <- ct2df(
  x     = CarsonDB,
  bid1  = "T1",
  bid2h = "TU",
  bid2l = "TL",
  yy    = "yy",
  yn    = "yn",
  ny    = "ny",
  nn    = "nn",
  type  = "double")
head(CarsonDB.dat)
summary(turnbull.db(R1 + R2 ~ bid1 + bid2, data = CarsonDB.dat))
}
\keyword{DCchoice}
\keyword{datasets}
