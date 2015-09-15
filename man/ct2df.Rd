\name{ct2df}
\alias{ct2df}
\title{
Convert a data frame in contingency-table format into a simple data frame of individual observations
}
\description{
A convinience function converting a data frame object in contingency-table format of bid(s) and responses of dichotomous choice CV into a simple data frame of individual observations. The outcome is suitable for the analysis using functions in the package.
}
\usage{
ct2df(x, bid1 = "bid1", bid2h = "bidh", bid2l = "bidl", 
    yy = "yy", yn = "yn", ny = "ny", nn = "nn", y = "y", n = "n", type = "double")
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
  \item{yy}{
a character string showing a number of respondents accepting both the first and the second bids
}
  \item{yn}{
a character string showing a number of respondents accepting only the first bid
}
  \item{ny}{
a character string showing a number of respondents accepting only the second bid
}
  \item{nn}{
a character string showing a number of respondents not accepting the first nor the second bids
}
  \item{y}{
a character string showing a number of respondents accepting the bid
}
  \item{n}{
a character string showing a number of respondents not accepting the bid
}
  \item{type}{
a character string setting the elicitation format, which takes one of "single" (single-bounded dichotomous choice format) or "double" (double-bounded dichotomous choice format)
}
}
\details{
The function \code{\link{ct2df}} implements a conversion of a data frame 
containing bid(s) and responses regarding dichotomous choice CV in 
contingency-table format into a data frame suitable for use by the functions 
\code{\link{sbchoice}}, \code{\link{dbchoice}}, \code{\link{kristrom}}, 
\code{\link{turnbull.sb}}, and \code{\link{turnbull.db}}. 
See \code{\link{CarsonSB}} and \code{\link{CarsonDB}} for dataset 
in contingency-table format. Arguments from \code{bid2h} to \code{nn} are 
only valid for double-bounded dichotomous choice format, while arguments 
\code{y} and \code{n} are only valid for single-bounded dichotomous choice 
format.

See the examples, for usage in detail.
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
}
\references{
Aizaki H, Nakatani T, Sato K (2014). Stated Preference Methods Using R. CRC Press, Boca Raton, FL.
}


\seealso{
\code{\link{sbchoice}}, \code{\link{dbchoice}}, \code{\link{kristrom}}, \code{\link{turnbull.sb}}, \code{\link{turnbull.db}}
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
