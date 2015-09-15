\name{Carson}
\alias{CarsonSB}
\alias{CarsonDB}
\docType{data}
\title{Exxon Valdez Oil Spill CVM data}
\description{
    Contingency tables for the suggested bids and the number of 
    respondents saying yes or no to the bids in the Exxon Valdez 
    Oil Spill CVM survey. 
}

\usage{
data(CarsonSB)
data(CarsonDB)
}

\format{
Both \code{CarsonSB} and \code{CarsonDB} are data frame objects 
of contingency tables. 

For \code{CarsonSB}, 
    \describe{
      \item{T1}{a bid expressd in USD.}
      \item{Y}{a number of respondents accepting the bid.}
      \item{N}{a number of respondents not accepting the bid.}
    }

For \code{CarsonDB}, 
    \describe{
      \item{T1}{a first stage bid expressd in USD.}
      \item{TU}{a second stage bid increased from the first one, expressd in USD.}
      \item{TL}{a second stage bid decreased from the first one, expressd in USD.}
      \item{yy}{a number of respondents accepting both the first and the second bids.}
      \item{yn}{a number of respondents accepting only the first bid.}
      \item{ny}{a number of respondents accepting only the second bid.}
      \item{nn}{a number of respondents not accepting the first nor the second bids.}
    }
}

\details{
Out of \code{CarsonSB} and \code{CarsonSB}, one may reconstruct the original yes/no 
type of data for 1043 observations. See the example for \code{CarsonSB}. 
}

\source{
    \code{CarsonSB} and \code{CarsonDB} are reproduced from Tables A-15, 
    A-16 and A-17 in Appendix C.1 of Carson et.al (1992). 

    The data are included under kind permission from Professor Richard 
    T. Carson of University of California, San Diego. 
}


\references{

Carson RT, Mitchell RC, Hanemann WM, Kopp RJ, Presser S, Ruud PA (1992).
\dQuote{A Contingent Valuation Study of Lost Passive Use Values
  Resulting from the Exxon Valdez Oil Spill.}
\emph{Technical Report Report to the Attorney General of the State of
  Alaska}, Natural Resource Damage Assessment Inc.
\url{http://mpra.ub.uni-muenchen.de/6984/}.

}

\seealso{
\code{\link{ct2df}}
}

\examples{
## The following lines of code reconstruct the original yes/no type of data 
## for 1043 observations. A data frame object sb.data consists of two variables, 
## namely, bid1 and R1. The conversion into a simole data frame object can be 
## done either manually or by using the \code{ct2df} function. 
data(CarsonSB)

## Using the \code{ct2df} function
CarsonSB.dat <- ct2df(
  x    = CarsonSB,
  bid1 = "T1",
  y    = "Y",
  n    = "N",
  type = "single")
head(CarsonSB.dat)

# Manual conversion
n <- rowSums(CarsonSB[, -1])
sb.data <- data.frame(
  bid = c(rep(CarsonSB$T1[1], n[1]), 
          rep(CarsonSB$T1[2], n[2]), 
          rep(CarsonSB$T1[3], n[3]), 
          rep(CarsonSB$T1[4], n[4])), 
  R1 = c(rep(1, CarsonSB$Y[1]), rep(0, CarsonSB$N[1]),
         rep(1, CarsonSB$Y[2]), rep(0, CarsonSB$N[2]),
         rep(1, CarsonSB$Y[3]), rep(0, CarsonSB$N[3]),
         rep(1, CarsonSB$Y[4]), rep(0, CarsonSB$N[4]))
)
dim(sb.data)
head(sb.data)

## Double-bounded dichotomous choice CV format.
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

## An example of manual conversion is omitted. 
## See Appendix 2.A of Aizaki, et. al. (2014). 

}


\keyword{DCchoice}
\keyword{datasets}
