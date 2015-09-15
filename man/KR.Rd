\name{KR}
\alias{KR}
\docType{data}
\title{Kristr\"{o}m's single-bounded dichotomous choice CVM data}
\description{
    A single-bounded dichotomous choice CVM data analyzed in
    Kristr\"{o}m (1990).
}

\usage{
data(KR)
}

\format{
    A data frame of single-bounded dichotomous choice contingent 
    valuation survey data. The number of observations is 562. 

    \describe{
      \item{bid1}{a vector of bids expressd in SEK.}
      \item{R1}{a vector of binary dummies equal to 1 if the bid is accepted and 0 otherwise.}
    }
}

\details{
    The data consist of the responses to the single-bounded dichotomous 
    choice survey for a sample of 562 Swedes regarding preservation of 
    the eleven virgin forests in Sweden. See Kristr\"{o}m (1990) 
    for more details.
}

\source{
    The data are used in Kristr\"{o}m (1990). 
    
    The data are bundled in this package under kind permission from 
    Professor Bengt Kristr\"{o}m, Swedish University of Agricultural Sciences.
}


\references{

Kristr\"{o}m B (1990). 
\dQuote{A Non-Parametric Approach to the Estimation of Welfare 
Measures in Discrete Response Valuation Studies.}
\emph{Land Economics}, \bold{66}(2), 135--139.

}


\keyword{DCchoice}
\keyword{datasets}
