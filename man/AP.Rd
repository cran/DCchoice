\name{AP}
\alias{AP}
\docType{data}
\title{Albemarle-Pamlico sounds CVM data}
\description{
    Double-bounded dichotomous choice survey data for quality inprovements in 
    the Albemarle-Pamlico sounds, North Carolina.
}

\usage{
data(AP)
}

\format{
a data frame containing 721 observations.

    \describe{
      \item{bid1}{a vector of bids expressd in USD.}
      \item{bid2}{a vector of bids expressd in USD.}
      \item{R1}{a vector of binary dummies equal to 1 if the bid is accepted and 0 otherwise.}
      \item{R2}{a vector of binary dummies equal to 1 if the bid is accepted and 0 otherwise.}
      \item{income}{a numeric vrector containing the annual household income in 1995 USD of the respondent}
      \item{work}{a vector of binary dummies equal to 1 if the respondent is employed in fulltime and 0 otherwise.}
      \item{age}{a numeric vector containing the age of the respondent.}
      \item{female}{a vector of binary dummies equal to 1 if the respondent is female and 0 otherwise.}
      \item{married}{a vector of binary dummies equal to 1 if the respondent is married and 0 otherwise.}
    }
}

\details{
The original data are based on a telephone survey regarding quality improvements 
in the Albemarle-Pamlico sounds, North Carolina. The data were intensively analyzed, 
for instance, in Whitehead (1995) and Whitehead, et. al. (1998) in different contexts. Details of the survey can be 
found in the web site (see \emph{Source} for the URL). 

The original data have $1077$ observations and include the bids and the responses 
of the double-bounded dichotomous choice survey. Various socio-demographic 
characteristics are also collected by the survey. 

A subset of the data consisting of the responses to the CVM questions 
as well as minimum number of socio-demographic characteristics. In addition, 
observations with missing values were removed from the subset. 
}

\source{
The complete data and details can be obtained in the online version of 
Whitehead (2015) at 
\url{http://dx.doi.org/10.1016/j.dib.2015.01.006}.

The data are included here under kind permission from Professor 
John Whitehead, Appalachian State University.
}

\references{
Whitehead JC (1995).
\dQuote{Willingness to Pay for Quality Improvements: Comparative
  Statics and Interpretation of Contingent Valuation Results.}
\emph{Land Economics}, \bold{71}(2), 207--215.

Whitehead JC (2015).
\dQuote{Albemarle-Pamlico Sounds Revealed and Stated Preference Data.}
\emph{Data in Brief}, \bold{3}, 90--94.
\url{http://dx.doi.org/10.1016/j.dib.2015.01.006}

Whitehead JC, Haab TC, Huang JC (1998).
\dQuote{Part-Whole Bias in Contingent Valuation: Will Scope Effects 
be Detected with Inexpensive Survey Methods?}
\emph{Southern Economic Journal}, \bold{65}(1), 160--168.
}

\keyword{DCchoice}
\keyword{datasets}
