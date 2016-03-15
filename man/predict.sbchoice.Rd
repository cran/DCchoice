\name{predict.sbchoice}
\alias{predict.sbchoice}

\title{
Predicting model for sbchoice
}

\description{
Predicting method for objects of class "sbchoice".
}

\usage{
\method{predict}{sbchoice}(object, newdata = NULL, type = c("utility", "probability"), 
...)
}

\arguments{
  \item{object}{
an object of class "sbchoice".
}
  \item{newdata}{
a data frame containing new data to predict. If \code{NULL}, the original data is used.
}
  \item{type}{
type of prediction (utility or probability).
}
  \item{\dots}{
optional arguments. Currently not in use.
}
}

\details{
The function \code{predict()} for S3 object "sbchoice" calculates predicted values 
according to the fitted model that is included in \code{object}. The values are 
predicted with the original data used for fitting the model if \code{newdata = NULL}, 
otherwise with a new data assigned to \code{newdata}.

The current function does not estimate standard errors of predicted values.
}

\value{
When \code{newdata = NULL} and \code{type = utility}, a vector containing utility 
values of choosing Yes under bid values is returned. When \code{newdata = NULL} 
and \code{type = probability}, a vector containing probabilities of choosing Yes
is returned. When a new data is assigned to \code{newdata}, a vector containing 
utility values of choosing Yes (\code{type = utility}) or probability of 
choosing Yes (\code{type = probability}) under the new data is returned.
}

\seealso{
\code{\link{sbchoice}}
}

\examples{
## See Examples in sbchoice.
}

\keyword{DCchoice}
\keyword{single-bounded}
\keyword{methods}
