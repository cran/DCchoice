\name{predict.dbchoice}
\alias{predict.dbchoice}

\title{
Predicting model for dbchoice
}

\description{
Predicting method for objects of class "dbchoice".
}

\usage{
\method{predict}{dbchoice}(object, newdata = NULL, type = c("utility", "probability"), 
bid = NULL, ...)
}

\arguments{
  \item{object}{
an object of class "dbchoice".
}
  \item{newdata}{
a data frame containing new data to predict. If \code{NULL}, the original data is used.
}
  \item{type}{
type of prediction (utility or probability).
}
  \item{bid}{
a bid design needed to predict with original data.
}
  \item{\dots}{
optional arguments. Currently not in use.
}
}

\details{
The function \code{predict()} for S3 object "dbchoice" calculates predicted values 
according to the fitted model that is included in \code{object}. The values are 
predicted with the original data used for fitting the model if \code{newdata = NULL}, 
otherwise with a new data assigned to \code{newdata}. There are two notes for 
\code{dbchoice()}: a bid design used for the fit must be assigned to \code{bid}, 
when predicting with the original data; the predicted values are calculated according 
to the relevant single-bounded dichotomous choice model, when predicting with a new 
data. See examples for details.

The current function does not estimate standard errors of predicted values.
}

\value{
When \code{newdata = NULL} and \code{type = utility}, a matrix containing utility 
values under first (\code{f}), second upper (\code{u}), and second lower (\code{l}) bids 
is returned. When \code{newdata = NULL} and \code{type = probability}, a matrix 
containing probabilities of choosing Yes-Yes (\code{yy}), No-No (\code{nn}), Yes-No (\code{yn}), 
and No-Yes (\code{ny}) is returned. When a new data is assigned to \code{newdata},
predictions are calculated according to the relevant single-bounded dichotomous choice 
model, and a vector containing utility values of choosing yes (\code{type = utility}) 
or probability of choosing yes (\code{type = probability}) is returned.
}

\seealso{
\code{\link{dbchoice}}
}

\examples{
## See Examples in dbchoice.
}

\keyword{DCchoice}
\keyword{double-bounded}
\keyword{methods}
