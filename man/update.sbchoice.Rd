\name{update.sbchoice}
\alias{update.sbchoice}

\title{
Updating and refitting model for sbchoice
}

\description{
Updating and refitting method for object of class "sbchoice".
}

\usage{
\method{update}{sbchoice}(object, new, evaluate = TRUE, ...)
}

\arguments{
  \item{object}{
an object of class "sbchoice".
}
  \item{new}{
a new call.
}
  \item{evaluate}{
If \code{TRUE}, refit the updated model, otherwise return the updated call.
}
  \item{...}{
optional arguments. Currently not in use.
}
}

\details{
The function \code{update()} for S3 object "sbchoice" updates a model used for 
the fit that is included in \code{object} according to a new call assigned to 
\code{new}, and then refits the updated model. The function returns the refitted 
model object when \code{evaluate = TRUE}, otherwise the updated model call.
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
