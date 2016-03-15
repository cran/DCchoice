\name{update.dbchoice}
\alias{update.dbchoice}

\title{
Updating and refitting model for dbchoice
}

\description{
Updating and refitting method for object of class "dbchoice".
}

\usage{
\method{update}{dbchoice}(object, new, evaluate = TRUE, ...)
}

\arguments{
  \item{object}{
an object of class "dbchoice".
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
The function \code{update()} for S3 object "dbchoice" updates a model used for 
the fit that is included in \code{object} according to a new call assigned to 
\code{new}, and then refits the updated model. The function returns the refitted 
model object when \code{evaluate = TRUE}, otherwise the updated model call.
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
