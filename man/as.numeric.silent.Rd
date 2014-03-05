\name{as.numeric.silent}
\alias{as.numeric.silent}
\title{Silent wrapper for coercing a vector to numeric}
\usage{
\method{as.numeric}{silent}(x)
}
\arguments{
  \item{x}{vector of any type}
}
\value{
If \code{as.numeric(x)} produces an error or warning,
\code{x} is returned unchanged.  Otherwise,
\code{as.numeric(x)} is returned.
}
\description{
Silent wrapper for coercing a vector to numeric
}
\examples{
as.numeric.silent(c("this","that"))
as.numeric.silent(c("2893.9","9423.48"))
}
\seealso{
\code{\link{as.numeric}}
}

