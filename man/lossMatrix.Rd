\name{lossMatrix}
\alias{lossMatrix}
\alias{print.lossMat}
\title{Build a loss matrix}
\usage{
lossMatrix(truthLabels, predLabels, lossValues)

\method{print}{lossMat}(x, ...)
}
\arguments{
  \item{truthLabels}{character vector of truth labels}

  \item{predLabels}{character vector of corresponding
  predicted labels, which must be the same length as
  \code{truthLabels}}

  \item{lossValues}{numeric vector of corresponding loss
  values, which must be the same length as
  \code{truthLabels} and \code{predLabels}.}

  \item{x}{An object of class \code{lossMat}}

  \item{\dots}{Additional arguments to
  \link{\code{print.default}}}
}
\value{
An object of class \code{lossMat}: a dataframe that
contains all the information of the loss matrix, to be used
by \link{\code{calcLoss}}
}
\description{
Build an arbitrary loss matrix for discrete classification
}
\details{
This function checks the inputs and binds the three
arguments columnwise into a dataframe.
}
\examples{
# A 2x2 symmetric loss matrix
lossMatrix(c("a","a","b","b"), c("a","b","a","b"), c(0, 1, 5, 0))

# An unbalanced loss matrix (with a missing element)
lossMatrix(c("a","a","b"), c("a","b","b"), c(0, 1, 0))

# A 3x2 asymmetric loss matrix
lossMatrix(rep(letters[1:3], each = 2), rep(letters[4:5], 3),
           c(0, 3, 2, 0, 1, 0))
}
\author{
Landon Sego
}

