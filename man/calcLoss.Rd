\name{calcLoss}
\alias{calcLoss}
\title{Calculate the loss of a classifier using discrete loss matrix}
\usage{
calcLoss(truthLabels, predLabels, lossMat, weight = rep(1,
  length(truthLabels)), aggregate = TRUE)
}
\arguments{
  \item{truthLabels}{A factor vector containing the true
  classes for each instance (or observation)}

  \item{predClass}{A factor vector containing the predicted
  classes for each instance (or observation)}

  \item{lossMat}{A object of class \code{lossMat}
  (constructed using \link{\code{lossMatrix}}) that
  indicates the loss for all possible classifications}

  \item{weight}{A non-negative numeric vector of weights
  that will be used to weight each observation in
  calculating the loss.}

  \item{aggregate}{A logical indicating whether the
  aggregate loss or the loss of individual observations
  should be returned (see the 'Value' section).}
}
\value{
The loss, either as a single numeric value caculated as the
weighted average of the individual loss values
(\code{aggregate = TRUE}) or as a dataframe showing the
loss for each observation (\code{aggregate = FALSE}).
}
\description{
Calculates the loss of a classifer using a discrete loss
matrix, where each observation may be weighted as desired.
}
\examples{
# Build the loss matrix
lMat <- lossMatrix(rep(letters[1:3], 3), rep(letters[1:3], each = 3),
                   c(0, 1, 2, 1, 0, 1, 2, 1, 0))
lMat

# Create a vector of labels, simulating instances
tClass <- factor(rep(letters[1:3], each = 5))
pClass <- sample(tClass)


# Calculate the loss
calcLoss(tClass, pClass, lMat, aggregate = FALSE)
calcLoss(tClass, pClass, lMat)
}
\author{
Landon Sego
}

