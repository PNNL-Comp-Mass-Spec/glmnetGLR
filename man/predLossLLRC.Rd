\name{predLossLLRC}
\alias{predLossLLRC}
\title{Prediction function for the LLRC for a variety of thresholds and calculate
the associated loss}
\usage{
predLossLLRC(glmnetFit, newData, truthLabels, lossMat, tauVec = seq(0.1, 0.9,
  by = 0.1), weight = rep(1, NROW(newData)))
}
\arguments{
  \item{glmnetFit}{The glmnet fitted object (returned from
  \link{\code{glmnet}} that inherits from the \code{lognet}
  and \code{glmnet} classes.}

  \item{newData}{A matrix or dataframe of new instances
  that match the predictors in \code{glmnetFit}}

  \item{truthDataLabels}{A factor vector containing the
  corresponding truth labels in \code{newData}.}

  \item{lossMat}{A loss matrix of class \code{lossMat},
  returned by \link{\code{lossMatrix}}}

  \item{tauVec}{A numeric sequence of threshold values for
  the binary classification.}

  \item{weight}{A numeric vector indicating the relative
  weight to ascribe to each row of \code{newData}}
}
\value{
A data frame containing \code{weightedSumLoss} (the sum of
the product of the weights and the loss) and
\code{sumWeights} (the sum of the weights) for each value
of \code{tau} and \code{lambda}.  The loss (for
\code{newData}) is given by \code{weightedSumLoss} divided
by \code{sumWeights}.
}
\description{
Let \code{z} represent class denoted by the last level of
\code{truthLabels}. Then the probability returned by
\code{predict(glmnetFit, newData, type = 'response')} is
the probability that the observations belong to class z.
If this probabilty exceeds \code{tau}, we will classify the
observation as belonging to \code{z}.
}
\author{
Landon Sego
}

