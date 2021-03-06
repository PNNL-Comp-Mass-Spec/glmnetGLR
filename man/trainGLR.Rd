\name{trainGLR}
\alias{predict.glmnetGLR}
\alias{print.glmnetGLR}
\alias{trainGLR}
\title{Train an elastic net logistic regression classifier}
\usage{
trainGLR(truthLabels, predictors, lossMat, weight = rep(1, NROW(predictors)),
  alphaVec = seq(0, 1, by = 0.2), tauVec = seq(0.1, 0.9, by = 0.05),
  naFilter = 0.6, cvFolds = 5, seed = 1, verbose = FALSE)

\method{print}{glmnetGLR}(glmnetGLRobject, ...)

\method{predict}{glmnetGLR}(glmnetGLRobject, newdata, truthCol = NULL,
  keepCols = NULL)
}
\arguments{
  \item{truthLabels}{the factor/character regression
  response variable}

  \item{predictors}{matrix whose columns are the
  explanatory regression variables}

  \item{lossMat}{a loss matrix specifying the penalties for
  classification errors}

  \item{weight}{the observation weights. The default value
  is 1 for each observation. Refer to \code{glmnet} for
  further information}

  \item{alphaVec}{The elastic net mixing parameter. Refer
  to \code{glmnet} for further information.}

  \item{tauVec}{A sequence of tau threshold values for the
  logistic regression classifier.}

  \item{naFilter}{the proportion of data for a given
  predictor (column) that does not consist of missing data.
  If the proportion of sample observations which are NA >
  naFilter, then those predictors are removed from the
  regression analysis.}

  \item{cvFolds}{the number of cross validation folds}

  \item{seed}{the random seed for sampling during cross
  validation}

  \item{verbose}{set to \code{TRUE} to receive messages
  regarding the progress of the training algorithm.}

  \item{glmnetGLRobject}{a train elastic net logistic
  classifier}

  \item{\dots}{Additional arguments to methods (ignored)}

  \item{glmnetGLRobject}{a trained elastic net logistic
  regression classifier}

  \item{newdata}{the new set of observations to be
  predicted. New data must be an array with the same column
  names as the training data}

  \item{\dots}{Additional arguments to predict method of
  \code{glmnet} objects are ignored, as we use the optimal
  values of \code{lambda}, \code{tau}, and \code{alpha}
  that were generated by \code{trainGLR()}.}
}
\value{
The elastic net logistic regression model that minimized
the expected loss over the set of \code{alpha},
\code{lambda}, and \code{tau} parameters.
}
\description{
Train an elastic net logistic regression classifier
}
\examples{
# Load the VOrbitrap Shewanella QC data
data(traindata)

# Here we select the predictor variables
predictors <- as.matrix(traindata[,9:96])

# The logistic regression model requires a binary response
# variable.
resp <- traindata[,"response"]

# Specify the loss matrix. The "Poor" class is the target of interest.
# The penalty for misclassifying a "Poor" item as "Good" results in a
# loss of 5.

lM <- lossMatrix(c("Good","Good","Poor","Poor"),
                c("Good","Poor","Good","Poor"),
                c(     0,     1,     5,     0))

# Train the elastic net classifier
elasticNet <- trainGLR(truthLabels = resp,
                       predictors = predictors,
                       lossMat = lM)

# Observe the optimal alpha, lambda, and tau values that produced
# this elastic net logistic regression classifier.
print(elasticNet)

# Load the new observations
data(testdata)

# Use an elastic net regression classifier to make predictions about
# new observations for the response variable.
predict(elasticNet, testdata)
}
\author{
Landon Sego

Alex Venzin
}

