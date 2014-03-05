\name{summary.glmnetGLRpred}
\alias{summary.glmnetGLRpred}
\title{Calculate a set of metrics to assess the performance of the LLRC}
\usage{
\method{summary}{glmnetGLRpred}(object, ...)
}
\arguments{
  \item{object}{the object produced by \code{predict}}

  \item{...}{additional arguments to other methods
  (currently ignored)}
}
\value{
\code{data.frame} with sensitivity, specificity, false
negative rate, false positive rate, and accuracy computed
by class.
}
\description{
This function allows for a quick computation of multiple
metrics to assess the performance of the elastic net
classififer assuming the ground truth is available.
}
\note{
The ground truth for the new data set is required to
generate these metrics.
}
\author{
Alex Venzin
}

