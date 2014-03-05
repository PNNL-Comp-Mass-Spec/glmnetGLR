\name{confusion}
\alias{confusion}
\title{Computes the confusion summary for a vector of classifications and a ground
truth vector.}
\usage{
confusion(truthClass, predictedClass)
}
\arguments{
  \item{truthClass}{vector of ground truth classification
  labels}

  \item{predictedClass}{vector of predicted classification
  labels}
}
\value{
list with the results of confusion matrix results for each
class.
}
\description{
For a vector of classifications and truth labels, we create
a confusion matrix. We allow binary and multi-class
classifications and compute the following four measures for
each class:
}
\details{
\itemize{ \item True Positives (TP) \item True Negatives
(TN) \item False Positives (FP) \item False Negatives (FN)
}

For multi-class classification, we consider each class in a
binary context. For example, suppose that we have the three
food condiment classes: ketchup, mustard, and other. When
calculating the TP, TN, FP, and FN values for ketchup, we
consider each observation as either 'ketchup' or 'not
ketchup.' Similarly, for mustard, we would consider
'mustard' and 'not mustard', and for other, we would
consider 'other' and 'not other.'

With the above counts for each class, we can quickly
calculate a variety of class-specific and aggregate
classification accuracy measures.
}
\examples{
data(prediction_values)

confusion(prediction_values[,"Curated_Quality"], prediction_values[,"PredictClass"])
}

