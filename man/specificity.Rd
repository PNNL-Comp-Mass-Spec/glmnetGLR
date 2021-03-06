\name{specificity}
\alias{specificity}
\title{Computes classification specificity from the confusion matrix summary based on
a set of predicted and truth classes for a signature.}
\usage{
specificity(confusionSummary, aggregate = c("micro", "macro"))
}
\arguments{
  \item{confusionSummary}{list containing the confusion
  summary for a set of classifications}

  \item{aggregate}{string that indicates the type of
  aggregation; by default, micro. See details.}
}
\value{
list with the accuracy measure for each class as well as
the macro- and micro-averages (aggregate measures across
all classes).
}
\description{
For each class, we calculate the classification specificity
in order to summarize its performance for the signature. We
compute one of two aggregate scores, to summarize the
overall performance of the signature.
}
\details{
To estimate specificity_j for the jth class, we compute

(TN_j) / (TN_j + FP_j),

where TN_j and FP_j are the true negatives and false
positives, respectively. More specifically, TN_j is the
number of observations that we correctly classified into
other classes than the jth class, and FP_j is the number of
observations that we have incorrectly classified into class
j.

The two aggregate score options are the macro- and
micro-aggregate (average) scores. The macro-aggregate score
is the arithmetic mean of the binary scores for each class.
The micro-aggregate score is a weighted average of each
class' binary score, where the weights are determined by
the sample sizes for each class. By default, we use the
micro-aggregate score because it is more robust, but the
macro-aggregate score might be more intuitive to some
users.

Notice that the specificity is equal to the TNR.

The specificity measure ranges from 0 to 1 with 1 being the
optimal value.
}
\examples{
data(prediction_values)

confmat <- confusion(prediction_values[,"Curated_Quality"], prediction_values[,"PredictClass"])

specificity(confmat)
}

