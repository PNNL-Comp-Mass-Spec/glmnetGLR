\name{sensitivity}
\alias{sensitivity}
\title{Computes classification sensitivity from the confusion matrix summary based on
a set of predicted and truth classes for a signature.}
\usage{
sensitivity(confusionSummary, aggregate = c("micro", "macro"))
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
For each class, we calculate the classification sensitivity
in order to summarize its performance for the signature. We
compute one of two aggregate scores, to summarize the
overall performance of the signature.
}
\details{
We suppose that an observation can be classified into one
(and only one) of K classes. For the jth class (j = 1, ...,
K), we define the sensitivity as the conditional
probability

Sensitivity_j = Pr(y_hat = j | y = j),

where y_hat and y are the empirical and true
classifications, respectively.

To estimate sensitivity_j for the jth class, we compute

(TP_j) / (TP_j + FN_j),

where TP_j and FN_j are the true positives and false
negatives, respectively. More specifically, TP_j is the
number of observations that we correctly classified into
the jth class, and FN_j is is the number of observations
that we should have classified into class j but failed to
do so.

The two aggregate score options are the macro- and
micro-aggregate (average) scores. The macro-aggregate score
is the arithmetic mean of the binary scores for each class.
The micro-aggregate score is a weighted average of each
class' binary score, where the weights are determined by
the sample sizes for each class. By default, we use the
micro-aggregate score because it is more robust, but the
macro-aggregate score might be more intuitive to some
users.

In statistical terms, notice that in the binary case (K =
2), the sensitivity is the recall.

Also, note that the sensitivity is equal to the TPR.

The sensitivity measure ranges from 0 to 1 with 1 being the
optimal value.
}
\examples{
data(prediction_values)

confmat <- confusion(prediction_values[,"Curated_Quality"], prediction_values[,"PredictClass"])

sensitivity(confmat)
}

