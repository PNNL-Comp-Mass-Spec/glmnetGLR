\name{accuracy}
\alias{accuracy}
\title{Computes classification accuracy from the confusion matrix summary based on a
set of predicted and truth classes for a signature.}
\usage{
accuracy(confusionSummary, aggregate = c("micro", "macro"))
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
For each class, we calculate the classification accuracy in
order to summarize its performance for the signature. We
compute one of two aggregate scores, to summarize the
overall performance of the signature.
}
\details{
We define the accuracy as the proportion of correct
classifications.

The two aggregate score options are the macro- and
micro-aggregate (average) scores. The macro-aggregate score
is the arithmetic mean of the binary scores for each class.
The micro-aggregate score is a weighted average of each
class' binary score, where the weights are determined by
the sample sizes for each class. By default, we use the
micro-aggregate score because it is more robust, but the
macro-aggregate score might be more intuitive to some
users.

Note that the macro- and micro-aggregate scores are the
same for classification accuracy.

The accuracy measure ranges from 0 to 1 with 1 being the
optimal value.
}
\examples{
data(prediction_values)

# Create the confusion matrix
confmat <- confusion(prediction_values[,"Curated_Quality"], prediction_values[,"PredictClass"])

accuracy(confmat)
}

