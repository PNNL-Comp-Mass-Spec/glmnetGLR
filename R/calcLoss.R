##' Calculate the loss of a classifier using discrete loss matrix
##'
##' Calculates the loss of a classifer using a discrete loss matrix,
##' where each observation may be weighted as desired.
##'
##' @export
##' 
##' @param truthLabels  A factor vector containing the true classes for each
##' instance (or observation)
##'
##' @param predClass A factor vector containing the predicted classes for
##' each instance (or observation)
##'
##' @param lossMat A object of class \code{lossMat} (constructed using
##' \link{\code{lossMatrix}}) that indicates the loss for all possible
##' classifications
##'
##' @param weight A non-negative numeric vector of weights that will
##' be used to weight each observation in calculating the loss.
##'
##' @param aggregate A logical indicating whether the aggregate loss
##' or the loss of individual observations should be returned
##' (see the 'Value' section).
##'
##' @return The loss, either as a single numeric value caculated as the
##' weighted average of the individual loss values
##' (\code{aggregate = TRUE}) or
##' as a dataframe showing the loss for each observation
##' (\code{aggregate = FALSE}).
##'
##' @author Landon Sego
##' 
##' @examples
##'
##' # Build the loss matrix
##' lMat <- lossMatrix(rep(letters[1:3], 3), rep(letters[1:3], each = 3),
##'                    c(0, 1, 2, 1, 0, 1, 2, 1, 0))
##' lMat
##' 
##' # Create a vector of labels, simulating instances
##' tClass <- factor(rep(letters[1:3], each = 5))
##' pClass <- sample(tClass)
##' 
##' 
##' # Calculate the loss 
##' calcLoss(tClass, pClass, lMat, aggregate = FALSE)
##' calcLoss(tClass, pClass, lMat)

calcLoss <- function(truthLabels, predLabels, lossMat,
                     weight = rep(1, length(truthLabels)),
                     aggregate = TRUE) {

  # Checks
  stopifnot(is.factor(truthLabels),
            is.factor(predLabels),
            inherits(lossMat, "lossMat"),
            length(truthLabels) == length(predLabels),
            length(truthLabels) == length(weight),
            is.numeric(weight),
            sum(weight) > 0,
            all(weight >= 0))

  # Convert the factors to characters and paste together
  tpData <- data.frame(truthLabels = truthLabels,
                       predLabels = predLabels,
                       weight = weight,
                       index = 1:length(truthLabels))

  # Assign the loss values to the corresponding predicted and actual labels
  tpD <- merge(tpData, lossMat, all.x = TRUE, sort = FALSE)

  # Verify there weren't any data for which we couldn't calculate the loss
  if (!all(cc <- complete.cases(tpD))) 
    warning("One or more observations could not be matched to a coresponding\n",
            "truth/predicted label pair in the loss matrix")

  # Retun the loss as requested
  if (!aggregate) {
    return(tpD[order(tpD$index),
               c("truthLabels", "predLabels", "weight", "loss")])
  }

  else {
    tpDcc <- tpD[cc, c("loss", "weight")]
    return(list(weightedSumLoss = sum(tpDcc$loss * tpDcc$weight),
                sumWeights = sum(tpDcc$weight)))
  }
           

} # calcLoss
