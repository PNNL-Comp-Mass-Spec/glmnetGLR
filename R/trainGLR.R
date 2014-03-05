##' Train an elastic net logistic regression classifier
##'
##' @author Landon Sego
##' @author Alex Venzin
##'
##' @rdname trainGLR
##'
##' @export
##' 
##' @param truthLabels the factor/character regression
##' response variable
##' 
##' @param predictors matrix whose columns are the
##' explanatory regression variables
##' 
##' @param lossMat a loss matrix specifying the penalties for
##' classification errors 
##' 
##' @param weight the observation weights. The default value
##' is 1 for each observation. Refer to \code{glmnet} for further
##' information
##'
##' @param alphaVec The elastic net mixing parameter. Refer to
##' \code{glmnet} for further information.
##' 
##' @param tauVec A sequence of tau threshold values for the 
##' logistic regression classifier.
##' 
##' @param naFilter the proportion of data for a given predictor
##' (column) that does not consist of missing data. If the proportion of
##' sample observations which are NA > naFilter, then those predictors 
##' are removed from the regression analysis.
##' 
##' @param cvFolds the number of cross validation folds
##' 
##' @param seed the random seed for sampling during cross validation
##' 
##' @param verbose set to \code{TRUE} to receive messages regarding
##' the progress of the training algorithm.
##' 
##' @return The elastic net logistic regression model that minimized
##' the expected loss over the set of \code{alpha},
##' \code{lambda}, and \code{tau} parameters. 
##' 
##' @examples 
##' # Load the VOrbitrap Shewanella QC data
##' data(traindata)
##'
##' # Here we select the predictor variables
##' predictors <- as.matrix(traindata[,9:96])
##'
##' # The logistic regression model requires a binary response
##' # variable.
##' resp <- traindata[,"response"]
##'
##' # Specify the loss matrix. The "Poor" class is the target of interest.
##' # The penalty for misclassifying a "Poor" item as "Good" results in a
##' # loss of 5.
##'
##' lM <- lossMatrix(c("Good","Good","Poor","Poor"),
##'                 c("Good","Poor","Good","Poor"),
##'                 c(     0,     1,     5,     0))
##'
##' # Train the elastic net classifier
##' elasticNet <- trainGLR(truthLabels = resp,
##'                        predictors = predictors,
##'                        lossMat = lM)
##'
##' # Observe the optimal alpha, lambda, and tau values that produced
##' # this elastic net logistic regression classifier.
##' print(elasticNet)
##'
##' # Load the new observations
##' data(testdata)
##'
##' # Use an elastic net regression classifier to make predictions about
##' # new observations for the response variable.
##' predict(elasticNet, testdata)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
trainGLR <- function(truthLabels, predictors, lossMat,
                      weight = rep(1, NROW(predictors)),
                      alphaVec = seq(0, 1, by = 0.2),
                      tauVec = seq(0.1, 0.9, by = 0.05),
                      naFilter = 0.6,
                      cvFolds = 5, 
                      seed = 1,
                      verbose = FALSE) {

  # Checks on inputs
  stopifnot(NCOL(predictors) > 1,
            is.matrix(predictors),
            is.numeric(predictors),
            length(unique(truthLabels)) == 2,
            length(truthLabels) == NROW(predictors),
            is.numeric(alphaVec),
            all(alphaVec <= 1) & all(alphaVec >= 0),            
            is.numeric(tauVec),
            all(tauVec < 1) & all(tauVec > 0),
            length(naFilter) == 1,
            is.numeric(naFilter),
            naFilter < 1 & naFilter > 0,
            length(seed) == 1,
            is.numeric(seed))
  
  # Force the evaluation of the weight object immediately--this is IMPORTANT
  force(weight)
  
  # If truthLabels is not a factor, make it one
  if (!is.factor(truthLabels))

    truthLabels <- as.factor(truthLabels)

  # Remove predictors that have fewer than naFilter% observations
  # TRUE means that predictor has enough data
  selPreds <- apply(predictors, 2,
                    function(x) {
                     (sum(!is.na(x)) / NROW(predictors)) > naFilter
                    })

  if (!all(selPreds)) {

    if (verbose)
      cat("The following predictors were removed because less than\n",
          naFilter * 100, "% of their observations were non-missing:\n",
          "'", paste(colnames(predictors)[!selPreds], collapse = "', '"), "'\n",
          sep = "")
    
    predictors <- predictors[,selPreds]

  }

  # Identify the obs that don't have any missing predictors or responses
  cc <- complete.cases(cbind(truthLabels, predictors))

  # Remove obs that don't have complete cases
  predictors <- predictors[cc,]
  truthLabels <- truthLabels[cc]
  weight <- weight[cc]

  
  if (verbose & (!all(cc)))
    cat(sum(!cc), "observations (rows) were removed because one or more of",
        "their values were missing\n")
  
  
  # number of observations after filtering for NAs
  n <- sum(cc)

  # Parse the observations in predictors into sets for cross validation
  trainFolds <- lapply(parseJob(n, cvFolds, random.seed = seed), sort)

  # Function to train and test over the folds
  trainTest <- function(x, a = 1, lambdaV = NULL) {

    # x is the set of traning indexes from 'trainFolds'
    # a is alpha
    glmnetFit <- glmnet(predictors[x,], truthLabels[x], weights = weight[x],
                        family = "binomial", lambda = lambdaV, alpha = a)

    # Find the complement of the training set
    testSet <- sort(setdiff(1:n, x))

    # Now test it on the complement of the observations
    out <- predLossLLRC(glmnetFit, predictors[testSet,], truthLabels[testSet],
                        lossMat, tauVec = tauVec, weight = weight[testSet])

    # Update the cvFold
#    if (verbose) {
#      pvar(cvFold)
#      cvFold <<- cvFold + 1
#    }

    return(out)

  } # trainTest
  
  # Run the cross validation for a particular alpha
  cvForAlpha <- function(alpha) {

    if (verbose) {
      pvar(alpha)
#      cvFold <<- 1
    }

    # Get the lambdaVec for this particular alpha using all the data
    lambdaVec <- glmnet(predictors, truthLabels, weights = weight,
                        family = "binomial", alpha = alpha)$lambda

    # Run the cross validation
    if (.Platform$OS.type != "windows")
      testAll <- list2df(mclapply(trainFolds, trainTest, a = alpha,
                                  lambdaV = lambdaVec,
                                  mc.cores = min(cvFolds, detectCores() - 1)))
    else
      testAll <- list2df(lapply(trainFolds, trainTest, a = alpha,
                                lambdaV = lambdaVec))
    
    # Add in the alpha
    testAll$alpha <- alpha

    # Remove the cvFold variable
#    if (verbose)
#      rm(cvFold, pos = ".GlobalEnv")

    return(testAll)
    
 
  } # cvForAlpha
  
  # Now run over the alphas
  completeTest <- list2df(lapply(alphaVec, cvForAlpha))
  completeTest$seed <- seed

  # Now summarize the loss over the cv folds, with a loss value for each
  # alpha, lambda, and tau combination for a given seed
  # (recommend using the 'by' function)
  
  dfData <- list2df(dlply(completeTest,
                          .variables = c('alpha', 'lambda','tau','seed'),
                          .fun = function(x){
                          # x = K x K data.frame of values for the K folds with 
                          # same (alpha, lambda, tau, seed) parameter values.
                          Eloss <- sum(x$weightedSumLoss) / sum(x$sumWeights)
                          return(list('ExpectedLoss' = Eloss,
                                      'alpha' = unique(x$alpha), 
                                      'tau' = unique(x$tau),
                                      'lambda' = unique(x$lambda)))}))

  # Otherwise it looks gross
  rownames(dfData) <- NULL

  # The values of alpha, lambda, seed, and tau are the columns of dfData.
  # Using na.rm => no crashing if NAs present, but also need to warn that
  # NAs are being omitted
  paramsdf <- data.frame(dfData[,'alpha'], dfData[,'tau'], dfData[,'lambda'])
  colnames(paramsdf) <- c('alpha','tau','lambda')

  if(any(is.na(paramsdf)))
    warning('NA values are being omitted.')

  # These would be used to eventually calculate the mean and sd of
  # the parameters across multiple bootstrap replicates.  For now,
  # they only serve to provide the 'sd' for creating the final
  # lambdaVec for the fit.
  Eparams <- apply(paramsdf, 2, function(x) mean(x, na.rm = TRUE))
  sdParams <- apply(paramsdf, 2, function(x) sd(x, na.rm = TRUE))

  # Searching for the minimum by sorting
  # Smaller expected loss is preferred
  # In the event of a tie, smaller sqErrorTau is prferred (tau closer to 0.5)
  # If still tied, larger values of lambda are prefered because they reduce the
  # number of predictors to create a more parsimonous model with fewer predictors
  dfData$sqErrorTau <- (dfData$tau - 0.5)^2
  gridMinimum <- sort.data.frame(dfData, ~ExpectedLoss + sqErrorTau - lambda)[1:10,]

  # Construct a lambda sequence that is centered at the minimum lambda
  lambdaFinal <- sort(c(gridMinimum[1, "lambda"],
                        seq(gridMinimum[1, "lambda"] + sdParams['lambda'], 
                            max(gridMinimum[1, "lambda"] - sdParams['lambda'], 1e-04),
                            length = 15)),
                      decreasing = TRUE)

  # Fit to all the data
  glmnetFinal <- glmnet(predictors, truthLabels, weights = weight,
                        family = "binomial", lambda = lambdaFinal,
                        alpha = gridMinimum[1, "alpha"])

  # Add in the optimal parameters
  glmnetFinal$optimalParms <- gridMinimum

  # Assign the class
  class(glmnetFinal) <- c("glmnetGLR", class(glmnetFinal))

  
  return(glmnetFinal)
  

} # trainGLR 



##' @rdname trainGLR
##' @method print glmnetGLR
##' @param glmnetGLRobject a train elastic net logistic classifier
##' @param \dots Additional arguments to methods (ignored) 
##' @S3method print glmnetGLR

print.glmnetGLR <- function(glmnetGLRobject, ...) {

  cat("Top", NROW(glmnetGLRobject$optimalParms), "optimal parameter values for the elastic net logistic regression fit:\n\n")
  rownames(glmnetGLRobject$optimalParms) <- 1:NROW(glmnetGLRobject$optimalParms)
  print(glmnetGLRobject$optimalParms)

  invisible(NULL)
  
} # print.glmnetGLR



##' @rdname trainGLR
##' @method predict glmnetGLR
##' @param glmnetGLRobject a trained elastic net logistic regression classifier
##' @param newdata the new set of observations to be predicted. New data must be
##' an array with the same column names as the training data
##' @param \dots Additional arguments to 
##' predict method of \code{glmnet} objects are ignored, as we use the
##' optimal values of \code{lambda}, \code{tau}, and \code{alpha} that were
##' generated by \code{trainGLR()}.
##' @S3method predict glmnetGLR

predict.glmnetGLR <- function(glmnetGLRobject, newdata, truthCol = NULL,
                              keepCols = NULL) {

  # Switching from column numbers to column names if necessary
  if (!is.null(truthCol) & is.numeric(truthCol)) {
     truthCol <- colnames(newdata)[truthCol]
  }

  if (!is.null(keepCols) & is.numeric(keepCols)) {
     keepCols <- colnames(newdata)[keepCols]
  }


  # Verify the levels of truthCol match the class names in the glmnetGLRobject
  if (!is.null(truthCol)) {

    # It needs to be a factor
    newdata[,truthCol] <- as.factor(newdata[,truthCol])

    if (!setequal(levels(newdata[,truthCol]), glmnetGLRobject$classnames))
      warning("The class labels in the 'truthCol' do not match those ",
              "in the 'glmnetGLRobject'")

  }

  # Get the predictor names expected by the glmnetGLRobject
  predictorNames <- glmnetGLRobject$beta@Dimnames[[1]]

  # Make sure all the predictor names are in the newdata
  if (!all(predictorNames %in% colnames(newdata)))
   stop("The following predictors are expected by 'glmnetGLRobject' but are not\n",
        "present in 'newdata'\n'",
        paste(setdiff(predictorNames, colnames(newdata)), collapse = "', '"), "'\n")

  # Prepare newdata for prediction
  nd <- as.matrix(newdata[,predictorNames])
  
  if (!is.numeric(nd))
    stop("One or more of the predictor columns in 'newdata' is/are not numeric")
  
  # Get the original glmnet glmnetGLRobject
  glmnetObject <- glmnetGLRobject[-which(names(glmnetGLRobject) == "optimalParms")]
  class(glmnetObject) <- setdiff(class(glmnetGLRobject), "glmnetGLR")

  
  # Get the numeric (probability) predictions from predict.glmnet
  preds <- predict(glmnetObject, nd,
                   s = glmnetGLRobject$optimalParms[1, "lambda"], type = "response")

  # Dichotomize the prediction using the optimal tau    
  predLabels <- factor(preds > glmnetGLRobject$optimalParms[1,"tau"],
                       levels = c(FALSE, TRUE),
                       labels = glmnetObject$classnames)

  # If there were rownames in newdata, add them in
  if (!is.null(rn <- rownames(newdata)))
    names(predLabels) <- rn

  # Combine new data
  output <- cbind(predLabels, newdata[,truthCol], newdata[,keepCols])
  colnames(output) <- c("PredictClass", truthCol, keepCols)
  
  # Assign the class if a truth column was provided
  if (!is.null(truthCol)) {

    class(output) <- c("glmnetGLRpred", class(output))
    
    attributes(output) <- c(attributes(output),
                            list(truthCol = truthCol,
                                 optimalParms = glmnetGLRobject$optimalParms[1,],
                                 classNames = glmnetObject$classnames))


  }

  else {
    attributes(output) <- c(attributes(output),
                            list(optimalParms = glmnetGLRobject$optimalParms[1,],
                                 classNames = glmnetObject$classnames))
  }

  return(output)

} # predict.glmnetGLR


##' @rdname glmnetGLRpredSummary
##' @method summary glmnetGLRpred
##' @S3method summary glmnetGLRpred
##' @title Calculate a set of metrics to assess the performance of the LLRC
##' @description This function allows for a quick computation of multiple metrics 
##' to assess the performance of the elastic net classififer 
##' assuming the ground truth is available.
##' 
##' @param object the object produced by \code{predict}
##' @param ... additional arguments to other methods (currently ignored)
##' @return \code{data.frame} with sensitivity, specificity, false negative rate, false positive rate, and 
##' accuracy computed by class.
##' @author Alex Venzin
##' @note The ground truth for the new data set is required to generate these metrics.

summary.glmnetGLRpred <- function(object, ...) {

  truthCol <- attributes(object)$truthCol
  
  # If there are any missing data in the PredictClass or the Truth class,
  # remove them
  cc <- complete.cases(object[, c("PredictClass", truthCol)])

  if (any(!cc)) {
    
    warning("'object' has ", sum(!cc),
            " missing observation(s) which will be removed\n")

    object <- object[cc,]

  }
  
  # Calculate the confusion matrix  
  confusion_matrix <- confusion(object[, truthCol],
                                object[,"PredictClass"])
  
  # calculate sensitivity
  sens <- sensitivity(confusion_matrix)
  
  # calculate specificity
  spec <- specificity(confusion_matrix)
  
  # calculate accuracy
  acc <- accuracy(confusion_matrix)
  
  # false negative rate = 1 - TPR = 1 - sensitivity
  FNR <- function(confusionSummary, aggregate = c('micro', 'macro')){
    sens <- sensitivity(confusionSummary)
    return(1 - sens$byClass)
  }
  
  # false positive rate = 1 - TNR = 1 - specificity
  FPR <- function(confusionSummary, aggregate = c('micro', 'macro')){
    spec <- specificity(confusionSummary)
    return(1 - spec$byClass)
  }
  
  # calculate fnr
  fnr <- FNR(confusion_matrix)
  
  # calculate fpr
  fpr <- FPR(confusion_matrix)
  
  # put it into a data.frame format
  summ_out <- rbind(sens$byClass, spec$byClass, fnr, fpr, acc$byClass)

  # Only show the second class--the target of the prediction
  summ_out <- matrix(summ_out[,attributes(object)$classNames[2]], ncol = 1)

  colnames(summ_out) <- attributes(object)$classNames[2]
  rownames(summ_out) <- c("sensitivity", "specificity",
                          "false negative rate", "false positive rate",
                          "accuracy")
  
  return(as.data.frame(summ_out))
  
} # end summary.glmnetGLRpred
