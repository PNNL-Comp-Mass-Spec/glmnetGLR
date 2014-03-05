##' Train a LASSO logistic regression classifier
##'
##' This function uses lasso logistic regression and 
##' cross-validation to determine the optimal LASSO parameters 
##' that minimize the probability of committing classification
##' errors.  If ties exist in the optimal solution, then a tau closest to .5
##' is selected.
##'
##' 
##' 
##' @rdname trainingLLRC
##'
##' @param resp the binary response vector 
##'
##' @param predictors a matrix containing all the predictor variables on the columns
##'
##' @param Kappa the LASSO loss function parameter relating the severity of committing
##' a false negative error relative to a false positive error (i.e., if a false negative erorr
##' is twice as costly as a false positive error, then \code{Kappa} = 2) 
##'
##' @param lambdaVec a sequence of lambda values for which the lasso model is fit
##' 
##' @param tauVec a sequence of candiate tau threshold values for the logistic regression classifier
##'
##' @param cvFolds the number of cross validation folds
##'
##' @param naFilter is proportion of data for a given predictor (column) that isn't NA, if there is too many NAs then predictor is dropped
##'
##' @param seed set the random seed for sampling during cross validation
##'
##' @return 
##' \item{llrc}{The optimal regression model and associated parameters}
##' \item{output}{loss function values for all combinations of the parameters lambda and tau}
##'
##' @author Brett Amidan
##'
##' @examples  MODIFY as necessary
##' # Load the training dataset
##' data(traindata)
##'
##' # The numeric columns form my predictor variables
##' predictors <- traindata[,9:96]
##'
##' #   The regression variable Y is the binomial response of interest
##' resp <- traindata[,"BinomResp"]
##'
##' # We want to perform 5 fold cross validation and we wish to penalize 
##' # false negative errors with the same prejudice as false positive errors
##' 
##' lambdaVec=seq(0.05,0.15,by=0.02)
##'
##' llrc.out <- trainingLLRC(resp = resp,
##'                         predictors = predictors,
##'                         Kappa = 5,
##'                         lambdaVec = lambdaVec,
##'                         cvFolds = 5,
##'                         tauVec = seq(0.01, 0.99, by = 0.01),
##'                         NAfilter = .6,
##'                         seed = 1)
##'                         
##' plot(llrc.out,main="Test")


trainingLLRC <- function(resp, predictors, Kappa = 1, 
                         lambdaVec = seq(0.05, 0.15, by = 0.02),
                         tauVec = seq(0.01, 0.99, by = 0.01),
                         naFilter = .6,
                         cvFolds = 5, 
                         seed = 1) {

  ## Give warning if data isn't correct class
  if (class(predictors) != "data.frame") {
    stop("predictors is not a numeric data.frame")
  } 
  
  # If the response variable is a factor, convert it to 0/1 values
  if(class(resp) == 'factor'){
    
    binary.vars <- levels(resp)
    
    if(length(levels(resp)) != 2){
      stop(paste("The number of levels implies that", 
                 resp, "is not a binary response variable", sep = ""))
    }
    
    resp <- ifelse(resp == as.character(binary.vars[1]), 0, 1)
    
    print(paste("The value 0 has been assigned to", binary.vars[1], "and the value 1 has 
                been assigned to", binary.vars[2], ". If the opposite value designation is desired,
                then rearrange the levels of the response variable.", sep = ""))
    
  }
  
  ## calculate proportions of 0 and 1 from truth
  prop1 <- sum(resp)/length(resp)
  prop0 <- 1 - prop1

  ## create temp output
  output <- matrix(NA,nrow=0,ncol=4)
  colnames(output)<- c("lambda","kappa","tau","lossValue")

  ###  Loop thru each lambda
  for (lambda in lambdaVec) {
    
    ### Do the cross validation logistic regression using Lasso
    lr.out <- CVlogisticReg2(x=predictors,resp=resp,lambda=lambda,
      cvFolds=cvFolds, naFilter = naFilter, seed=seed)

    ### Find the optimal cut given kappa
    lossVals <- NULL
    
    ## loop thru each tau
    for (i in tauVec) {
      ## count # predicting poor, but actually good
      ind1 <- lr.out$pred >= i & lr.out$actual == 0
      ## count # predicting good, but actually poor * kappa
      ind2 <- lr.out$pred < i & lr.out$actual == 1
      ## score using the loss function
      lossVals <- c(lossVals,(sum(ind1)+sum(ind2)*Kappa)/
        length(lr.out$pred))
    }
    ### store the results
    temp1 <- cbind(rep(lambda,length(tauVec)),rep(Kappa,length(tauVec)),
      tauVec,lossVals)
    output <- rbind(output,temp1)
  }
  
  ### Find the minimum loss value
  indy2 <- output[,"lossValue"] == min(output[,"lossValue"])
  minvals <- output[indy2,]
  ### if ties, look for the closest to tau=.5
  if (sum(indy2)>1) {
    tt <- abs(minvals[,"tau"]-.5)
    indy3 <- tt == min(tt)
    minvals <- minvals[indy3,]
    ## if still ties pick lowest lambda
    if (sum(indy3)>1) {
      indy4 <- minvals[,"lambda"]==min(minvals[,"lambda"])
      minvals <- minvals[indy4,]
    } # ends if
  } # ends if
  
  ###### Perform the optimal LLRC
  tempdata <- cbind(resp,predictors)
  colnames(tempdata)[1] <- "Response"
  
  ## remove columns of data with too many NAs
  ind <- colSums(!is.na(tempdata))/nrow(tempdata) > naFilter
  tempdata <- tempdata[,ind]
  
  ## remove any rows that are NAs
  ind2 <- !is.na(rowSums(tempdata))
  tempdata <- tempdata[ind2,]
  
  ## fit the model
  llrc.out <- logisticReg(data = tempdata, respName="Response", lambda = minvals["lambda"])
  llrc.out <- c(llrc.out,minvals)
  
  #######  Prepare Output
  
  outlist <- list(llrc = llrc.out, output = output)
  
  class(outlist) <- unique(c("qcdm", class(outlist)))

  return(outlist)

} # ends function

################################################################################
### Plotting, prediction, and summary methods
################################################################################

##' @rdname trainingLLRC
##' @method plot qcdm
##' @param x the second element in the list of returned items from 
##' \code{trainingLLRC}
##' @param \dots additional arguments to pass to the plot function
##' @S3method plot qcdm

# Put x back in for this 
# Plot methods
plot.qcdm <- function(x, ...) {

  # Extract the data for the plots
  plotdata <- x$output

  # Get plotting parameters
  plotParms <- list(...)
  
  # Default plot parameters, to be used if not provided
  # If xlab not provided, then use the names of saFuns...
  if(!("xlab" %in% names(plotParms)))
    plotParms$xlab <- expression(tau)
  
  if(!("ylab" %in% names(plotParms)))
    plotParms$ylab <- "Loss Values"
  
  if(!("type" %in% names(plotParms)))
    plotParms$type <- "n"
  
  if(!('sub' %in% names(plotParms)))
    plotParms$sub <- paste("(Kappa = ",plotdata[1,"kappa"],")",sep="")
  
  ## Add in the data to be plotted
  plotParms$x <- plotdata[,"tau"]
  plotParms$y <- plotdata[,"lossValue"]

  # Make the plot
  do.call("plot", plotParms)

  ## get list of all lambda's used
  lambda.vec <- unique(plotdata[,"lambda"])
  
  ## add lines
  for (i in rev(1:length(lambda.vec))) {
    indy <- plotdata[,"lambda"] == lambda.vec[i]
    lines(plotdata[indy,"tau"], plotdata[indy,"lossValue"], col=i,lwd=2)

  }
  ## add legend
  legend(x="topright",legend=lambda.vec,col=1:length(lambda.vec),
    lwd=2,cex=.8)

  invisible()
} # ends function

##' @rdname qcdmPred
##' @method predict qcdm
##' @S3method predict qcdm
##' @title Predict the class of a new observation
##' @description Use the LLRC to predict the classes for a new set of observations
##' @export
##' 
##' @param object the regression model produced by \code{trainingLLRC}
##' @param newdata the new observations that need class predictions
##' @param truthCol the column name containing the true class values
##' @param keepCols the set of predictor variables that get passed along with the LLRC data
##' @param ... additional arguments 
##' 
##' @return 
##' \item{Results}{a \code{data.frame} with logistic function value, predicted class, 
##' truthCol, and keepCols}
##' 
##' @author Brett Amidan
##' 
##' @note Logistic regression in this package operates on two classes. 
##' If the user is interested in modeling more than two classes, then 
##' multinomial logistic regression is required (not available in this package).
##' Any columns named "Probability" or "PredictClass" will be overwritten
##' using the \code{predict} function.
##' 
##' @examples
##' 
##'# Load the training dataset
##' data(traindata)
##'
##' # The numeric columns form my predictor variables
##' predictors <- traindata[,9:96]
##'
##' #   The regression variable Y is the binomial response of interest
##' resp <- traindata[,"BinomResp"]
##'
##' # We want to perform 5 fold cross validation and we wish to penalize 
##' # false negative errors with the same prejudice as false positive errors
##' 
##' lambdaVec=seq(0.05,0.15,by=0.02)
##'
##' llrc.out <- trainingLLRC(resp = resp,
##'                         predictors = predictors,
##'                         Kappa = 1,
##'                         lambdaVec = lambdaVec,
##'                         cvFolds = 5,
##'                         tauVec = seq(0.01, 0.99, by = 0.01),
##'                         NAfilter = .6,
##'                         seed = 1)
##'                                                  
##' # We will use the trained LLRC to predict the class of new observations in testdata
##' data(testdata)
##'
##' prediction_values <- predict(object = llrc.out,
##'                              newdata = testdata,
##'                              truthCol = "Curated_Quality")

predict.qcdm <- function(object, newdata, truthCol = NULL, keepCols = NULL, ...) {
  
  stopifnot(is.data.frame(newdata))

  # Switching from column numbers to column names if necessary
  if (!is.null(truthCol) & is.numeric(truthCol)) {
     truthCol <- colnames(newdata)[truthCol]
  }

  if (!is.null(keepCols) & is.numeric(keepCols)) {
     keepCols <- colnames(newdata)[keepCols]
  }
  
  # If truth column provided, convert the two responses into 0/1 variables 
  if(!is.null(truthCol)){
    # Need to pull out the two variables and assign 0/1 designation 
  
    factor_tc <- factor(newdata[,truthCol])
  
    if(length(levels(factor_tc)) != 2){
      stop(paste("After conversion to factor, the number of levels implies that", 
                  truthCol, "is not a binary response variable", sep = ""))
    }
  
    classes <- levels(factor_tc)
  
    truth <- ifelse(newdata[,truthCol] == as.character(classes[1]), 0, 1)
    newdata[,truthCol] <- truth
  }
  
  # Drop the columns named "Probability" and "PredictClass" from newdata
  
  if("Probability" %in% colnames(newdata)){
    newdata <- subset(newdata, select = -Probability)
    warning("The column named Probability has been overwritten.")
  }
  
  if("PredictClass" %in% colnames(newdata)){
    newdata <- subset(newdata, select = -PredictClass)
    warning("The column named PredictClass has been overwritten.")
  }
  
  # keepCols are column names or numbers in the 'newdata' that will tag along
  # with the prediction
  
  # Extract the data for the predictions
  llrc <- object$llrc

  ## get the right data
  coefs <- names(llrc$coef)
  coefs <- coefs[2:length(coefs)]
  
  ind <- is.element(coefs,colnames(newdata))
  
  if (all(coefs %in% colnames(newdata))) {
    #### make the predictions from the newdata
    newdata2 <- cbind(rep(1,nrow(newdata)),newdata[,coefs])
    coef.mat <- matrix(llrc$coef,nrow=nrow(newdata2),ncol=ncol(newdata2),byrow=TRUE)
    scores <- rowSums(newdata2*coef.mat)
    ind <- scores > 100
    scores[ind] <- 100
    probs <- round(exp(scores)/(exp(scores)+1),2)
    ## predicted Class
    ind2 <- probs > llrc$tau
    pclass <- rep(0,length(probs))
    pclass[ind2] <- 1
    ## create output
    output <- cbind(probs,pclass)
    rownames(output) <- rownames(newdata2)
  } else {
    stop("newdata column names do not match variables in LLRC")
  }
  
  # [output : the columns you wish to keep] (concatenation)  
  
  output <- cbind(output, newdata[,truthCol], newdata[,keepCols])
  colnames(output) <- c("Probability", "PredictClass", truthCol, keepCols)
  
  # Assign the class
  class(output) <- c("qcdmPred", class(output))

  attributes(output) <- c(attributes(output),
                          truthCol = truthCol,
#                          lambdaHat = lambdaHat,
#                          tauHat = tauHat,
                          classNames = as.character(classes))
  
  # data.frame [probability : predict class : ...]
  return(output)
  
} # end of predict.qcdm

##' @rdname qcdmSummary
##' @method summary qcdmPred
##' @S3method summary qcdmPred
##' @title Calculate a set of metrics to assess the performance of the LLRC
##' @description This function allows for a quick computation of multiple metrics to assess the performance
##' of the LLRC given that the ground truth is available.
##' @export
##' @param object the object produced by \code{predict}
##' @param ... additional arguments
##' @return \code{data.frame} with sensitivity, specificity, false negative rate, false positive rate, and 
##' accuracy computed by class.
##' @author Alex Venzin
##' @note The ground truth for the new data set is required to generate these metrics.
##' @examples 
##' 
##' # Load the data produced during prediction
##' data(prediction_values)
##' 
##' summary(prediction_values)

# Summary methods

summary.qcdmPred <- function(object, ...) {
  
  classNames <- c(attributes(object)$classNames1, attributes(object)$classNames2)
  
  if (is.null(truthCol <- attributes(object)$truthCol))
    stop("The truth column was not indicated in the call to the ",
          "predict.qcdm method,\ni.e., 'truthCol = NULL'")
  
  confusion_matrix <- confusion(object[,truthCol],
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
  
  colnames(summ_out) <- classNames
  rownames(summ_out) <- c("sensitivity", "specificity",
                          "false negative rate", "false positive rate",
                          "accuracy")
  
  return(as.data.frame(summ_out))
  
} # end summary.qcdm


