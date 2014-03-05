### R code from vignette source 'glmnetGLR.Rnw'

###################################################
### code chunk number 1: install
###################################################
require(glmnetGLR)


###################################################
### code chunk number 2: data
###################################################
# Load the VOrbitrap Shewanella QC data
data(traindata)

# Here are the first few observations of the datasets
traindata[1:5, 1:7]

# Here we select the predictor variables
predictors <- as.matrix(traindata[,9:96])

# The logistic regression model requires a binary response
# variable. 
resp <- traindata[,"response"]

# Set the loss matrix
lM <- lossMatrix(c("Good","Good","Poor","Poor"),
                 c("Good","Poor","Good","Poor"),
                 c(     0,     1,     5,     0))

# Train the elastic net classifier
elasticNet <- trainGLR(truthLabels = resp,
                      predictors = predictors,
                      lossMat = lM)


###################################################
### code chunk number 3: checkdat
###################################################
elasticNet


###################################################
### code chunk number 4: predict
###################################################
# load the data for testing
data(testdata)

prediction_values <- predict(glmnetGLRobject = elasticNet, 
                             newdata = testdata,
                             truthCol = 'Curated_Quality',
                             keepCols = 12:14)


###################################################
### code chunk number 5: summary
###################################################
summary(prediction_values)


