#library("pROC")
library("randomForest")
library("caret")

#'Train a classifier 
#'
#'This function trains a random forest classifier using different types of features
#'
#'@param FeatureMatrix a matrix of feature values plus the last column indicating positive or negative snps
#'@return a random forest classifier
#'@export
trainMod <- function(FeatureMatrix){
  x <- as.matrix(FeatureMatrix[,1:(dim(FeatureMatrix)[2] - 1)])
  y <- as.factor(as.character(FeatureMatrix[,dim(FeatureMatrix)[2]]))#last column indicates their labels
  train_control <- trainControl(method="cv", number=5,summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
  RFmodel <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(10,10))
  print("The random forest model has been trained!")
  return(RFmodel)
}

#'Compute prediction score 
#'
#'This function predicts the probablity score using a trained random forest classifier
#'
#'@param FeatureMatrix a matrix of feature values
#'@param RFmodel a random forest classifier
#'@return prediction score matrix
#'@export
predMod <- function(FeatureMatrix, RFmodel){
  xtest <- as.matrix(FeatureMatrix)
  pred <- predict(RFmodel, xtest, probability=TRUE)
  prob <- predict(RFmodel, xtest, type="prob")
  return(prob)
}

