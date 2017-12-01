#library("pROC")
#library("randomForest")
#library("caret")

#'Train a classifier 
#'
#'This function ...
#'
#'@param FeatureMatrix a matrix of feature values
#'@param pos a vector specifying the positive set of elements
#'@param neg a vector specifying the control set of elements
#'@return a random forest classifier
trainMod <- function(FeatureMatrix){
  x <- as.matrix(FeatureMatrix[,1:(dim(FeatureMatrix)[2] - 1)])
  y <- as.factor(as.character(FeatureMatrix[,dim(FeatureMatrix)[2]]))#last column indicates their labels
  train_control <- trainControl(method="cv", number=5,summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
  RFmodel <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
  print("The random forest model has been trained!")
  return(RFmodel)
}

#'Compute prediction score 
#'
#'This function ...
#'
#'@param FeatureMatrix a matrix of feature values
#'@param RFmodel a random forest classifier
#'@return a random forest classifier
predMod <- function(FeatureMatrix, RFmodel){
  return(predScore)
}


# #############################################
# Features_all <- read.table("Combined_Features_avg.txt", sep = "\t", head=T)
# Features_all[is.na(Features_all)] <- 0
# #Features_all <- Features_all[,-c(4)]
# #Features_all <- Features_all[, c(1:7, dim(Features_all)[2]-1, dim(Features_all)[2])]#network feature alone
# #Features_all <- Features_all[, c(1,8:181,dim(Features_all)[2]-1, dim(Features_all)[2])]#gwawa
# #Features_all <- Features_all[, c(1, 182:dim(Features_all)[2])]#funseq
# 
# snp_dis <- subset(Features_all, Features_all[,dim(Features_all)[2]] == "positive")
# snp_rnd <- subset(Features_all, Features_all[,dim(Features_all)[2]-1] == "RND")
# snp_tss <- subset(Features_all, Features_all[,dim(Features_all)[2]-1] == "TSS")
# snp_1k <- subset(Features_all, Features_all[,dim(Features_all)[2]-1] == "1K")
# snp_1k[, dim(Features_all)[2]-1] <- "X1K"
# 
# name_dis <- snp_dis[,1]
# name_rnd <- snp_rnd[,1]
# name_tss <- snp_tss[,1]
# name_1k <- snp_1k[,1]
# snp_dis <- snp_dis[,-c(1,dim(snp_dis)[2] - 1)]
# snp_rnd <- snp_rnd[,-c(1,dim(snp_rnd)[2])]
# snp_tss <- snp_tss[,-c(1,dim(snp_tss)[2])]
# snp_1k <- snp_1k[,-c(1, dim(snp_1k)[2])]
# colnames(snp_dis)[dim(snp_dis)[2]] <- "category"
# colnames(snp_rnd)[dim(snp_rnd)[2]] <- "category"
# colnames(snp_tss)[dim(snp_tss)[2]] <- "category"
# colnames(snp_1k)[dim(snp_1k)[2]] <- "category"
# 
# com_rnd <- rbind(snp_dis, snp_rnd)
# com_tss <- rbind(snp_dis, snp_tss)
# com_1k <- rbind(snp_dis, snp_1k)
# 
# x <- as.matrix(com_rnd[,1:(dim(com_rnd)[2] - 1)])
# y <- as.factor(as.character(com_rnd[,dim(com_rnd)[2]]))
# train_control <- trainControl(method="cv", number=5,summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
# rfFit_Com_RND <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# rfFit <- rfFit_Com_RND
# #rfFit_Net_RND <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(100,100))
# #rfFit <- rfFit_Net_RND
# #rfFit_gwawa_RND <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_gwawa_RND
# #rfFit_funseq_RND <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_funseq_RND
# print(rfFit)
# selectedIndices <- rfFit$pred$mtry == 6
# ROC_Com_RND <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_Net_RND <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_gwawa_RND <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_funseq_RND <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# plot(ROC)
# 
# 
# x <- as.matrix(com_tss[,1:(dim(com_tss)[2] - 1)])
# y <- as.factor(as.character(com_tss[,dim(com_tss)[2]]))
# train_control <- trainControl(method="cv", number=10,summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
# #train_control <- trainControl(method="repeatedcv", number=10, repeats=3,summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
# rfFit_Com_TSS <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# rfFit <- rfFit_Com_TSS
# #rfFit_Net_TSS <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_Net_TSS
# #rfFit_gwawa_TSS <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_gwawa_TSS
# #rfFit_funseq_TSS <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_funseq_TSS
# print(rfFit)
# selectedIndices <- rfFit$pred$mtry == 6
# ROC_Com_TSS <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_Net_TSS <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_gwawa_TSS <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_funseq_TSS <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# plot(ROC_funseq_TSS)
# 
# 
# x <- as.matrix(com_1k[,1:(dim(com_1k)[2] - 1)])
# y <- as.factor(as.character(com_1k[,dim(com_1k)[2]]))
# train_control <- trainControl(method="cv", number=10,summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
# #train_control <- trainControl(method="repeatedcv", number=10, repeats=3,summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
# rfFit_Com_1k <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# rfFit <- rfFit_Com_1k
# #rfFit_Net_1k <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_Net_1k
# #rfFit_gwawa_1k <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_gwawa_1k
# #rfFit_funseq_1k <- train(x, y, trControl=train_control, method="rf", strata=y, sampsize=c(50,50))
# #rfFit <- rfFit_funseq_1k
# print(rfFit)
# selectedIndices <- rfFit$pred$mtry == 186
# ROC_Com_1k <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_Net_1k <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_gwawa_1k <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# #ROC_funseq_1k <- roc(rfFit$pred$obs[selectedIndices], rfFit$pred$positive[selectedIndices])
# plot(ROC)
# 
# 
# 
# #######################################################################################
# load("Aug22_OPT_classification.RData")
# print(rfFit_Com_RND_opt)
# selectedIndices <- rfFit_Com_RND_opt$pred$mtry == 40
# ROC_Com_RND <- roc(rfFit_Com_RND_opt$pred$obs[selectedIndices], rfFit_Com_RND_opt$pred$positive[selectedIndices])
# 
# print(rfFit_Com_TSS_opt)
# selectedIndices <- rfFit_Com_TSS_opt$pred$mtry == 18
# ROC_Com_TSS <- roc(rfFit_Com_TSS_opt$pred$obs[selectedIndices], rfFit_Com_TSS_opt$pred$positive[selectedIndices])
# 
# 
# postscript("../Figures/ROC_cv10_RF_RND.eps")
# plot(smooth(ROC_Com_RND),col="red")
# lines(smooth(ROC_Net_RND), col="green")
# lines(smooth(ROC_gwawa_RND), col="orange")
# lines(smooth(ROC_funseq_RND), col="blue")
# dev.off()
# 
# postscript("../Figures/ROC_cv10_RF_TSS.eps")
# plot(smooth(ROC_Com_TSS),col="red")
# lines(smooth(ROC_Net_TSS), col="green")
# lines(smooth(ROC_gwawa_TSS), col="orange")
# lines(smooth(ROC_funseq_TSS), col="blue")
# dev.off()
# 
# postscript("../Figures/ROC_cv10_RF_1k.eps")
# plot(smooth(ROC_Com_1k),col="red")
# lines(smooth(ROC_Net_1k), col="green")
# lines(smooth(ROC_gwawa_1k), col="orange")
# lines(smooth(ROC_funseq_1k), col="blue")
# dev.off()
# 
# smooth(ROC_Com_RND)
# smooth(ROC_Net_RND)
# smooth(ROC_gwawa_RND)
# smooth(ROC_funseq_RND)
# CI_Com_RND <- ci.auc(smooth(ROC_Com_RND))
# CI_Net_RND <- ci.auc(smooth(ROC_Net_RND))
# CI_gwawa_RND <- ci.auc(smooth(ROC_gwawa_RND))
# CI_funseq_RND <- ci.auc(smooth(ROC_funseq_RND))
# 
# CI_Com_TSS <- ci.auc(smooth(ROC_Com_TSS))
# CI_Net_TSS <- ci.auc(smooth(ROC_Net_TSS))
# CI_gwawa_TSS <- ci.auc(smooth(ROC_gwawa_TSS))
# CI_funseq_TSS <- ci.auc(smooth(ROC_funseq_TSS))
# 
# CI_Com_1k <- ci.auc(smooth(ROC_Com_1k))
# CI_Net_1k <- ci.auc(smooth(ROC_Net_1k))
# CI_gwawa_1k <- ci.auc(smooth(ROC_gwawa_1k))
# CI_funseq_1k <- ci.auc(smooth(ROC_funseq_1k))
# 
# Com_auc_set <- c(smooth(ROC_Com_RND)$auc, smooth(ROC_Com_TSS)$auc, smooth(ROC_Com_1k)$auc)
# Net_auc_set <- c(smooth(ROC_Net_RND)$auc, smooth(ROC_Net_TSS)$auc, smooth(ROC_Net_1k)$auc)
# gwawa_auc_set <- c(smooth(ROC_gwawa_RND)$auc, smooth(ROC_gwawa_TSS)$auc, smooth(ROC_gwawa_1k)$auc)
# funseq_auc_set <- c(smooth(ROC_funseq_RND)$auc, smooth(ROC_funseq_TSS)$auc, smooth(ROC_funseq_1k)$auc)
# 
# Com_up_set <- c(max(CI_Com_RND),max(CI_Com_TSS),max(CI_Com_1k))
# Net_up_set <- c(max(CI_Net_RND),max(CI_Net_TSS),max(CI_Net_1k))
# gwawa_up_set <- c(max(CI_gwawa_RND),max(CI_gwawa_TSS),max(CI_gwawa_1k))
# funseq_up_set <- c(max(CI_funseq_RND),max(CI_funseq_TSS),max(CI_funseq_1k))
# 
# com_mean <- cbind(funseq_auc_set, gwawa_auc_set, Net_auc_set, Com_auc_set)
# com_se <- cbind(funseq_up_set, gwawa_up_set, Net_up_set, Com_up_set)
# postscript("../Figures/AUC_cv10_RF.eps")
# barCenters <- barplot(com_mean, ylim=c(0,1), xaxt="n", beside=TRUE)
# labCenters <- barCenters[2,]
# axis(1, at=labCenters,labels=c("FunSeq", "GWAWA", "Network","Combined"))
# arrows(barCenters, com_mean, barCenters,
#        com_se, lwd = 1.5, angle = 90,
#        code = 2, length = 0.1)
# dev.off()
# 
# #######################
# #feature rank
# feature_rank_rnd <- varImp(rfFit_Com_RND)$importance$Overall
# feature_rank_tss <- varImp(rfFit_Com_TSS)$importance$Overall
# feature_rank_rnd <- varImp(rfFit_Com_1k)$importance$Overall
# temp_rank_rnd <- rank(-feature_rank_rnd)
# temp_rank_tss <- rank(-feature_rank_tss)
# temp_rank_1k <- rank(-feature_rank_1k)
# com_rank <- (temp_rank_rnd + temp_rank_tss + temp_rank_1k)/3
# names(com_rank) <- row.names(feature_rank_rnd)
# com_rank_frame <- cbind(temp_rank_rnd, temp_rank_tss,temp_rank_1k, com_rank)
# row.names(com_rank_frame) <- row.names(varImp(rfFit_Com_RND)$importance)
# write.table(com_rank_frame,"All_Features_importance_rank.txt",quote=F, sep="\t")
# 
# save(list=ls(),file="Aug22_CV_classification.RData")
# 
# 
# 
# ###########################################################
# #test different ROC curves
# roc.test(smooth(ROC_Com_RND), smooth(ROC_Net_RND))#8.652e-07
# roc.test(smooth(ROC_Com_RND), smooth(ROC_gwawa_RND))#0.01687
# roc.test(smooth(ROC_Com_RND), smooth(ROC_funseq_RND))#8.603e-10
# 
# roc.test(smooth(ROC_Com_TSS), smooth(ROC_Net_TSS))#0.4463
# roc.test(smooth(ROC_Com_TSS), smooth(ROC_gwawa_TSS))#2.985e-11
# roc.test(smooth(ROC_Com_TSS), smooth(ROC_funseq_TSS))#2.2e-16
# 
# roc.test(smooth(ROC_Com_1k), smooth(ROC_Net_1k))#0.9693
# roc.test(smooth(ROC_Com_1k), smooth(ROC_gwawa_1k))#2.625e-13
# roc.test(smooth(ROC_Com_1k), smooth(ROC_funseq_1k))#2.2e-16
