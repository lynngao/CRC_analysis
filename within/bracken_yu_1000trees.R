library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

known_yu = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_yu.csv",col_names = TRUE)
known_yu = as.data.frame(known_yu)
rownames(known_yu) = known_yu[,1]
known_yu = known_yu[,-1]
metadata_yu <- read_csv("/home1/yilingao/proj/abundance/metadata_yu.csv",col_names = FALSE)
metadata_yu$X16[metadata_yu$X16 == "CTR"] = "control"
known_yu = known_yu[metadata_yu$X1,]
known_yu$status = metadata_yu$X16

#renormalization
filtering_percent = 0.9
known_yu = known_yu[,colMeans(known_yu[,1:ncol(known_yu)]==0)<filtering_percent]
known_yu[,1:(ncol(known_yu)-1)] = known_yu[,1:(ncol(known_yu)-1)]/rowSums(known_yu[,1:(ncol(known_yu)-1)]) 

set.seed(101)
######################################################################################################
####################################################################################
auc_known = as.data.frame(matrix(NA,ncol = 1, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = "yu"

for (i in 1:30) {
  sample_crc = sample.int(n=nrow(known_yu[known_yu$status == "CRC",]),size=floor(0.8*nrow(known_yu[known_yu$status == "CRC",])),replace=F)
  sample_control = sample.int(n=nrow(known_yu[known_yu$status == "control",]),size=floor(0.8*nrow(known_yu[known_yu$status == "control",])),replace=F)
  known_training_crc = known_yu[known_yu$status == "CRC",][sample_crc,]
  known_training_control = known_yu[known_yu$status == "control",][sample_control,]
  known_training = rbind(known_training_crc,known_training_control)
  known_test = known_yu[rownames(known_yu)%!in%rownames(known_training),]
  
  known_rf <- train(
    status ~ ., data = known_training, method = "rf", ntree = 1000, metric = "ROC",
    trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  )

  known_pred <- predict(known_rf, known_test, type="prob")
  auc_known[i,1] <- auc(known_test$status, known_pred[,1])
}

auc_known[31,1] <- mean(auc_known[1:30,1])
auc_known[32,1] <- sd(auc_known[1:30,1])


write.csv(auc_known, "/home1/yilingao/proj/abundance/auc_within_datasets/bracken/yu_known_1000trees.csv")
