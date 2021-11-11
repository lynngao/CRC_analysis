library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

known_thomas = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_thomas.csv",col_names = TRUE)
known_thomas = as.data.frame(known_thomas)
rownames(known_thomas) = known_thomas[,1]
known_thomas = known_thomas[,-1]
metadata_thomas <- read_csv("/home1/yilingao/proj/abundance/metadata_thomas.csv",col_names = F)
metadata_thomas$X1[metadata_thomas$X1 == "CTR"] = "control"
known_thomas = known_thomas[metadata_thomas$X3,]
known_thomas$status = metadata_thomas$X1

#renormalization
filtering_percent = 0.9
known_thomas = known_thomas[,colMeans(known_thomas[,1:ncol(known_thomas)]==0)<filtering_percent]
known_thomas[,1:(ncol(known_thomas)-1)] = known_thomas[,1:(ncol(known_thomas)-1)]/rowSums(known_thomas[,1:(ncol(known_thomas)-1)]) 

set.seed(101)
######################################################################################################
####################################################################################
auc_known = as.data.frame(matrix(NA,ncol = 1, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = "thomas"

for (i in 1:30) {
  sample_crc = sample.int(n=nrow(known_thomas[known_thomas$status == "CRC",]),size=floor(0.8*nrow(known_thomas[known_thomas$status == "CRC",])),replace=F)
  sample_control = sample.int(n=nrow(known_thomas[known_thomas$status == "control",]),size=floor(0.8*nrow(known_thomas[known_thomas$status == "control",])),replace=F)
  known_training_crc = known_thomas[known_thomas$status == "CRC",][sample_crc,]
  known_training_control = known_thomas[known_thomas$status == "control",][sample_control,]
  known_training = rbind(known_training_crc,known_training_control)
  known_test = known_thomas[rownames(known_thomas)%!in%rownames(known_training),]
  
  known_rf <- train(
    status ~ ., data = known_training, method = "rf", ntree = 1000, metric = "ROC",
    trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  )

  known_pred <- predict(known_rf, known_test, type="prob")
  auc_known[i,1] <- auc(known_test$status, known_pred[,1])
}

auc_known[31,1] <- mean(auc_known[1:30,1])
auc_known[32,1] <- sd(auc_known[1:30,1])

write.csv(auc_known, "/home1/yilingao/proj/abundance/auc_within_datasets/bracken/thomas_known_1000trees.csv")




