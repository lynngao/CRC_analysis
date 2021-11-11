library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

known_hannigan = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_hannigan.csv",col_names = TRUE)
known_hannigan = as.data.frame(known_hannigan)
rownames(known_hannigan) = known_hannigan[,1]
known_hannigan = known_hannigan[,-1]
metadata_hannigan <- read_csv("/home1/yilingao/proj/abundance/metadata_hannigan.csv",col_names = FALSE)
metadata_hannigan$X5[metadata_hannigan$X5 == "Cancer"] = "CRC"
metadata_hannigan$X5[metadata_hannigan$X5 == "Healthy"] = "control"
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665060",]
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665121",]
known_hannigan = known_hannigan[metadata_hannigan$X3,]
known_hannigan$status = metadata_hannigan$X5
##exclude adenoma samples
known_hannigan = subset(known_hannigan,known_hannigan$status!=c("Adenoma"))

#renormalization
filtering_percent = 0.9
known_hannigan = known_hannigan[,colMeans(known_hannigan[,1:ncol(known_hannigan)]==0)<filtering_percent]
known_hannigan[,1:(ncol(known_hannigan)-1)] = known_hannigan[,1:(ncol(known_hannigan)-1)]/rowSums(known_hannigan[,1:(ncol(known_hannigan)-1)]) 

set.seed(101)
######################################################################################################
####################################################################################
auc_known = as.data.frame(matrix(NA,ncol = 1, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = "hannigan"

for (i in 1:30) {
  sample_crc = sample.int(n=nrow(known_hannigan[known_hannigan$status == "CRC",]),size=floor(0.8*nrow(known_hannigan[known_hannigan$status == "CRC",])),replace=F)
  sample_control = sample.int(n=nrow(known_hannigan[known_hannigan$status == "control",]),size=floor(0.8*nrow(known_hannigan[known_hannigan$status == "control",])),replace=F)
  known_training_crc = known_hannigan[known_hannigan$status == "CRC",][sample_crc,]
  known_training_control = known_hannigan[known_hannigan$status == "control",][sample_control,]
  known_training = rbind(known_training_crc,known_training_control)
  known_test = known_hannigan[rownames(known_hannigan)%!in%rownames(known_training),]
  
  known_rf <- train(
    status ~ ., data = known_training, method = "rf", ntree = 1000, metric = "ROC",
    trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  )

  
  known_pred <- predict(known_rf, known_test, type="prob")
  auc_known[i,1] <- auc(known_test$status, known_pred[,1])
}

auc_known[31,1] <- mean(auc_known[1:30,1])
auc_known[32,1] <- sd(auc_known[1:30,1])

write.csv(auc_known, "/home1/yilingao/proj/abundance/auc_within_datasets/bracken/hannigan_known_1000trees.csv")


