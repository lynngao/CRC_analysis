library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

combined_zeller = readRDS("/home1/yilingao/proj/abundance/combined_abundance_zeller.rds")
unknown_zeller = readRDS("/home1/yilingao/proj/abundance/unknown_abundance_zeller.rds")
unknown_zeller = as.data.frame(unknown_zeller)
known_zeller = readRDS("/home1/yilingao/proj/abundance/centrifuge_abundance_zeller.rds")
metadata_zeller = read_csv("/home1/yilingao/proj/abundance/Zeller_CRC.csv")
combined_zeller = combined_zeller[order(rownames(combined_zeller)),]
unknown_zeller = unknown_zeller[order(rownames(unknown_zeller)),]
known_zeller = known_zeller[order(rownames(known_zeller)),]
metadata_zeller = metadata_zeller[order(metadata_zeller$No),]
combined_zeller$status = metadata_zeller$study_condition
known_zeller$status = metadata_zeller$study_condition
colnames(combined_zeller)[(ncol(known_zeller)):(ncol(combined_zeller)-1)] = paste0('u',colnames(combined_zeller)[(ncol(known_zeller)):(ncol(combined_zeller)-1)])
colnames(unknown_zeller)[1:ncol(unknown_zeller)] = paste0('u',colnames(unknown_zeller)[1:ncol(unknown_zeller)])
unknown_zeller$status = metadata_zeller$study_condition

#renormalization
known_zeller[,1:(ncol(known_zeller)-1)] = known_zeller[,1:(ncol(known_zeller)-1)]/rowSums(known_zeller[,1:(ncol(known_zeller)-1)]) 

set.seed(101)
######################################################################################################
####################################################################################
auc_known = as.data.frame(matrix(NA,ncol = 1, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = "zeller"

for (i in 1:30) {
  sample_crc = sample.int(n=nrow(known_zeller[known_zeller$status == "CRC",]),size=floor(0.8*nrow(known_zeller[known_zeller$status == "CRC",])),replace=F)
  sample_control = sample.int(n=nrow(known_zeller[known_zeller$status == "control",]),size=floor(0.8*nrow(known_zeller[known_zeller$status == "control",])),replace=F)
  known_training_crc = known_zeller[known_zeller$status == "CRC",][sample_crc,]
  known_training_control = known_zeller[known_zeller$status == "control",][sample_control,]
  known_training = rbind(known_training_crc,known_training_control)
  known_test = known_zeller[rownames(known_zeller)%!in%rownames(known_training),]
 
  known_rf <- train(
    status ~ ., data = known_training, method = "rf", ntree = 1000, metric = "ROC",
    trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  )

  known_pred <- predict(known_rf, known_test, type="prob")
  auc_known[i,1] <- auc(known_test$status, known_pred[,1])
}

auc_known[31,1] <- mean(auc_known[1:30,1])
auc_known[32,1] <- sd(auc_known[1:30,1])

write.csv(auc_known, "/home1/yilingao/proj/abundance/auc_within_datasets/zeller_known_1000trees.csv")

