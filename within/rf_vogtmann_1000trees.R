library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

known_vogtmann = readRDS("/home1/yilingao/proj/abundance/centrifuge_abundance_vogtmann.rds")
unknown_vogtmann = readRDS("/home1/yilingao/proj/abundance/unknown_abundance_vogtmann.rds")
unknown_vogtmann = as.data.frame(unknown_vogtmann)
combined_vogtmann = readRDS("/home1/yilingao/proj/abundance/combined_abundance_vogtmann.rds")
metadata_vogtmann <- read_csv("/home1/yilingao/proj/abundance/metadata_vogtmann.csv",col_names = T)
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 0] = "control"
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 1] = "CRC"
colnames(combined_vogtmann)[(ncol(known_vogtmann)+1):(ncol(combined_vogtmann))] = paste0('u',colnames(combined_vogtmann)[(ncol(known_vogtmann)+1):(ncol(combined_vogtmann))])
colnames(unknown_vogtmann) = paste0('u',colnames(unknown_vogtmann))
known_vogtmann$status = metadata_vogtmann$casectl
unknown_vogtmann$status = metadata_vogtmann$casectl
combined_vogtmann$status = metadata_vogtmann$casectl

#renormalization
known_vogtmann[,1:(ncol(known_vogtmann)-1)] = known_vogtmann[,1:(ncol(known_vogtmann)-1)]/rowSums(known_vogtmann[,1:(ncol(known_vogtmann)-1)]) 

set.seed(101)
######################################################################################################
####################################################################################
auc_known = as.data.frame(matrix(NA,ncol = 1, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = "vogtmann"

for (i in 1:30) {
  sample_crc = sample.int(n=nrow(known_vogtmann[known_vogtmann$status == "CRC",]),size=floor(0.8*nrow(known_vogtmann[known_vogtmann$status == "CRC",])),replace=F)
  sample_control = sample.int(n=nrow(known_vogtmann[known_vogtmann$status == "control",]),size=floor(0.8*nrow(known_vogtmann[known_vogtmann$status == "control",])),replace=F)
  known_training_crc = known_vogtmann[known_vogtmann$status == "CRC",][sample_crc,]
  known_training_control = known_vogtmann[known_vogtmann$status == "control",][sample_control,]
  known_training = rbind(known_training_crc,known_training_control)
  known_test = known_vogtmann[rownames(known_vogtmann)%!in%rownames(known_training),]
  
  known_rf <- train(
    status ~ ., data = known_training, method = "rf", ntree = 1000, metric = "ROC",
    trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  )

  known_pred <- predict(known_rf, known_test, type="prob")
  auc_known[i,1] <- auc(known_test$status, known_pred[,1])
}

auc_known[31,1] <- mean(auc_known[1:30,1])
auc_known[32,1] <- sd(auc_known[1:30,1])

write.csv(auc_known, "/home1/yilingao/proj/abundance/auc_within_datasets/vogtmann_known_1000trees.csv")


