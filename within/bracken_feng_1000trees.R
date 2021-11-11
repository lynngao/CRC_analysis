library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

known_feng = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_feng.csv",col_names = TRUE)
known_feng = as.data.frame(known_feng)
rownames(known_feng) = known_feng[,1]
known_feng = known_feng[,-1]
metadata_feng <- read_csv("/home1/yilingao/proj/abundance/metadata_feng.csv",col_names = T)
metadata_feng$X4[metadata_feng$X4 == "carcinoma"] = "CRC"
metadata_feng$X4[metadata_feng$X4 == "controls"] = "control"
known_feng = known_feng[metadata_feng$fastq,]
known_feng$status = metadata_feng$X4
##exclude adenoma samples
known_feng = subset(known_feng,known_feng$status!=c("advanced adenoma"))


#renormalization
filtering_percent = 0.9
known_feng = known_feng[,colMeans(known_feng[,1:ncol(known_feng)]==0)<filtering_percent] 
known_feng[,1:(ncol(known_feng)-1)] = known_feng[,1:(ncol(known_feng)-1)]/rowSums(known_feng[,1:(ncol(known_feng)-1)]) 

set.seed(101)
######################################################################################################
####################################################################################
auc_known = as.data.frame(matrix(NA,ncol = 1, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = "feng"

for (i in 1:30) {
  sample_crc = sample.int(n=nrow(known_feng[known_feng$status == "CRC",]),size=floor(0.8*nrow(known_feng[known_feng$status == "CRC",])),replace=F)
  sample_control = sample.int(n=nrow(known_feng[known_feng$status == "control",]),size=floor(0.8*nrow(known_feng[known_feng$status == "control",])),replace=F)
  known_training_crc = known_feng[known_feng$status == "CRC",][sample_crc,]
  known_training_control = known_feng[known_feng$status == "control",][sample_control,]
  known_training = rbind(known_training_crc,known_training_control)
  known_test = known_feng[rownames(known_feng)%!in%rownames(known_training),]
  
  known_rf <- train(
    status ~ ., data = known_training, method = "rf", ntree = 1000, metric = "ROC",
    trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  )
  known_pred <- predict(known_rf, known_test, type="prob")
  auc_known[i,1] <- auc(known_test$status, known_pred[,1])
}

auc_known[31,1] <- mean(auc_known[1:30,1])
auc_known[32,1] <- sd(auc_known[1:30,1])

write.csv(auc_known, "/home1/yilingao/proj/abundance/auc_within_datasets/bracken/feng_known_1000trees.csv")

