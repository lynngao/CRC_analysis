library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

known_yu = readRDS("/home1/yilingao/proj/abundance/centrifuge_abundance_yu.rds")
unknown_yu = readRDS("/home1/yilingao/proj/abundance/unknown_abundance_yu.rds")
combined_yu = readRDS("/home1/yilingao/proj/abundance/combined_abundance_yu.rds")
metadata_yu <- read_csv("/home1/yilingao/proj/abundance/metadata_yu.csv",col_names = FALSE)
metadata_yu$X16[metadata_yu$X16 == "CTR"] = "control"
colnames(combined_yu)[(ncol(known_yu)+1):(ncol(combined_yu))] = paste0('u',colnames(combined_yu)[(ncol(known_yu)+1):(ncol(combined_yu))])
unknown_yu = as.data.frame(unknown_yu)
colnames(unknown_yu) = paste0('u',colnames(unknown_yu))
known_yu$status = metadata_yu$X16
combined_yu$status = metadata_yu$X16
unknown_yu$status = metadata_yu$X16

known_hannigan = readRDS("/home1/yilingao/proj/abundance/centrifuge_abundance_hannigan.rds")
unknown_hannigan = readRDS("/home1/yilingao/proj/abundance/unknown_abundance_hannigan.rds")
unknown_hannigan = as.data.frame(unknown_hannigan)
combined_hannigan = readRDS("/home1/yilingao/proj/abundance/combined_abundance_hannigan.rds")
metadata_hannigan <- read_csv("/home1/yilingao/proj/abundance/metadata_hannigan.csv",col_names = FALSE)
metadata_hannigan$X5[metadata_hannigan$X5 == "Cancer"] = "CRC"
metadata_hannigan$X5[metadata_hannigan$X5 == "Healthy"] = "control"
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665060",]
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665121",]
colnames(combined_hannigan)[(ncol(known_hannigan)+1):(ncol(combined_hannigan))] = paste0('u',colnames(combined_hannigan)[(ncol(known_hannigan)+1):(ncol(combined_hannigan))])
colnames(unknown_hannigan) = paste0('u',colnames(unknown_hannigan))
known_hannigan$status = metadata_hannigan$X5
combined_hannigan$status = metadata_hannigan$X5
unknown_hannigan$status = metadata_hannigan$X5
##exclude adenoma samples
known_hannigan = subset(known_hannigan,known_hannigan$status!=c("Adenoma"))
combined_hannigan = subset(combined_hannigan,combined_hannigan$status!=c("Adenoma"))
unknown_hannigan = subset(unknown_hannigan,unknown_hannigan$status!=c("Adenoma"))


known_feng = readRDS("/home1/yilingao/proj/abundance/centrifuge_abundance_feng.rds")
unknown_feng = readRDS("/home1/yilingao/proj/abundance/unknown_abundance_feng.rds")
unknown_feng = as.data.frame(unknown_feng)
combined_feng = readRDS("/home1/yilingao/proj/abundance/combined_abundance_feng.rds")
metadata_feng <- read_csv("/home1/yilingao/proj/abundance/metadata_feng.csv",col_names = T)
metadata_feng$X4[metadata_feng$X4 == "carcinoma"] = "CRC"
metadata_feng$X4[metadata_feng$X4 == "controls"] = "control"
colnames(combined_feng)[(ncol(known_feng)+1):(ncol(combined_feng))] = paste0('u',colnames(combined_feng)[(ncol(known_feng)+1):(ncol(combined_feng))])
colnames(unknown_feng) = paste0('u',colnames(unknown_feng))
known_feng$status = metadata_feng$X4
unknown_feng$status = metadata_feng$X4
combined_feng$status = metadata_feng$X4
##exclude adenoma samples
known_feng = subset(known_feng,known_feng$status!=c("advanced adenoma"))
unknown_feng = subset(unknown_feng,unknown_feng$status!=c("advanced adenoma"))
combined_feng = subset(combined_feng,combined_feng$status!=c("advanced adenoma"))

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

known_thomas = readRDS("/home1/yilingao/proj/abundance/centrifuge_abundance_thomas.rds")
unknown_thomas = readRDS("/home1/yilingao/proj/abundance/unknown_abundance_thomas.rds")
unknown_thomas = as.data.frame(unknown_thomas)
combined_thomas = readRDS("/home1/yilingao/proj/abundance/combined_abundance_thomas.rds")
metadata_thomas <- read_csv("/home1/yilingao/proj/abundance/metadata_thomas.csv",col_names = F)
metadata_thomas$X1[metadata_thomas$X1 == "CTR"] = "control"
colnames(combined_thomas)[(ncol(known_thomas)+1):(ncol(combined_thomas))] = paste0('u',colnames(combined_thomas)[(ncol(known_thomas)+1):(ncol(combined_thomas))])
colnames(unknown_thomas) = paste0('u',colnames(unknown_thomas))
known_thomas$status = metadata_thomas$X1
unknown_thomas$status = metadata_thomas$X1
combined_thomas$status = metadata_thomas$X1

set.seed(101)

#############################cross-cohort known auc#########################################################
auc_cross_cohort = as.data.frame(matrix(NA,ncol = 5, nrow = 32))
colnames(auc_cross_cohort) = c("training:thomas, test:yu", "training:thomas, test:hannigan", "training:thomas, test:feng","training:thomas, test:vogtmann","training:thomas, test:zeller")
rownames(auc_cross_cohort) = c(1:30,"mean","sd")
###################training:thomas###########################
training = known_thomas
training[,1:(ncol(training)-1)] = training[,1:(ncol(training)-1)]/rowSums(training[,1:(ncol(training)-1)]) #renormalize

##1.test:yu
test = known_yu
for (i in 1:ncol(training)) {
  if(colnames(training)[i] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[i]
  }
}
test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
avg_pred = list(0)
for (i in 1:30){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 1000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred <- predict(rf, test, type="prob")
avg_pred = avg_pred + pred
auc_cross_cohort[i,1] <- auc(test$status, pred[,1])
}
avg_pred = avg_pred/30
#write.csv(avg_pred, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/thomas/known_pred_test_yu.csv")

##2.test:hannigan
test = known_hannigan
for (i in 1:ncol(training)) {
  if(colnames(training)[i] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[i]
  }
}
test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
avg_pred = list(0)
for (i in 1:30){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 1000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred <- predict(rf, test, type="prob")
avg_pred = avg_pred + pred
auc_cross_cohort[i,2] <- auc(test$status, pred[,1])
}
avg_pred = avg_pred/30
#write.csv(avg_pred, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/thomas/known_pred_test_hannigan.csv")

##3.test:feng
test = known_feng
for (i in 1:ncol(training)) {
  if(colnames(training)[i] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[i]
  }
}
test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
avg_pred = list(0)
for (i in 1:30){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 1000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred <- predict(rf, test, type="prob")
avg_pred = avg_pred + pred
auc_cross_cohort[i,3] <- auc(test$status, pred[,1])
}
avg_pred = avg_pred/30
#write.csv(avg_pred, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/thomas/known_pred_test_feng.csv")

##4.test:vogtmann
test = known_vogtmann
for (i in 1:ncol(training)) {
  if(colnames(training)[i] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[i]
  }
}
test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
avg_pred = list(0)
for (i in 1:30){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 1000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred <- predict(rf, test, type="prob")
avg_pred = avg_pred + pred
auc_cross_cohort[i,4] <- auc(test$status, pred[,1])
}
avg_pred = avg_pred/30
#write.csv(avg_pred, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/thomas/known_pred_test_vogtmann.csv")

##5.test:zeller
test = known_zeller
for (i in 1:ncol(training)) {
  if(colnames(training)[i] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[i]
  }
}
test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
avg_pred = list(0)
for (i in 1:30){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 1000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred <- predict(rf, test, type="prob")
avg_pred = avg_pred + pred
auc_cross_cohort[i,5] <- auc(test$status, pred[,1])
}
avg_pred = avg_pred/30
#write.csv(avg_pred, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/thomas/known_pred_test_zeller.csv")

auc_cross_cohort[31,1] <- mean(auc_cross_cohort[1:30,1])
auc_cross_cohort[32,1] <- sd(auc_cross_cohort[1:30,1])
auc_cross_cohort[31,2] <- mean(auc_cross_cohort[1:30,2])
auc_cross_cohort[32,2] <- sd(auc_cross_cohort[1:30,2])
auc_cross_cohort[31,3] <- mean(auc_cross_cohort[1:30,3])
auc_cross_cohort[32,3] <- sd(auc_cross_cohort[1:30,3])
auc_cross_cohort[31,4] <- mean(auc_cross_cohort[1:30,4])
auc_cross_cohort[32,4] <- sd(auc_cross_cohort[1:30,4])
auc_cross_cohort[31,5] <- mean(auc_cross_cohort[1:30,5])
auc_cross_cohort[32,5] <- sd(auc_cross_cohort[1:30,5])

write.csv(auc_cross_cohort, "/home1/yilingao/proj/abundance/auc_cross_datasets/1000trees/cross_datasets_known_thomas_1000trees.csv")

