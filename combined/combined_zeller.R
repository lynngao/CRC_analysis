library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)

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
colnames(metadata_feng)[4] = "X4"
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

##################################combined auc#################################
z_test_vogtmann_combined <- readRDS("/home1/yilingao/proj/abundance/training:zeller/z_test_vogtmann_combined.rds")
colnames(z_test_vogtmann_combined)[(ncol(known_vogtmann)):(ncol(z_test_vogtmann_combined))] = paste0('u',colnames(z_test_vogtmann_combined)[(ncol(known_vogtmann)):(ncol(z_test_vogtmann_combined))])
z_test_vogtmann_combined$status = combined_vogtmann$status
z_test_feng_combined <- readRDS("/home1/yilingao/proj/abundance/training:zeller/z_test_feng_combined.rds")
colnames(z_test_feng_combined)[(ncol(known_feng)):(ncol(z_test_feng_combined))] = paste0('u',colnames(z_test_feng_combined)[(ncol(known_feng)):(ncol(z_test_feng_combined))])
z_test_feng_combined = z_test_feng_combined[rownames(z_test_feng_combined)%in%rownames(combined_feng),]
z_test_feng_combined$status = combined_feng$status
z_test_hannigan_combined <- readRDS("/home1/yilingao/proj/abundance/training:zeller/z_test_hannigan_combined.rds")
colnames(z_test_hannigan_combined)[(ncol(known_hannigan)):(ncol(z_test_hannigan_combined))] = paste0('u',colnames(z_test_hannigan_combined)[(ncol(known_hannigan)):(ncol(z_test_hannigan_combined))])
z_test_hannigan_combined = z_test_hannigan_combined[rownames(z_test_hannigan_combined)%in%rownames(combined_hannigan),]
z_test_hannigan_combined$status = combined_hannigan$status
z_test_thomas_combined <- readRDS("/home1/yilingao/proj/abundance/training:zeller/z_test_thomas_combined.rds")
colnames(z_test_thomas_combined)[(ncol(known_thomas)):(ncol(z_test_thomas_combined))] = paste0('u',colnames(z_test_thomas_combined)[(ncol(known_thomas)):(ncol(z_test_thomas_combined))])
z_test_thomas_combined$status = combined_thomas$status
z_test_yu_combined <- readRDS("/home1/yilingao/proj/abundance/training:zeller/z_test_yu_combined.rds")
colnames(z_test_yu_combined)[(ncol(known_yu)):(ncol(z_test_yu_combined))] = paste0('u',colnames(z_test_yu_combined)[(ncol(known_yu)):(ncol(z_test_yu_combined))])
z_test_yu_combined$status = combined_yu$status

auc_cross_cohort = as.data.frame(matrix(NA,ncol = 5, nrow = 32))
colnames(auc_cross_cohort) = c("training:zeller, test:yu", "training:zeller, test:hannigan", "training:zeller, test:feng","training:zeller, test:vogtmann","training:zeller, test:thomas")
rownames(auc_cross_cohort) = c(1:30,"mean","sd")

###################training:zeller###########################
combined_zeller[,1:(ncol(combined_zeller)-1)] = combined_zeller[,1:(ncol(combined_zeller)-1)]/rowSums(combined_zeller[,1:(ncol(combined_zeller)-1)]) 
training = combined_zeller

features <- function(training, test){
  test_status = test$status
  test = test[,-ncol(test)]
  for (i in 1:ncol(training)) {
  if(colnames(training)[i] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[i]
  }
}
test = test/rowSums(test)
test$status = test_status
return(test)
}

rf_model <- function(training,test1,test2,test3,test4,test5,auc_cross_cohort) {  
test1 = features(training,test1)
test2 = features(training,test2)
test3 = features(training,test3)
test4 = features(training,test4)
test5 = features(training,test5)

avg_pred1 = list(0)
avg_pred2 = list(0)
avg_pred3 = list(0)
avg_pred4 = list(0)
avg_pred5 = list(0)
for (i in 1:30){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 5000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred1 <- predict(rf, test1, type="prob")
pred2 <- predict(rf, test2, type="prob")
pred3 <- predict(rf, test3, type="prob")
pred4 <- predict(rf, test4, type="prob")
pred5 <- predict(rf, test5, type="prob")
avg_pred1 = avg_pred1 + pred1
avg_pred2 = avg_pred2 + pred2
avg_pred3 = avg_pred3 + pred3
avg_pred4 = avg_pred4 + pred4
avg_pred5 = avg_pred5 + pred5
auc_cross_cohort[i,1] <- auc(test1$status, pred1[,1])
auc_cross_cohort[i,2] <- auc(test2$status, pred2[,1])
auc_cross_cohort[i,3] <- auc(test3$status, pred3[,1])
auc_cross_cohort[i,4] <- auc(test4$status, pred4[,1])
auc_cross_cohort[i,5] <- auc(test5$status, pred5[,1])
}
avg_pred1 = avg_pred1/30
avg_pred2 = avg_pred2/30
avg_pred3 = avg_pred3/30
avg_pred4 = avg_pred4/30
avg_pred5 = avg_pred5/30
write.csv(avg_pred1, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/zeller/combined_pred_test_yu.csv")
write.csv(avg_pred2, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/zeller/combined_pred_test_hannigan.csv")
write.csv(avg_pred3, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/zeller/combined_pred_test_feng.csv")
write.csv(avg_pred4, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/zeller/combined_pred_test_vogtmann.csv")
write.csv(avg_pred5, "/home1/yilingao/proj/abundance/auc_cross_datasets/pred/zeller/combined_pred_test_thomas.csv")

return(auc_cross_cohort)
}


test1 = z_test_yu_combined
test2 = z_test_hannigan_combined
test3 = z_test_feng_combined
test4 = z_test_vogtmann_combined
test5 = z_test_thomas_combined

auc_cross_cohort <- rf_model(training,test1,test2,test3,test4,test5,auc_cross_cohort)

write.csv(auc_cross_cohort, "/home1/yilingao/proj/abundance/auc_cross_datasets/cross_datasets_combined_zeller_5000trees.csv")
