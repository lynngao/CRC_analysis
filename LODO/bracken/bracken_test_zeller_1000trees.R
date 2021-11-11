library(readr)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)

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
known_hannigan = known_hannigan[,colMeans(known_hannigan[,1:ncol(known_hannigan)]==0)<filtering_percent]
known_hannigan[,1:(ncol(known_hannigan)-1)] = known_hannigan[,1:(ncol(known_hannigan)-1)]/rowSums(known_hannigan[,1:(ncol(known_hannigan)-1)]) 


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
known_feng = known_feng[,colMeans(known_feng[,1:ncol(known_feng)]==0)<filtering_percent] 
known_feng[,1:(ncol(known_feng)-1)] = known_feng[,1:(ncol(known_feng)-1)]/rowSums(known_feng[,1:(ncol(known_feng)-1)]) 


known_vogtmann = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_vogtmann.csv",col_names = TRUE)
known_vogtmann = as.data.frame(known_vogtmann)
rownames(known_vogtmann) = known_vogtmann[,1]
known_vogtmann = known_vogtmann[,-1]
metadata_vogtmann <- read_csv("/home1/yilingao/proj/abundance/metadata_vogtmann.csv",col_names = T)
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 0] = "control"
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 1] = "CRC"
known_vogtmann = known_vogtmann[metadata_vogtmann$embl_id,]
known_vogtmann$status = metadata_vogtmann$casectl
#renormalization
known_vogtmann = known_vogtmann[,colMeans(known_vogtmann[,1:ncol(known_vogtmann)]==0)<filtering_percent]
known_vogtmann[,1:(ncol(known_vogtmann)-1)] = known_vogtmann[,1:(ncol(known_vogtmann)-1)]/rowSums(known_vogtmann[,1:(ncol(known_vogtmann)-1)]) 


known_zeller = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_zeller.csv",col_names = TRUE)
known_zeller = as.data.frame(known_zeller)
rownames(known_zeller) = known_zeller[,1]
known_zeller = known_zeller[,-1]
metadata_zeller = read_csv("/home1/yilingao/proj/abundance/Zeller_CRC.csv")
known_zeller = known_zeller[order(rownames(known_zeller)),]
metadata_zeller = metadata_zeller[order(metadata_zeller$No),]
known_zeller$status = metadata_zeller$study_condition
#renormalization
known_zeller = known_zeller[,colMeans(known_zeller[,1:ncol(known_zeller)]==0)<filtering_percent]
known_zeller[,1:(ncol(known_zeller)-1)] = known_zeller[,1:(ncol(known_zeller)-1)]/rowSums(known_zeller[,1:(ncol(known_zeller)-1)]) 


known_thomas = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_thomas.csv",col_names = TRUE)
known_thomas = as.data.frame(known_thomas)
rownames(known_thomas) = known_thomas[,1]
known_thomas = known_thomas[,-1]
metadata_thomas <- read_csv("/home1/yilingao/proj/abundance/metadata_thomas.csv",col_names = F)
metadata_thomas$X1[metadata_thomas$X1 == "CTR"] = "control"
known_thomas = known_thomas[metadata_thomas$X3,]
known_thomas$status = metadata_thomas$X1
#renormalization
known_thomas = known_thomas[,colMeans(known_thomas[,1:ncol(known_thomas)]==0)<filtering_percent]
known_thomas[,1:(ncol(known_thomas)-1)] = known_thomas[,1:(ncol(known_thomas)-1)]/rowSums(known_thomas[,1:(ncol(known_thomas)-1)]) 

set.seed(101)

###################################################################
auc_known_LODO = as.data.frame(matrix(NA,ncol = 1, nrow = 12))

features = unique(c(colnames(known_hannigan),colnames(known_feng),colnames(known_vogtmann),colnames(known_yu),colnames(known_thomas)))

###############################################################################
##test on zeller

###########function for selecting only training species
feature_selection = function(training){
for (i in 1:length(features)) {
  if(features[i] %!in% colnames(training)){
    training = cbind(training,0)
    colnames(training)[ncol(training)] = features[i]
  }
}
training = training[,order(colnames(training))]
return(training)
}

known_1 = feature_selection(known_yu)
known_2 = feature_selection(known_hannigan)
known_3 = feature_selection(known_feng)
known_4 = feature_selection(known_vogtmann)
known_5 = feature_selection(known_thomas)

training = rbind(known_1,known_2,known_3,known_4,known_5)
training[,1:(ncol(training)-1)] = training[,1:(ncol(training)-1)]/rowSums(training[,1:(ncol(training)-1)]) #renormalize
test = known_zeller
for (i in 1:ncol(training)) {
  if(colnames(training)[i] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[i]
  }
}
test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
for (i in 1:10){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 1000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred <- predict(rf, test, type="prob")
auc_known_LODO[i,1] <- auc(test$status, pred[,1])
}

auc_known_LODO[11,1] <- mean(auc_known_LODO[1:10,1])
auc_known_LODO[12,1] <- sd(auc_known_LODO[1:10,1])

write.csv(auc_known_LODO, "/home1/yilingao/proj/abundance/auc_known_LODO/bracken/bracken_test_zeller_1000trees.csv")

