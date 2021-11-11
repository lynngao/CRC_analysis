library(readr)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)

yu = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_yu.csv",col_names = TRUE)
yu = as.data.frame(yu)
rownames(yu) = yu[,1]
yu = yu[,-1]
metadata_yu <- read_csv("/home1/yilingao/proj/abundance/metadata_yu.csv",col_names = FALSE)
metadata_yu$X16[metadata_yu$X16 == "CTR"] = "control"
yu = yu[metadata_yu$X1,]
yu$status = metadata_yu$X16
#renormalization
filtering_percent = 0.9
yu = yu[,colMeans(yu[,1:ncol(yu)]==0)<filtering_percent]
yu[,1:(ncol(yu)-1)] = yu[,1:(ncol(yu)-1)]/rowSums(yu[,1:(ncol(yu)-1)]) 


hannigan = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_hannigan.csv",col_names = TRUE)
hannigan = as.data.frame(hannigan)
rownames(hannigan) = hannigan[,1]
hannigan = hannigan[,-1]
metadata_hannigan <- read_csv("/home1/yilingao/proj/abundance/metadata_hannigan.csv",col_names = FALSE)
metadata_hannigan$X5[metadata_hannigan$X5 == "Cancer"] = "CRC"
metadata_hannigan$X5[metadata_hannigan$X5 == "Healthy"] = "control"
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665060",]
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665121",]
hannigan = hannigan[metadata_hannigan$X3,]
hannigan$status = metadata_hannigan$X5
##exclude adenoma samples
hannigan = subset(hannigan,hannigan$status!=c("Adenoma"))
#renormalization
hannigan = hannigan[,colMeans(hannigan[,1:ncol(hannigan)]==0)<filtering_percent]
hannigan[,1:(ncol(hannigan)-1)] = hannigan[,1:(ncol(hannigan)-1)]/rowSums(hannigan[,1:(ncol(hannigan)-1)]) 


feng = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_feng.csv",col_names = TRUE)
feng = as.data.frame(feng)
rownames(feng) = feng[,1]
feng = feng[,-1]
metadata_feng <- read_csv("/home1/yilingao/proj/abundance/metadata_feng.csv",col_names = T)
metadata_feng$X4[metadata_feng$X4 == "carcinoma"] = "CRC"
metadata_feng$X4[metadata_feng$X4 == "controls"] = "control"
feng = feng[metadata_feng$fastq,]
feng$status = metadata_feng$X4
##exclude adenoma samples
feng = subset(feng,feng$status!=c("advanced adenoma"))
#renormalization
feng = feng[,colMeans(feng[,1:ncol(feng)]==0)<filtering_percent] 
feng[,1:(ncol(feng)-1)] = feng[,1:(ncol(feng)-1)]/rowSums(feng[,1:(ncol(feng)-1)]) 


vogtmann = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_vogtmann.csv",col_names = TRUE)
vogtmann = as.data.frame(vogtmann)
rownames(vogtmann) = vogtmann[,1]
vogtmann = vogtmann[,-1]
metadata_vogtmann <- read_csv("/home1/yilingao/proj/abundance/metadata_vogtmann.csv",col_names = T)
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 0] = "control"
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 1] = "CRC"
vogtmann = vogtmann[metadata_vogtmann$embl_id,]
vogtmann$status = metadata_vogtmann$casectl
#renormalization
vogtmann = vogtmann[,colMeans(vogtmann[,1:ncol(vogtmann)]==0)<filtering_percent]
vogtmann[,1:(ncol(vogtmann)-1)] = vogtmann[,1:(ncol(vogtmann)-1)]/rowSums(vogtmann[,1:(ncol(vogtmann)-1)]) 


zeller = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_zeller.csv",col_names = TRUE)
zeller = as.data.frame(zeller)
rownames(zeller) = zeller[,1]
zeller = zeller[,-1]
metadata_zeller = read_csv("/home1/yilingao/proj/abundance/Zeller_CRC.csv")
zeller = zeller[order(rownames(zeller)),]
metadata_zeller = metadata_zeller[order(metadata_zeller$No),]
zeller$status = metadata_zeller$study_condition
#renormalization
zeller = zeller[,colMeans(zeller[,1:ncol(zeller)]==0)<filtering_percent]
zeller[,1:(ncol(zeller)-1)] = zeller[,1:(ncol(zeller)-1)]/rowSums(zeller[,1:(ncol(zeller)-1)]) 


thomas = read_csv("/home1/yilingao/proj/abundance/bracken_abundance_table/bracken_thomas.csv",col_names = TRUE)
thomas = as.data.frame(thomas)
rownames(thomas) = thomas[,1]
thomas = thomas[,-1]
metadata_thomas <- read_csv("/home1/yilingao/proj/abundance/metadata_thomas.csv",col_names = F)
metadata_thomas$X1[metadata_thomas$X1 == "CTR"] = "control"
thomas = thomas[metadata_thomas$X3,]
thomas$status = metadata_thomas$X1
#renormalization
thomas = thomas[,colMeans(thomas[,1:ncol(thomas)]==0)<filtering_percent]
thomas[,1:(ncol(thomas)-1)] = thomas[,1:(ncol(thomas)-1)]/rowSums(thomas[,1:(ncol(thomas)-1)]) 

set.seed(101)

#############################cross-cohort known auc#########################################################
auc_cross_cohort = as.data.frame(matrix(NA,ncol = 5, nrow = 32))
colnames(auc_cross_cohort) = c("training:vogtmann, test:yu", "training:vogtmann, test:hannigan", "training:vogtmann, test:feng", "training:vogtmann, test:zeller","training:vogtmann, test:thomas")
rownames(auc_cross_cohort) = c(1:30,"mean","sd")

###########random forest function
rf_model = function(training, training_status, test, test_status, dataframe, col){

for (j in 1:ncol(training)) {
  if(colnames(training)[j] %!in% colnames(test)){
    test = cbind(test,0)
    colnames(test)[ncol(test)] = colnames(training)[j]
  }
}

training$status = training_status
test$status = test_status

for (i in 1:30){
rf <- train(
  status ~ ., data = training, method = "rf", ntree = 1000, metric = "ROC",
  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
)
pred <- predict(rf, test, type="prob")
dataframe[i,col] <- auc(test$status, pred[,1])
}
return(dataframe)
}


training = vogtmann[1:(ncol(vogtmann)-1)]
training_status = vogtmann$status
test1 = yu[1:(ncol(yu)-1)]
test1_status = yu$status
test2 = hannigan[1:(ncol(hannigan)-1)]
test2_status = hannigan$status
test3 = feng[1:(ncol(feng)-1)]
test3_status = feng$status
test4 = zeller[1:(ncol(zeller)-1)]
test4_status = zeller$status
test5 = thomas[1:(ncol(thomas)-1)]
test5_status = thomas$status

auc_cross_cohort <- rf_model(training, training_status, test1, test1_status, auc_cross_cohort, 1)
auc_cross_cohort <- rf_model(training, training_status, test2, test2_status, auc_cross_cohort, 2)
auc_cross_cohort <- rf_model(training, training_status, test3, test3_status, auc_cross_cohort, 3)
auc_cross_cohort <- rf_model(training, training_status, test4, test4_status, auc_cross_cohort, 4)
auc_cross_cohort <- rf_model(training, training_status, test5, test5_status, auc_cross_cohort, 5)


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


write.csv(auc_cross_cohort, "/home1/yilingao/proj/abundance/auc_cross_datasets/bracken/bracken_vogtmann_1000trees.csv")

