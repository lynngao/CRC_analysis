library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)

known_yu = readRDS("/Users/lynngao/Desktop/metadata/centrifuge_abundance_yu.rds")
metadata_yu <- read_csv("/Users/lynngao/Desktop/metadata/metadata_yu.csv",col_names = FALSE)
metadata_yu$X16[metadata_yu$X16 == "CTR"] = "control"
known_yu$status = metadata_yu$X16

known_hannigan = readRDS("/Users/lynngao/Desktop/metadata/centrifuge_abundance_hannigan.rds")
metadata_hannigan <- read_csv("/Users/lynngao/Desktop/metadata/metadata_hannigan.csv",col_names = FALSE)
metadata_hannigan$X5[metadata_hannigan$X5 == "Cancer"] = "CRC"
metadata_hannigan$X5[metadata_hannigan$X5 == "Healthy"] = "control"
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665060",]
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665121",]
known_hannigan$status = metadata_hannigan$X5
##exclude adenoma samples
known_hannigan = subset(known_hannigan,known_hannigan$status!=c("Adenoma"))

known_feng = readRDS("/Users/lynngao/Desktop/metadata/centrifuge_abundance_feng.rds")
metadata_feng <- read_csv("/Users/lynngao/Desktop/metadata/metadata_feng.csv",col_names = T)
colnames(metadata_feng)[4] = "X4"
metadata_feng$X4[metadata_feng$X4 == "carcinoma"] = "CRC"
metadata_feng$X4[metadata_feng$X4 == "controls"] = "control"
known_feng$status = metadata_feng$X4
##exclude adenoma samples
known_feng = subset(known_feng,known_feng$status!=c("advanced adenoma"))

known_vogtmann = readRDS("/Users/lynngao/Desktop/metadata/centrifuge_abundance_vogtmann.rds")
metadata_vogtmann <- read_csv("/Users/lynngao/Desktop/metadata/metadata_vogtmann.csv",col_names = T)
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 0] = "control"
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 1] = "CRC"
known_vogtmann$status = metadata_vogtmann$casectl

known_zeller = readRDS("/Users/lynngao/Desktop/metadata/centrifuge_abundance_zeller.rds")
metadata_zeller = read_csv("/Users/lynngao/Desktop/metadata/Zeller_CRC.csv")
known_zeller = known_zeller[order(rownames(known_zeller)),]
metadata_zeller = metadata_zeller[order(metadata_zeller$No),]
known_zeller$status = metadata_zeller$study_condition

auc_res <- function(dataset, ds) {
  for (i in 1:nrow(dataset)) {
    if (dataset[i,1]%in%rownames(ds)) {
      dataset$status[i] = ds$status[rownames(ds) == dataset[i,1]]
    }
  }
  auc(dataset$status, dataset$CRC)
}


AUC_LODO = as.data.frame(matrix(NA,nrow = 3, ncol = 6))
colnames(AUC_LODO) = c("yu","hannigan", "feng", "vogtmann","zeller", "thomas")
rownames(AUC_LODO) = c("known", "LOSO","known_cross")

##################AUC for known species from LODO prediction probabilities#################
##yu
res_yu = NA
for (i in 1:10){
  yu = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/LODO/pred_yu",i,".csv", sep="")))
  res_temp = auc_res(yu,known_yu)
  res_yu = c(res_yu, res_temp)
}
AUC_LODO[1,1] = mean(res_yu[-1])

##hannigan
res_hannigan = NA
for (i in 1:10){
  hannigan = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/LODO/pred_hannigan",i,".csv", sep="")))
  res_temp = auc_res(hannigan,known_hannigan)
  res_hannigan = c(res_hannigan, res_temp)
}
AUC_LODO[1,2] = mean(res_hannigan[-1])

##feng
res_feng = NA
for (i in 1:10){
  feng = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/LODO/pred_feng",i,".csv", sep="")))
  res_temp = auc_res(feng,known_feng)
  res_feng = c(res_feng, res_temp)
}
AUC_LODO[1,3] = mean(res_feng[-1])

##vogtmann
res_vogtmann = NA
for (i in 1:10){
  vogtmann = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/LODO/pred_vogtmann",i,".csv", sep="")))
  res_temp = auc_res(vogtmann,known_vogtmann)
  res_vogtmann = c(res_vogtmann, res_temp)
}
AUC_LODO[1,4] = mean(res_vogtmann[-1])

##zeller
res_zeller = NA
for (i in 1:10){
  zeller = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/LODO/pred_zeller",i,".csv", sep="")))
  res_temp = auc_res(zeller,known_zeller)
  res_zeller = c(res_zeller, res_temp)
}
AUC_LODO[1,5] = mean(res_zeller[-1])

##thomas
res_thomas = NA
for (i in 1:10){
  thomas = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/LODO/pred_thomas",i,".csv", sep="")))
  res_temp = auc_res(thomas,known_thomas)
  res_thomas = c(res_thomas, res_temp)
}
AUC_LODO[1,6] = mean(res_thomas[-1])
########################################

##########test:yu
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_yu.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_yu.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_yu.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_yu.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_yu.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_yu.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_yu.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_yu.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_yu.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_yu.csv"))

test = known_yu

hannigan_known$status = test$status
feng_known$status = test$status
vogtmann_known$status = test$status
zeller_known$status = test$status
thomas_known$status = test$status
hannigan_unknown$status = test$status
feng_unknown$status = test$status
vogtmann_unknown$status = test$status
zeller_unknown$status = test$status
thomas_unknown$status = test$status

##training:hannigan
hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(hannigan_AUC)){
	hannigan_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(hannigan_AUC)[-i]])
	hannigan_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(hannigan_AUC)[-i]])
}
hannigan_AUC[hannigan_AUC < 0.5] = 0.5

##training:feng
feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(feng_AUC) = c("known_AUC","unknown_AUC")
rownames(feng_AUC) = feng_known[,1]
for (i in 1:nrow(feng_AUC)){
	feng_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(feng_AUC)[-i]])
	feng_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(feng_AUC)[-i]])
}
feng_AUC[feng_AUC < 0.5] = 0.5

##training:vogtmann
vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(vogtmann_AUC)){
	vogtmann_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(vogtmann_AUC)[-i]])
	vogtmann_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(vogtmann_AUC)[-i]])
}
vogtmann_AUC[vogtmann_AUC < 0.5] = 0.5

##training:zeller
zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(zeller_AUC)){
	zeller_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(zeller_AUC)[-i]])
	zeller_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(zeller_AUC)[-i]])
}
zeller_AUC[zeller_AUC < 0.5] = 0.5

##training:thomas
thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(thomas_AUC)){
	thomas_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(thomas_AUC)[-i]])
	thomas_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(thomas_AUC)[-i]])
}
thomas_AUC[thomas_AUC < 0.5] = 0.5


pred_known = ((hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3]) / (hannigan_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] - 2.5) 
AUC_LODO[3,1] = auc(test$status, pred_known) #0.8574
pred_unknown = ((hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (hannigan_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 2.5) 
auc(test$status, pred_unknown) #0.8058
pred = ((hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (hannigan_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] + hannigan_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 5) 
AUC_LODO[2,1] = auc(test$status, pred) #0.8766




##########test:hannigan
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_hannigan.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_hannigan.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_hannigan.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_hannigan.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_hannigan.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_hannigan.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_hannigan.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_hannigan.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_hannigan.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_hannigan.csv"))

test = known_hannigan

##training:yu
yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(yu_AUC) = c("known_AUC","unknown_AUC")
rownames(yu_AUC) = yu_known[,1]
for (i in 1:nrow(yu_AUC)){
	yu_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(yu_AUC)[-i]])
	yu_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(yu_AUC)[-i]])
}
yu_AUC[yu_AUC < 0.5] = 0.5

##training:feng
feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(feng_AUC) = c("known_AUC","unknown_AUC")
rownames(feng_AUC) = feng_known[,1]
for (i in 1:nrow(feng_AUC)){
	feng_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(feng_AUC)[-i]])
	feng_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(feng_AUC)[-i]])
}
feng_AUC[feng_AUC < 0.5] = 0.5

##training:vogtmann
vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(vogtmann_AUC)){
	vogtmann_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(vogtmann_AUC)[-i]])
	vogtmann_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(vogtmann_AUC)[-i]])
}
vogtmann_AUC[vogtmann_AUC < 0.5] = 0.5

##training:zeller
zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(zeller_AUC)){
	zeller_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(zeller_AUC)[-i]])
	zeller_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(zeller_AUC)[-i]])
}
zeller_AUC[zeller_AUC < 0.5] = 0.5

##training:thomas
thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(thomas_AUC)){
	thomas_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(thomas_AUC)[-i]])
	thomas_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(thomas_AUC)[-i]])
}
thomas_AUC[thomas_AUC < 0.5] = 0.5


pred_known = ((yu_AUC[,1]-0.5)*yu_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3]) / (yu_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] - 2.5) 
AUC_LODO[3,2] = auc(test$status, pred_known) #0.6362
pred_unknown = ((yu_AUC[,2]-0.5)*yu_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 2.5) 
auc(test$status, pred_unknown) #0.5595
pred = ((yu_AUC[,1]-0.5)*yu_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3] + (yu_AUC[,2]-0.5)*yu_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] + yu_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 5) 
AUC_LODO[2,2] = auc(test$status, pred) #0.627



##########test:feng
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_feng.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_feng.csv"))
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_feng.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_feng.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_feng.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_feng.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_feng.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_feng.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_feng.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_feng.csv"))

test = known_feng

##training:yu
yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(yu_AUC) = c("known_AUC","unknown_AUC")
rownames(yu_AUC) = yu_known[,1]
for (i in 1:nrow(yu_AUC)){
	yu_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(yu_AUC)[-i]])
	yu_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(yu_AUC)[-i]])
}
yu_AUC[yu_AUC < 0.5] = 0.5

##training:hannigan
hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(hannigan_AUC)){
	hannigan_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(hannigan_AUC)[-i]])
	hannigan_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(hannigan_AUC)[-i]])
}
hannigan_AUC[hannigan_AUC < 0.5] = 0.5

##training:vogtmann
vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(vogtmann_AUC)){
	vogtmann_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(vogtmann_AUC)[-i]])
	vogtmann_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(vogtmann_AUC)[-i]])
}
vogtmann_AUC[vogtmann_AUC < 0.5] = 0.5

##training:zeller
zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(zeller_AUC)){
	zeller_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(zeller_AUC)[-i]])
	zeller_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(zeller_AUC)[-i]])
}
zeller_AUC[zeller_AUC < 0.5] = 0.5

##training:thomas
thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(thomas_AUC)){
	thomas_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(thomas_AUC)[-i]])
	thomas_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(thomas_AUC)[-i]])
}
thomas_AUC[thomas_AUC < 0.5] = 0.5


pred_known = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] - 2.5) 
AUC_LODO[3,3] = auc(test$status, pred_known) #0.893
pred_unknown = ((yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,2] + hannigan_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 2.5) 
auc(test$status, pred_unknown) #0.8506
pred = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3] + (yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] + yu_AUC[,2] + hannigan_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 5) 
AUC_LODO[2,3] = auc(test$status, pred) #0.9017


##########test:vogtmann
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_vogtmann.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_vogtmann.csv"))
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_vogtmann.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_vogtmann.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_vogtmann.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_vogtmann.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_vogtmann.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_vogtmann.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_vogtmann.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_vogtmann.csv"))

test = known_vogtmann

##training:yu
yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(yu_AUC) = c("known_AUC","unknown_AUC")
rownames(yu_AUC) = yu_known[,1]
for (i in 1:nrow(yu_AUC)){
	yu_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(yu_AUC)[-i]])
	yu_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(yu_AUC)[-i]])
}
yu_AUC[yu_AUC < 0.5] = 0.5

##training:hannigan
hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(hannigan_AUC)){
	hannigan_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(hannigan_AUC)[-i]])
	hannigan_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(hannigan_AUC)[-i]])
}
hannigan_AUC[hannigan_AUC < 0.5] = 0.5

##training:feng
feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(feng_AUC) = c("known_AUC","unknown_AUC")
rownames(feng_AUC) = feng_known[,1]
for (i in 1:nrow(feng_AUC)){
	feng_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(feng_AUC)[-i]])
	feng_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(feng_AUC)[-i]])
}
feng_AUC[feng_AUC < 0.5] = 0.5

##training:zeller
zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(zeller_AUC)){
	zeller_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(zeller_AUC)[-i]])
	zeller_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(zeller_AUC)[-i]])
}
zeller_AUC[zeller_AUC < 0.5] = 0.5

##training:thomas
thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(thomas_AUC)){
	thomas_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(thomas_AUC)[-i]])
	thomas_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(thomas_AUC)[-i]])
}
thomas_AUC[thomas_AUC < 0.5] = 0.5


pred_known = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + feng_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] - 2.5) 
AUC_LODO[3,4] = auc(test$status, pred_known) #0.7777
pred_unknown = ((yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,2] + hannigan_AUC[,2] + feng_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 2.5) 
auc(test$status, pred_unknown) #0.743
pred = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3] + (yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + feng_AUC[,1] + zeller_AUC[,1] + thomas_AUC[,1] + yu_AUC[,2] + hannigan_AUC[,2] + feng_AUC[,2] + zeller_AUC[,2] + thomas_AUC[,2] - 5) 
AUC_LODO[2,4] = auc(test$status, pred) #0.7999



##########test:zeller
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_zeller.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_zeller.csv"))
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_zeller.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_zeller.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_zeller.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_zeller.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_zeller.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_zeller.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_zeller.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_zeller.csv"))

test = known_zeller

##training:yu
yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(yu_AUC) = c("known_AUC","unknown_AUC")
rownames(yu_AUC) = yu_known[,1]
for (i in 1:nrow(yu_AUC)){
	yu_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(yu_AUC)[-i]])
	yu_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(yu_AUC)[-i]])
}
yu_AUC[yu_AUC < 0.5] = 0.5

##training:hannigan
hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(hannigan_AUC)){
	hannigan_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(hannigan_AUC)[-i]])
	hannigan_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(hannigan_AUC)[-i]])
}
hannigan_AUC[hannigan_AUC < 0.5] = 0.5

##training:feng
feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(feng_AUC) = c("known_AUC","unknown_AUC")
rownames(feng_AUC) = feng_known[,1]
for (i in 1:nrow(feng_AUC)){
	feng_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(feng_AUC)[-i]])
	feng_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(feng_AUC)[-i]])
}
feng_AUC[feng_AUC < 0.5] = 0.5

##training:vogtmann
vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(vogtmann_AUC)){
	vogtmann_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(vogtmann_AUC)[-i]])
	vogtmann_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(vogtmann_AUC)[-i]])
}
vogtmann_AUC[vogtmann_AUC < 0.5] = 0.5

##training:thomas
thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(thomas_AUC)){
	thomas_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(thomas_AUC)[-i]])
	thomas_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(thomas_AUC)[-i]])
}
thomas_AUC[thomas_AUC < 0.5] = 0.5


pred_known = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + thomas_AUC[,1] - 2.5) 
AUC_LODO[3,5] = auc(test$status, pred_known) #0.7844
pred_unknown = ((yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,2] + hannigan_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + thomas_AUC[,2] - 2.5) 
auc(test$status, pred_unknown) #0.6322
pred = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (thomas_AUC[,1]-0.5)*thomas_known[,3] + (yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (thomas_AUC[,2]-0.5)*thomas_unknown[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + thomas_AUC[,1] + yu_AUC[,2] + hannigan_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + thomas_AUC[,2] - 5) 
AUC_LODO[2,5] = auc(test$status, pred) #0.7976



##########test:thomas
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_thomas.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_thomas.csv"))
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_thomas.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_thomas.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_thomas.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_thomas.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_thomas.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_thomas.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_thomas.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_thomas.csv"))

test = known_thomas

##training:yu
yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(yu_AUC) = c("known_AUC","unknown_AUC")
rownames(yu_AUC) = yu_known[,1]
for (i in 1:nrow(yu_AUC)){
	yu_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(yu_AUC)[-i]])
	yu_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(yu_AUC)[-i]])
}
yu_AUC[yu_AUC < 0.5] = 0.5

##training:hannigan
hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(hannigan_AUC)){
	hannigan_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(hannigan_AUC)[-i]])
	hannigan_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(hannigan_AUC)[-i]])
}
hannigan_AUC[hannigan_AUC < 0.5] = 0.5

##training:feng
feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(feng_AUC) = c("known_AUC","unknown_AUC")
rownames(feng_AUC) = feng_known[,1]
for (i in 1:nrow(feng_AUC)){
	feng_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(feng_AUC)[-i]])
	feng_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(feng_AUC)[-i]])
}
feng_AUC[feng_AUC < 0.5] = 0.5

##training:vogtmann
vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(vogtmann_AUC)){
	vogtmann_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(vogtmann_AUC)[-i]])
	vogtmann_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(vogtmann_AUC)[-i]])
}
vogtmann_AUC[vogtmann_AUC < 0.5] = 0.5

##training:zeller
zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(zeller_AUC)){
	zeller_AUC[i,1] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(zeller_AUC)[-i]])
	zeller_AUC[i,2] = auc(test$status[rownames(test)%in%rownames(zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(zeller_AUC)[-i]])
}
zeller_AUC[zeller_AUC < 0.5] = 0.5


pred_known = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] - 2.5) 
AUC_LODO[3,6] = auc(test$status, pred_known) #0.7702
pred_unknown = ((yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3]) / (yu_AUC[,2] + hannigan_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] - 2.5) 
auc(test$status, pred_unknown) #0.7156
pred = ((yu_AUC[,1]-0.5)*yu_known[,3] + (hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (feng_AUC[,1]-0.5)*feng_known[,3] + (vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (zeller_AUC[,1]-0.5)*zeller_known[,3] + (yu_AUC[,2]-0.5)*yu_unknown[,3] + (hannigan_AUC[,2]-0.5)*hannigan_unknown[,3] + (feng_AUC[,2]-0.5)*feng_unknown[,3] + (vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3] + (zeller_AUC[,2]-0.5)*zeller_unknown[,3]) / (yu_AUC[,1] + hannigan_AUC[,1] + feng_AUC[,1] + vogtmann_AUC[,1] + zeller_AUC[,1] + yu_AUC[,2] + hannigan_AUC[,2] + feng_AUC[,2] + vogtmann_AUC[,2] + zeller_AUC[,2] - 5) 
AUC_LODO[2,6] = auc(test$status, pred) #0.7705


AUC_LODO = round(AUC_LODO, digits = 2)
write.csv(AUC_LODO, "/Users/lynngao/Desktop/results/AUC_LODO.csv")




