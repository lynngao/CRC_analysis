library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(PRROC)

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

known_thomas = readRDS("/Users/lynngao/Desktop/metadata/centrifuge_abundance_thomas.rds")
metadata_thomas <- read_csv("/Users/lynngao/Desktop/metadata/metadata_thomas.csv",col_names = F)
metadata_thomas$X1[metadata_thomas$X1 == "CTR"] = "control"
known_thomas$status = metadata_thomas$X1

auprc <- function(dataset, ds) {
  dataset[,"known_CRC"] = round(dataset[,"known_CRC"], digits = 2)
  dataset[,"unknown_CRC"] = round(dataset[,"unknown_CRC"], digits = 2)
  dataset[,"known_control"] = 1 - dataset[,"known_CRC"]
  dataset[,"unknown_control"] = 1 - dataset[,"unknown_CRC"]
  for (i in 1:nrow(dataset)) {
    if (dataset[i,1]%in%rownames(ds)) {
      dataset$status[i] = ds$status[rownames(ds) == dataset[i,1]]
    }
  }
  dataset$status = ifelse(dataset$status == "control", 0, 1)
  auprc_known = pr.curve(scores.class0 = dataset[,"known_CRC"], weights.class0 = dataset$status, curve = T)$auc.integral
  auprc_unknown = pr.curve(scores.class0 = dataset[,"unknown_CRC"], weights.class0 = dataset$status, curve = T)$auc.integral
  return(list(auprc_known, auprc_unknown))
}

auprc_combined <- function(dataset, ds) {
  dataset[,"CRC"] = round(dataset[,"CRC"], digits = 2)
  dataset[,"control"] = 1 - dataset[,"CRC"]
  for (i in 1:nrow(dataset)) {
    if (dataset[i,1]%in%rownames(ds)) {
      dataset$status[i] = ds$status[rownames(ds) == dataset[i,1]]
    }
  }
  dataset$status = ifelse(dataset$status == "control", 0, 1)
  auprc_combined = pr.curve(scores.class0 = dataset[,"CRC"], weights.class0 = dataset$status, curve = T)$auc.integral
  return(auprc_combined)
}


auprc_sample <- function(dataset, ds, ratio){
  dataset[,"known_CRC"] = round(dataset[,"known_CRC"], digits = 2)
  dataset[,"unknown_CRC"] = round(dataset[,"unknown_CRC"], digits = 2)
  dataset[,"known_control"] = 1 - dataset[,"known_CRC"]
  dataset[,"unknown_control"] = 1 - dataset[,"unknown_CRC"]
  for (j in 1:nrow(dataset)) {
    if (dataset[j,1]%in%rownames(ds)) {
      dataset$status[j] = ds$status[rownames(ds) == dataset[j,1]]
    }
  }
  dataset$status = ifelse(dataset$status == "control", 0, 1)
  temp_prob = as.data.frame(matrix(NA, ncol = 4, nrow = nrow(dataset)))
  colnames(temp_prob) = c("known_AUC","unknown_AUC", "avg_control","avg_CRC")
  rownames(temp_prob) = dataset[,1]
  for (k in 1:nrow(temp_prob)){
    temp_prob[k,1] = auc(dataset$status[dataset[,1] != rownames(ds)[k]],dataset[dataset[,1] != rownames(ds)[k],"known_CRC"])
    temp_prob[k,2] = auc(dataset$status[dataset[,1] != rownames(ds)[k]],dataset[dataset[,1] != rownames(ds)[k],"unknown_CRC"])
  }
  temp_prob[temp_prob < 0.5] = 0.5
  temp_prob$avg_control = (( temp_prob[,1]-0.5)*dataset[,"known_control"] + (temp_prob[,2]-0.5)*dataset[,"unknown_control"])/(temp_prob[,1]+temp_prob[,2]-2*0.5)
  temp_prob$avg_CRC = (( temp_prob[,1]-0.5)*dataset[,"known_CRC"] + (temp_prob[,2]-0.5)*dataset[,"unknown_CRC"])/(temp_prob[,1]+temp_prob[,2]-2*0.5)
  temp_prob$status = dataset$status
  pr.curve(scores.class0 = temp_prob[,"avg_CRC"], weights.class0 = temp_prob$status, curve = T)$auc.integral
}

AUPRC_within = as.data.frame(matrix(NA,nrow = 4, ncol = 6))
colnames(AUPRC_within) = c("yu","hannigan", "feng", "vogtmann","zeller", "thomas")
rownames(AUPRC_within) = c("known", "unknown", "LOSO", "combined")

#################within-dataset#################

##yu
res_yu_known = NA
res_yu_unknown = NA
res_yu_sample = NA
for (i in 1:30){
  yu = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_yu",i,".csv", sep="")))
  res_known_temp = unlist(auprc(yu,known_yu)[1])
  res_unknown_temp = unlist(auprc(yu,known_yu)[2])
  res_yu_known = c(res_yu_known, res_known_temp)
  res_yu_unknown = c(res_yu_unknown, res_unknown_temp)
  res_yu_temp_sample = auprc_sample(yu,known_yu, ratio_hannigan)
  res_yu_sample = c(res_yu_sample, res_yu_temp_sample)
}
AUPRC_within[1,1] = mean(res_yu_known[-1])
AUPRC_within[2,1]= mean(res_yu_unknown[-1])
AUPRC_within[3,1] = mean(res_yu_sample[-1])
#combined AUPRC
res_yu_combined = NA
for (i in 1:30){
  yu = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/combined/pred_yu",i,".csv", sep="")))
  res_combined_temp = unlist(auprc_combined(yu,known_yu)[1])
  res_yu_combined = c(res_yu_combined, res_combined_temp)
}
AUPRC_within[4,1] = mean(res_yu_combined[-1])

##hannigan
res_hannigan_known = NA
res_hannigan_unknown = NA
res_hannigan_sample = NA
for (i in 1:30){
  hannigan = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_hannigan",i,".csv", sep="")))
  res_known_temp = unlist(auprc(hannigan,known_hannigan)[1])
  res_unknown_temp = unlist(auprc(hannigan,known_hannigan)[2])
  res_hannigan_known = c(res_hannigan_known, res_known_temp)
  res_hannigan_unknown = c(res_hannigan_unknown, res_unknown_temp)
  res_hannigan_temp_sample = auprc_sample(hannigan,known_hannigan, ratio_hannigan)
  res_hannigan_sample = c(res_hannigan_sample, res_hannigan_temp_sample)
}
AUPRC_within[1,2] = mean(res_hannigan_known[-1])
AUPRC_within[2,2]= mean(res_hannigan_unknown[-1])
AUPRC_within[3,2] = mean(res_hannigan_sample[-1])
#combined AUPRC
res_hannigan_combined = NA
for (i in 1:30){
  hannigan = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/combined/pred_hannigan",i,".csv", sep="")))
  res_combined_temp = unlist(auprc_combined(hannigan,known_hannigan)[1])
  res_hannigan_combined = c(res_hannigan_combined, res_combined_temp)
}
AUPRC_within[4,2] = mean(res_hannigan_combined[-1])

##feng
res_feng_known = NA
res_feng_unknown = NA
res_feng_sample = NA
for (i in 1:30){
  feng = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_feng",i,".csv", sep="")))
  res_known_temp = unlist(auprc(feng,known_feng)[1])
  res_unknown_temp = unlist(auprc(feng,known_feng)[2])
  res_feng_known = c(res_feng_known, res_known_temp)
  res_feng_unknown = c(res_feng_unknown, res_unknown_temp)
  res_feng_temp_sample = auprc_sample(feng,known_feng, ratio_hannigan)
  res_feng_sample = c(res_feng_sample, res_feng_temp_sample)
}
AUPRC_within[1,3] = mean(res_feng_known[-1])
AUPRC_within[2,3]= mean(res_feng_unknown[-1])
AUPRC_within[3,3] = mean(res_feng_sample[-1])
#combined AUPRC
res_feng_combined = NA
for (i in 1:30){
  feng = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/combined/pred_feng",i,".csv", sep="")))
  res_combined_temp = unlist(auprc_combined(feng,known_feng)[1])
  res_feng_combined = c(res_feng_combined, res_combined_temp)
}
AUPRC_within[4,3] = mean(res_feng_combined[-1])

##vogtmann
res_vogtmann_known = NA
res_vogtmann_unknown = NA
res_vogtmann_sample = NA
for (i in 1:30){
  vogtmann = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_vogtmann",i,".csv", sep="")))
  res_known_temp = unlist(auprc(vogtmann,known_vogtmann)[1])
  res_unknown_temp = unlist(auprc(vogtmann,known_vogtmann)[2])
  res_vogtmann_known = c(res_vogtmann_known, res_known_temp)
  res_vogtmann_unknown = c(res_vogtmann_unknown, res_unknown_temp)
  res_vogtmann_temp_sample = auprc_sample(vogtmann,known_vogtmann, ratio_hannigan)
  res_vogtmann_sample = c(res_vogtmann_sample, res_vogtmann_temp_sample)
}
AUPRC_within[1,4] = mean(res_vogtmann_known[-1])
AUPRC_within[2,4]= mean(res_vogtmann_unknown[-1])
AUPRC_within[3,4] = mean(res_vogtmann_sample[-1])
#combined AUPRC
res_vogtmann_combined = NA
for (i in 1:30){
  vogtmann = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/combined/pred_vogtmann",i,".csv", sep="")))
  res_combined_temp = unlist(auprc_combined(vogtmann,known_vogtmann)[1])
  res_vogtmann_combined = c(res_vogtmann_combined, res_combined_temp)
}
AUPRC_within[4,4] = mean(res_vogtmann_combined[-1])

##zeller
res_zeller_known = NA
res_zeller_unknown = NA
res_zeller_sample = NA
for (i in 1:30){
  zeller = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_zeller",i,".csv", sep="")))
  res_known_temp = unlist(auprc(zeller,known_zeller)[1])
  res_unknown_temp = unlist(auprc(zeller,known_zeller)[2])
  res_zeller_known = c(res_zeller_known, res_known_temp)
  res_zeller_unknown = c(res_zeller_unknown, res_unknown_temp)
  res_zeller_temp_sample = auprc_sample(zeller,known_zeller, ratio_hannigan)
  res_zeller_sample = c(res_zeller_sample, res_zeller_temp_sample)
}
AUPRC_within[1,5] = mean(res_zeller_known[-1])
AUPRC_within[2,5]= mean(res_zeller_unknown[-1])
AUPRC_within[3,5] = mean(res_zeller_sample[-1])
#combined AUPRC
res_zeller_combined = NA
for (i in 1:30){
  zeller = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/combined/pred_zeller",i,".csv", sep="")))
  res_combined_temp = unlist(auprc_combined(zeller,known_zeller)[1])
  res_zeller_combined = c(res_zeller_combined, res_combined_temp)
}
AUPRC_within[4,5] = mean(res_zeller_combined[-1])

##thomas
res_thomas_known = NA
res_thomas_unknown = NA
res_thomas_sample = NA
for (i in 1:30){
  thomas = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_thomas",i,".csv", sep="")))
  res_known_temp = unlist(auprc(thomas,known_thomas)[1])
  res_unknown_temp = unlist(auprc(thomas,known_thomas)[2])
  res_thomas_known = c(res_thomas_known, res_known_temp)
  res_thomas_unknown = c(res_thomas_unknown, res_unknown_temp)
  res_thomas_temp_sample = auprc_sample(thomas,known_thomas, ratio_hannigan)
  res_thomas_sample = c(res_thomas_sample, res_thomas_temp_sample)
}
AUPRC_within[1,6] = mean(res_thomas_known[-1])
AUPRC_within[2,6]= mean(res_thomas_unknown[-1])
AUPRC_within[3,6] = mean(res_thomas_sample[-1])
#combined AUPRC
res_thomas_combined = NA
for (i in 1:30){
  thomas = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/combined/pred_thomas",i,".csv", sep="")))
  res_combined_temp = unlist(auprc_combined(thomas,known_thomas)[1])
  res_thomas_combined = c(res_thomas_combined, res_combined_temp)
}
AUPRC_within[4,6] = mean(res_thomas_combined[-1])

AUPRC_within = round(AUPRC_within, digits = 2)
write.csv(AUPRC_within, "/Users/lynngao/Desktop/results/AUPRC/AUPRC_within1.csv")
