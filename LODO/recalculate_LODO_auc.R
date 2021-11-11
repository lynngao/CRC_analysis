library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

known_yu = readRDS("/home/rcf-40/yilingao/cmb/abundance/centrifuge_abundance_yu.rds")
unknown_yu = readRDS("/home/rcf-40/yilingao/cmb/abundance/unknown_abundance_yu.rds")
combined_yu = readRDS("/home/rcf-40/yilingao/cmb/abundance/combined_abundance_yu.rds")
metadata_yu <- read_csv("/home/rcf-40/yilingao/cmb/abundance/metadata_yu.csv",col_names = FALSE)
metadata_yu$X16[metadata_yu$X16 == "CTR"] = "control"
colnames(combined_yu)[(ncol(known_yu)+1):(ncol(combined_yu))] = paste0('u',colnames(combined_yu)[(ncol(known_yu)+1):(ncol(combined_yu))])
unknown_yu = as.data.frame(unknown_yu)
colnames(unknown_yu) = paste0('u',colnames(unknown_yu))
known_yu$status = metadata_yu$X16
combined_yu$status = metadata_yu$X16
unknown_yu$status = metadata_yu$X16



known_hannigan = readRDS("/home/rcf-40/yilingao/cmb/abundance/centrifuge_abundance_hannigan.rds")
unknown_hannigan = readRDS("/home/rcf-40/yilingao/cmb/abundance/unknown_abundance_hannigan.rds")
unknown_hannigan = as.data.frame(unknown_hannigan)
combined_hannigan = readRDS("/home/rcf-40/yilingao/cmb/abundance/combined_abundance_hannigan.rds")
metadata_hannigan <- read_csv("/home/rcf-40/yilingao/cmb/abundance/metadata_hannigan.csv",col_names = FALSE)
metadata_hannigan$X5[metadata_hannigan$X5 == "Cancer"] = "CRC"
metadata_hannigan$X5[metadata_hannigan$X5 == "Healthy"] = "control"
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665060",]
metadata_hannigan = metadata_hannigan[metadata_hannigan$X3 != "SRR5665121",]
colnames(combined_hannigan)[(ncol(known_hannigan)+1):(ncol(combined_hannigan))] = paste0('u',colnames(combined_hannigan)[(ncol(known_hannigan)+1):(ncol(combined_hannigan))])
known_hannigan$status = metadata_hannigan$X5
combined_hannigan$status = metadata_hannigan$X5
unknown_hannigan$status = metadata_hannigan$X5
##exclude adenoma samples
known_hannigan = subset(known_hannigan,known_hannigan$status!=c("Adenoma"))
combined_hannigan = subset(combined_hannigan,combined_hannigan$status!=c("Adenoma"))
unknown_hannigan = subset(unknown_hannigan,unknown_hannigan$status!=c("Adenoma"))


known_feng = readRDS("/home/rcf-40/yilingao/cmb/abundance/centrifuge_abundance_feng.rds")
unknown_feng = readRDS("/home/rcf-40/yilingao/cmb/abundance/unknown_abundance_feng.rds")
unknown_feng = as.data.frame(unknown_feng)
combined_feng = readRDS("/home/rcf-40/yilingao/cmb/abundance/combined_abundance_feng.rds")
metadata_feng <- read_csv("/home/rcf-40/yilingao/cmb/abundance/metadata_feng.csv",col_names = T)
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

known_vogtmann = readRDS("/home/rcf-40/yilingao/cmb/abundance/centrifuge_abundance_vogtmann.rds")
unknown_vogtmann = readRDS("/home/rcf-40/yilingao/cmb/abundance/unknown_abundance_vogtmann.rds")
unknown_vogtmann = as.data.frame(unknown_vogtmann)
combined_vogtmann = readRDS("/home/rcf-40/yilingao/cmb/abundance/combined_abundance_vogtmann.rds")
metadata_vogtmann <- read_csv("/home/rcf-40/yilingao/cmb/abundance/metadata_vogtmann.csv",col_names = T)
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 0] = "control"
metadata_vogtmann$casectl[metadata_vogtmann$casectl == 1] = "CRC"
colnames(combined_vogtmann)[(ncol(known_vogtmann)+1):(ncol(combined_vogtmann))] = paste0('u',colnames(combined_vogtmann)[(ncol(known_vogtmann)+1):(ncol(combined_vogtmann))])
colnames(unknown_vogtmann) = paste0('u',colnames(unknown_vogtmann))
known_vogtmann$status = metadata_vogtmann$casectl
unknown_vogtmann$status = metadata_vogtmann$casectl
combined_vogtmann$status = metadata_vogtmann$casectl

combined_zeller = readRDS("/home/rcf-40/yilingao/cmb/abundance/combined_abundance_zeller.rds")
unknown_zeller = readRDS("/home/rcf-40/yilingao/cmb/abundance/unknown_abundance_zeller.rds")
unknown_zeller = as.data.frame(unknown_zeller)
known_zeller = readRDS("/home/rcf-40/yilingao/cmb/abundance/centrifuge_abundance_zeller.rds")
metadata_zeller = read_csv("/home/rcf-40/yilingao/cmb/abundance/Zeller_CRC.csv")
combined_zeller = combined_zeller[order(rownames(combined_zeller)),]
unknown_zeller = unknown_zeller[order(rownames(unknown_zeller)),]
known_zeller = known_zeller[order(rownames(known_zeller)),]
metadata_zeller = metadata_zeller[order(metadata_zeller$No),]
combined_zeller$status = metadata_zeller$study_condition
unknown_zeller$status = metadata_zeller$study_condition
known_zeller$status = metadata_zeller$study_condition
colnames(combined_zeller)[(ncol(known_zeller)):(ncol(combined_zeller)-1)] = paste0('u',colnames(combined_zeller)[(ncol(known_zeller)):(ncol(combined_zeller)-1)])
colnames(unknown_zeller)[1:1252] = paste0('u',colnames(unknown_zeller)[1:1252])

known_thomas = readRDS("/home/rcf-40/yilingao/cmb/abundance/centrifuge_abundance_thomas.rds")
unknown_thomas = readRDS("/home/rcf-40/yilingao/cmb/abundance/unknown_abundance_thomas.rds")
unknown_thomas = as.data.frame(unknown_thomas)
combined_thomas = readRDS("/home/rcf-40/yilingao/cmb/abundance/combined_abundance_thomas.rds")
metadata_thomas <- read_csv("/home/rcf-40/yilingao/cmb/abundance/metadata_thomas.csv",col_names = F)
metadata_thomas$X1[metadata_thomas$X1 == "CTR"] = "control"
colnames(combined_thomas)[(ncol(known_thomas)+1):(ncol(combined_thomas))] = paste0('u',colnames(combined_thomas)[(ncol(known_thomas)+1):(ncol(combined_thomas))])
colnames(unknown_thomas) = paste0('u',colnames(unknown_thomas))
known_thomas$status = metadata_thomas$X1
unknown_thomas$status = metadata_thomas$X1
combined_thomas$status = metadata_thomas$X1

#############################################################################


##########test:yu
feng = read_csv("/home/rcf-40/yilingao/staging/pred/test_yu/pred_training_feng.csv")
hannigan = read_csv("/home/rcf-40/yilingao/staging/pred/test_yu/pred_training_hannigan.csv")
vogtmann = read_csv("/home/rcf-40/yilingao/staging/pred/test_yu/pred_training_vogtmann.csv")
zeller = read_csv("/home/rcf-40/yilingao/staging/pred/test_yu/pred_training_zeller.csv")
thomas = read_csv("/home/rcf-40/yilingao/staging/pred/test_yu/pred_training_thomas.csv")

test = known_yu

##average
pred = (feng[-1] + hannigan[-1] + zeller[-1] + thomas[-1] + vogtmann[-1]) / 5
auc(test$status, pred[,1]) //0.8549 //nzv:0.8616 
##weights based on cross-dataset AUC
0.65+0.8+0.78+0.71+0.8 = 3.74
pred1 = (0.65*hannigan[-1] + 0.8*feng[-1] + 0.78*vogtmann[-1] + 0.71*zeller[-1] + 0.8*thomas[-1] ) / 3.74
auc(test$status, pred1[,1]) //0.8556 //nzv=0.8629


##########test:hannigan
feng = read_csv("/home/rcf-40/yilingao/staging/pred/test_hannigan/pred_training_feng.csv")
yu = read_csv("/home/rcf-40/yilingao/staging/pred/test_hannigan/pred_training_yu.csv")
vogtmann = read_csv("/home/rcf-40/yilingao/staging/pred/test_hannigan/pred_training_vogtmann.csv")
zeller = read_csv("/home/rcf-40/yilingao/staging/pred/test_hannigan/pred_training_zeller.csv")
thomas = read_csv("/home/rcf-40/yilingao/staging/pred/test_hannigan/pred_training_thomas.csv")

test = known_hannigan

##average
pred = (feng[-1] + yu[-1] + zeller[-1] + thomas[-1] + vogtmann[-1]) / 5
auc(test$status, pred[,1]) //0.664 
##weights based on cross-dataset AUC
0.67+0.65+0.59+0.64+0.58 = 3.13
pred1 = (0.67*yu[-1] + 0.65*feng[-1] + 0.59*vogtmann[-1] + 0.64*zeller[-1] + 0.58*thomas[-1] ) / 3.13
auc(test$status, pred1[,1]) //0.664
##weights based on logistic regression 
pred2 = (16.820687*yu[-1] + 5.582355*feng[-1]) / 22.403042
auc(test$status, pred2[,1]) //0.6839

##########test: feng
hannigan = read_csv("/home/rcf-40/yilingao/staging/pred/test_feng/pred_training_hannigan.csv")
yu = read_csv("/home/rcf-40/yilingao/staging/pred/test_feng/pred_training_yu.csv")
vogtmann = read_csv("/home/rcf-40/yilingao/staging/pred/test_feng/pred_training_vogtmann.csv")
zeller = read_csv("/home/rcf-40/yilingao/staging/pred/test_feng/pred_training_zeller.csv")
thomas = read_csv("/home/rcf-40/yilingao/staging/pred/test_feng/pred_training_thomas.csv")

test = known_feng

##average
pred = (hannigan[-1] + yu[-1] + zeller[-1] + thomas[-1] + vogtmann[-1]) / 5
auc(test$status, pred[,1]) //0.893
##weights based on cross-dataset AUC
0.85+0.60+0.75+0.80+0.82 = 3.82
pred1 = (0.85*yu[-1] + 0.60*hannigan[-1] + 0.75*vogtmann[-1] + 0.80*zeller[-1] + 0.82*thomas[-1] ) / 3.82
auc(test$status, pred1[,1]) //0.8954


##########test: vogtmann
hannigan = read_csv("/home/rcf-40/yilingao/staging/pred/test_vogtmann/pred_training_hannigan.csv")
yu = read_csv("/home/rcf-40/yilingao/staging/pred/test_vogtmann/pred_training_yu.csv")
feng = read_csv("/home/rcf-40/yilingao/staging/pred/test_vogtmann/pred_training_feng.csv")
zeller = read_csv("/home/rcf-40/yilingao/staging/pred/test_vogtmann/pred_training_zeller.csv")
thomas = read_csv("/home/rcf-40/yilingao/staging/pred/test_vogtmann/pred_training_thomas.csv")

test = known_vogtmann

##average
pred = (hannigan[-1] + yu[-1] + zeller[-1] + thomas[-1] + feng[-1]) / 5
auc(test$status, pred[,1]) //0.7777
##weights based on cross-dataset AUC
0.76+0.55+0.73+0.66+0.72 = 3.42
pred1 = (0.76*yu[-1] + 0.55*hannigan[-1] + 0.73*feng[-1] + 0.66*zeller[-1] + 0.72*thomas[-1] ) / 3.42
auc(test$status, pred1[,1]) //0.7807


##########test: zeller
hannigan = read_csv("/home/rcf-40/yilingao/staging/pred/test_zeller/pred_training_hannigan.csv")
yu = read_csv("/home/rcf-40/yilingao/staging/pred/test_zeller/pred_training_yu.csv")
feng = read_csv("/home/rcf-40/yilingao/staging/pred/test_zeller/pred_training_feng.csv")
vogtmann = read_csv("/home/rcf-40/yilingao/staging/pred/test_zeller/pred_training_vogtmann.csv")
thomas = read_csv("/home/rcf-40/yilingao/staging/pred/test_zeller/pred_training_thomas.csv")

test = known_zeller

##average
pred = (hannigan[-1] + yu[-1] + vogtmann[-1] + thomas[-1] + feng[-1]) / 5
auc(test$status, pred[,1]) //0.7806
##weights based on cross-dataset AUC
0.79+0.67+0.73+0.70+0.69 = 3.58
pred1 = (0.79*yu[-1] + 0.67*hannigan[-1] + 0.73*feng[-1] + 0.70*vogtmann[-1] + 0.69*thomas[-1] ) / 3.58
auc(test$status, pred1[,1]) //0.7801


##########test: thomas
hannigan = read_csv("/home/rcf-40/yilingao/staging/pred/test_thomas/pred_training_hannigan.csv")
yu = read_csv("/home/rcf-40/yilingao/staging/pred/test_thomas/pred_training_yu.csv")
feng = read_csv("/home/rcf-40/yilingao/staging/pred/test_thomas/pred_training_feng.csv")
vogtmann = read_csv("/home/rcf-40/yilingao/staging/pred/test_thomas/pred_training_vogtmann.csv")
zeller = read_csv("/home/rcf-40/yilingao/staging/pred/test_thomas/pred_training_zeller.csv")

test = known_thomas

##average
pred = (hannigan[-1] + yu[-1] + vogtmann[-1] + zeller[-1] + feng[-1]) / 5
auc(test$status, pred[,1]) //0.7749
##weights based on cross-dataset AUC
0.74+0.55+0.77+0.74+0.75 = 3.55
pred1 = (0.74*yu[-1] + 0.55*hannigan[-1] + 0.77*feng[-1] + 0.74*vogtmann[-1] + 0.75*zeller[-1] ) / 3.55
auc(test$status, pred1[,1]) //0.773

