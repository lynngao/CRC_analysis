library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(curatedMetagenomicData)

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


model_stacking = as.data.frame(matrix(NA, ncol = 6, nrow = 30))
colnames(model_stacking) = c("yu", "hannigan", "feng", "vogtmann", "zeller", "thomas")
for(i in 1:30) {
	pred_yu = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_yu",i,".csv")))
	rownames(pred_yu) = pred_yu[,1]
	pred_yu[,1] = NULL
	for (j in 1:nrow(pred_yu)){
		known_coef = auc(known_yu$status[rownames(known_yu)%in%rownames(pred_yu)[-j]], pred_yu$known_CRC[-j])
		unknown_coef = auc(known_yu$status[rownames(known_yu)%in%rownames(pred_yu)[-j]], pred_yu$unknown_CRC[-j])
		known_coef = ifelse(known_coef < 0.5, 0.5, known_coef)
		unknown_coef = ifelse(unknown_coef < 0.5, 0.5, unknown_coef)
		pred_yu$avg_control[j]= ((known_coef-0.5)*pred_yu[j,1] + (unknown_coef-0.5)*pred_yu[j,3])/(known_coef+unknown_coef-1)
		pred_yu$avg_CRC[j]= ((known_coef-0.5)*pred_yu[j,2] + (unknown_coef-0.5)*pred_yu[j,4])/(known_coef+unknown_coef-1)
	}
	model_stacking[i,1] = auc(known_yu$status[rownames(known_yu)%in%rownames(pred_yu)], pred_yu$avg_CRC)

	pred_hannigan = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_hannigan",i,".csv")))
	rownames(pred_hannigan) = pred_hannigan[,1]
	pred_hannigan[,1] = NULL
	for (j in 1:nrow(pred_hannigan)){
		known_coef = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(pred_hannigan)[-j]], pred_hannigan$known_CRC[-j])
		unknown_coef = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(pred_hannigan)[-j]], pred_hannigan$unknown_CRC[-j])
		known_coef = ifelse(known_coef < 0.5, 0.5, known_coef)
		unknown_coef = ifelse(unknown_coef < 0.5, 0.5, unknown_coef)
		pred_hannigan$avg_control[j]= ((known_coef-0.5)*pred_hannigan[j,1] + (unknown_coef-0.5)*pred_hannigan[j,3])/(known_coef+unknown_coef-1)
		pred_hannigan$avg_CRC[j]= ((known_coef-0.5)*pred_hannigan[j,2] + (unknown_coef-0.5)*pred_hannigan[j,4])/(known_coef+unknown_coef-1)
	}
	model_stacking[i,2] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(pred_hannigan)], pred_hannigan$avg_CRC)

	pred_feng = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_feng",i,".csv")))
	rownames(pred_feng) = pred_feng[,1]
	pred_feng[,1] = NULL
	for (j in 1:nrow(pred_feng)){
		known_coef = auc(known_feng$status[rownames(known_feng)%in%rownames(pred_feng)[-j]], pred_feng$known_CRC[-j])
		unknown_coef = auc(known_feng$status[rownames(known_feng)%in%rownames(pred_feng)[-j]], pred_feng$unknown_CRC[-j])
		known_coef = ifelse(known_coef < 0.5, 0.5, known_coef)
		unknown_coef = ifelse(unknown_coef < 0.5, 0.5, unknown_coef)
		pred_feng$avg_control[j]= ((known_coef-0.5)*pred_feng[j,1] + (unknown_coef-0.5)*pred_feng[j,3])/(known_coef+unknown_coef-1)
		pred_feng$avg_CRC[j]= ((known_coef-0.5)*pred_feng[j,2] + (unknown_coef-0.5)*pred_feng[j,4])/(known_coef+unknown_coef-1)
	}
	model_stacking[i,3] = auc(known_feng$status[rownames(known_feng)%in%rownames(pred_feng)], pred_feng$avg_CRC)

	pred_vogtmann = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_vogtmann",i,".csv")))
	rownames(pred_vogtmann) = pred_vogtmann[,1]
	pred_vogtmann[,1] = NULL
	for (j in 1:nrow(pred_vogtmann)){
		known_coef = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(pred_vogtmann)[-j]], pred_vogtmann$known_CRC[-j])
		unknown_coef = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(pred_vogtmann)[-j]], pred_vogtmann$unknown_CRC[-j])
		known_coef = ifelse(known_coef < 0.5, 0.5, known_coef)
		unknown_coef = ifelse(unknown_coef < 0.5, 0.5, unknown_coef)
		pred_vogtmann$avg_control[j]= ((known_coef-0.5)*pred_vogtmann[j,1] + (unknown_coef-0.5)*pred_vogtmann[j,3])/(known_coef+unknown_coef-1)
		pred_vogtmann$avg_CRC[j]= ((known_coef-0.5)*pred_vogtmann[j,2] + (unknown_coef-0.5)*pred_vogtmann[j,4])/(known_coef+unknown_coef-1)
	}
	model_stacking[i,4] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(pred_vogtmann)], pred_vogtmann$avg_CRC)

	pred_zeller = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_zeller",i,".csv")))
	rownames(pred_zeller) = pred_zeller[,1]
	pred_zeller[,1] = NULL
	for (j in 1:nrow(pred_zeller)){
		known_coef = auc(known_zeller$status[rownames(known_zeller)%in%rownames(pred_zeller)[-j]], pred_zeller$known_CRC[-j])
		unknown_coef = auc(known_zeller$status[rownames(known_zeller)%in%rownames(pred_zeller)[-j]], pred_zeller$unknown_CRC[-j])
		known_coef = ifelse(known_coef < 0.5, 0.5, known_coef)
		unknown_coef = ifelse(unknown_coef < 0.5, 0.5, unknown_coef)
		pred_zeller$avg_control[j]= ((known_coef-0.5)*pred_zeller[j,1] + (unknown_coef-0.5)*pred_zeller[j,3])/(known_coef+unknown_coef-1)
		pred_zeller$avg_CRC[j]= ((known_coef-0.5)*pred_zeller[j,2] + (unknown_coef-0.5)*pred_zeller[j,4])/(known_coef+unknown_coef-1)
	}
	model_stacking[i,5] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(pred_zeller)], pred_zeller$avg_CRC)

	pred_thomas = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_thomas",i,".csv")))
	rownames(pred_thomas) = pred_thomas[,1]
	pred_thomas[,1] = NULL
	for (j in 1:nrow(pred_thomas)){
		known_coef = auc(known_thomas$status[rownames(known_thomas)%in%rownames(pred_thomas)[-j]], pred_thomas$known_CRC[-j])
		unknown_coef = auc(known_thomas$status[rownames(known_thomas)%in%rownames(pred_thomas)[-j]], pred_thomas$unknown_CRC[-j])
		known_coef = ifelse(known_coef < 0.5, 0.5, known_coef)
		unknown_coef = ifelse(unknown_coef < 0.5, 0.5, unknown_coef)
		pred_thomas$avg_control[j]= ((known_coef-0.5)*pred_thomas[j,1] + (unknown_coef-0.5)*pred_thomas[j,3])/(known_coef+unknown_coef-1)
		pred_thomas$avg_CRC[j]= ((known_coef-0.5)*pred_thomas[j,2] + (unknown_coef-0.5)*pred_thomas[j,4])/(known_coef+unknown_coef-1)
	}
	model_stacking[i,6] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(pred_thomas)], pred_thomas$avg_CRC)
}

mean(model_stacking[,1]) #0.9026263
mean(model_stacking[,2]) #0.692037
mean(model_stacking[,3]) #0.9451282
mean(model_stacking[,4]) #0.7057307
mean(model_stacking[,5]) #0.9042475
mean(model_stacking[,6]) #0.7606061

sd(model_stacking[,1]) #0.060101
sd(model_stacking[,2]) #0.1162205
sd(model_stacking[,3]) #0.04653833
sd(model_stacking[,4]) #0.09281038
sd(model_stacking[,5]) #0.04549201
sd(model_stacking[,6]) #0.1089047


##fit logistic regression for model stacking
for(i in 1:30) {
	pred_yu = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_yu",i,".csv")))
	rownames(pred_yu) = pred_yu[,1]
	pred_yu[,1] = NULL
	pred_yu$status = known_yu$status[rownames(known_yu)%in%rownames(pred_yu)]
	pred_yu$status = factor(pred_yu$status)
	mylogit <- glm(status ~ 0 + ., data = pred_yu[,-c(1,3)], family = "binomial"(link = "logit"))
	#mylogit <- glmnet(x=as.matrix(pred_yu[,c(2,4)]), y=pred_yu$status, family = "binomial", lower.limits = 0, lambda = 0, standardize=TRUE, intercept=FALSE)
    x = coef(mylogit)
    pred_yu$avg_control = x[1]*pred_yu[,1] + x[2]*pred_yu[,3]
	pred_yu$avg_CRC = x[1]*pred_yu[,2] + x[2]*pred_yu[,4]
	pred_yu$status = ifelse(pred_yu$status == "control", 0 , 1)
	model_stacking[i,1] = auc(pred_yu$status, pred_yu$avg_CRC)

	pred_hannigan = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_hannigan",i,".csv")))
	rownames(pred_hannigan) = pred_hannigan[,1]
	pred_hannigan[,1] = NULL
	pred_hannigan$status = known_hannigan$status[rownames(known_hannigan)%in%rownames(pred_hannigan)]
	pred_hannigan$status = factor(pred_hannigan$status)
	mylogit <- glm(status ~ 0 + ., data = pred_hannigan[,-c(1,3)], family = "binomial"(link = "logit"))
	#mylogit <- glmnet(x=as.matrix(pred_hannigan[,c(2,4)]), y=pred_hannigan$status, family = "binomial", lower.limits = 0, lambda = 0, standardize=TRUE, intercept=FALSE)
    x = coef(mylogit)
    pred_hannigan$avg_control = x[1]*pred_hannigan[,1] + x[2]*pred_hannigan[,3]
	pred_hannigan$avg_CRC = x[1]*pred_hannigan[,2] + x[2]*pred_hannigan[,4]
	pred_hannigan$status = ifelse(pred_hannigan$status == "control", 0 , 1)
	model_stacking[i,2] = auc(pred_hannigan$status, pred_hannigan$avg_CRC)

	pred_feng = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_feng",i,".csv")))
	rownames(pred_feng) = pred_feng[,1]
	pred_feng[,1] = NULL
	pred_feng$status = known_feng$status[rownames(known_feng)%in%rownames(pred_feng)]
	pred_feng$status = factor(pred_feng$status)
	mylogit <- glm(status ~ 0 + ., data = pred_feng[,-c(1,3)], family = "binomial"(link = "logit"))
	#mylogit <- glmnet(x=as.matrix(pred_feng[,c(2,4)]), y=pred_feng$status, family = "binomial", lower.limits = 0, lambda = 0, standardize=TRUE, intercept=FALSE)
    x = coef(mylogit)
    pred_feng$avg_control = x[1]*pred_feng[,1] + x[2]*pred_feng[,3]
	pred_feng$avg_CRC = x[1]*pred_feng[,2] + x[2]*pred_feng[,4]
	pred_feng$status = ifelse(pred_feng$status == "control", 0 , 1)
	model_stacking[i,3] = auc(pred_feng$status, pred_feng$avg_CRC)

	pred_vogtmann = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_vogtmann",i,".csv")))
	rownames(pred_vogtmann) = pred_vogtmann[,1]
	pred_vogtmann[,1] = NULL
	pred_vogtmann$status = known_vogtmann$status[rownames(known_vogtmann)%in%rownames(pred_vogtmann)]
	pred_vogtmann$status = factor(pred_vogtmann$status)
	mylogit <- glm(status ~ 0 + ., data = pred_vogtmann[,-c(1,3)], family = "binomial"(link = "logit"))
	#mylogit <- glmnet(x=as.matrix(pred_vogtmann[,c(2,4)]), y=pred_vogtmann$status, family = "binomial", lower.limits = 0, lambda = 0, standardize=TRUE, intercept=FALSE)
    x = coef(mylogit)
    pred_vogtmann$avg_control = x[1]*pred_vogtmann[,1] + x[2]*pred_vogtmann[,3]
	pred_vogtmann$avg_CRC = x[1]*pred_vogtmann[,2] + x[2]*pred_vogtmann[,4]
	pred_vogtmann$status = ifelse(pred_vogtmann$status == "control", 0 , 1)
	model_stacking[i,4] = auc(pred_vogtmann$status, pred_vogtmann$avg_CRC)

	pred_zeller = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_zeller",i,".csv")))
	rownames(pred_zeller) = pred_zeller[,1]
	pred_zeller[,1] = NULL
	pred_zeller$status = known_zeller$status[rownames(known_zeller)%in%rownames(pred_zeller)]
	pred_zeller$status = factor(pred_zeller$status)
	mylogit <- glm(status ~ 0 + ., data = pred_zeller[,-c(1,3)], family = "binomial"(link = "logit"))
	#mylogit <- glmnet(x=as.matrix(pred_zeller[,c(2,4)]), y=pred_zeller$status, family = "binomial", lower.limits = 0, lambda = 0, standardize=TRUE, intercept=FALSE)
    x = coef(mylogit)
    pred_zeller$avg_control = x[1]*pred_zeller[,1] + x[2]*pred_zeller[,3]
	pred_zeller$avg_CRC = x[1]*pred_zeller[,2] + x[2]*pred_zeller[,4]
	pred_zeller$status = ifelse(pred_zeller$status == "control", 0 , 1)
	model_stacking[i,5] = auc(pred_zeller$status, pred_zeller$avg_CRC)

	pred_thomas = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/within_dataset/pred_thomas",i,".csv")))
	rownames(pred_thomas) = pred_thomas[,1]
	pred_thomas[,1] = NULL
	pred_thomas$status = known_thomas$status[rownames(known_thomas)%in%rownames(pred_thomas)]
	pred_thomas$status = factor(pred_thomas$status)
	mylogit <- glm(status ~ 0 + ., data = pred_thomas[,-c(1,3)], family = "binomial"(link = "logit"))
	#mylogit <- glmnet(x=as.matrix(pred_thomas[,c(2,4)]), y=pred_thomas$status, family = "binomial", lower.limits = 0, lambda = 0, standardize=TRUE, intercept=FALSE)
    x = coef(mylogit)
    pred_thomas$avg_control = x[1]*pred_thomas[,1] + x[2]*pred_thomas[,3]
	pred_thomas$avg_CRC = x[1]*pred_thomas[,2] + x[2]*pred_thomas[,4]
	pred_thomas$status = ifelse(pred_thomas$status == "control", 0 , 1)
	model_stacking[i,6] = auc(pred_thomas$status, pred_thomas$avg_CRC)
}

mean(model_stacking[,1])
mean(model_stacking[,2])
mean(model_stacking[,3])
mean(model_stacking[,4])
mean(model_stacking[,5])
mean(model_stacking[,6])

sd(model_stacking[,1])
sd(model_stacking[,2])
sd(model_stacking[,3])
sd(model_stacking[,4])
sd(model_stacking[,5])
sd(model_stacking[,6])


