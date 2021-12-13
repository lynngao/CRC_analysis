library(readr)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(glmnet)
library(PRROC)

####################################README##################################
####################require species abundance table#########################
##require prediction probabilties for calculating LOSO model stacking AUC###
############################################################################

source("LOSO_model_stacking.R")

#################################within-dataset functions####################################
auprc_within <- function(ds, prob_ds) {
  for (i in 1:nrow(prob_ds)) {
    if (prob_ds[i,1]%in%rownames(ds)) {
      prob_ds$status[i] = ds$status[rownames(ds) == prob_ds[i,1]]
    }
  }
  prob_ds$status = ifelse(prob_ds$status == "control", 0, 1)
  auprc_known = pr.curve(scores.class0 = prob_ds[,"known_CRC"], weights.class0 = prob_ds$status, curve = T)$auc.integral
  auprc_unknown = pr.curve(scores.class0 = prob_ds[,"unknown_CRC"], weights.class0 = prob_ds$status, curve = T)$auc.integral
  return(list(auprc_known, auprc_unknown))
}

auprc_within_combined <- function(ds, prob_ds) {
  for (i in 1:nrow(prob_ds)) {
    if (prob_ds[i,1]%in%rownames(ds)) {
      prob_ds$status[i] = ds$status[rownames(ds) == prob_ds[i,1]]
    }
  }
  prob_ds$status = ifelse(prob_ds$status == "control", 0, 1)
  auprc_combined = pr.curve(scores.class0 = prob_ds[,"CRC"], weights.class0 = prob_ds$status, curve = T)$auc.integral
  return(auprc_combined)
}

auprc_LOSO_within <- function(ds, prob_ds){
  for (j in 1:nrow(prob_ds)) {
    if (prob_ds[j,1]%in%rownames(ds)) {
      prob_ds$status[j] = ds$status[rownames(ds) == prob_ds[j,1]]
    }
  }
  prob_ds$status = ifelse(prob_ds$status == "control", 0, 1)
  temp_prob = as.data.frame(matrix(NA, ncol = 4, nrow = nrow(prob_ds)))
  colnames(temp_prob) = c("known_AUC","unknown_AUC", "avg_control","avg_CRC")
  rownames(temp_prob) = prob_ds[,1]
  for (k in 1:nrow(temp_prob)){
    temp_prob[k,1] = auc(prob_ds$status[prob_ds[,1] != rownames(ds)[k]],prob_ds[prob_ds[,1] != rownames(ds)[k],"known_CRC"])
    temp_prob[k,2] = auc(prob_ds$status[prob_ds[,1] != rownames(ds)[k]],prob_ds[prob_ds[,1] != rownames(ds)[k],"unknown_CRC"])
  }
  temp_prob[temp_prob < 0.5] = 0.5
  temp_prob$avg_control = ((temp_prob[,1]-0.5)*prob_ds[,"known_control"] + (temp_prob[,2]-0.5)*prob_ds[,"unknown_control"])/(temp_prob[,1]+temp_prob[,2]-2*0.5)
  temp_prob$avg_CRC = ((temp_prob[,1]-0.5)*prob_ds[,"known_CRC"] + (temp_prob[,2]-0.5)*prob_ds[,"unknown_CRC"])/(temp_prob[,1]+temp_prob[,2]-2*0.5)
  temp_prob$status = prob_ds$status
  pr.curve(scores.class0 = temp_prob[,"avg_CRC"], weights.class0 = temp_prob$status, curve = T)$auc.integral
}


#################################cross-dataset function####################################
auprc <- function(ds, prob_ds) {
  for (i in 1:nrow(prob_ds)) {
    if (prob_ds[i,1]%in%rownames(ds)) {
      prob_ds$status[i] = ds$status[rownames(ds) == prob_ds[i,1]]
    }
  }
  prob_ds$status = ifelse(prob_ds$status == "control", 0, 1)
  pr.curve(scores.class0 = prob_ds[,3], weights.class0 = prob_ds$status, curve = T)$auc.integral
  
}

#################################LODO function####################################
auprc_LODO <- function(test, prob_known1, prob_unknown1, prob_known2, prob_unknown2, 
	prob_known3, prob_unknown3,prob_known4, prob_unknown4, prob_known5, prob_unknown5, return_pred=TRUE){
	pred = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(test)))
	rownames(pred) = rownames(test)
	colnames(pred) = c("CRC", "status")
	pred$CRC = LODO_model_stacking(test, prob_known1, prob_unknown1, prob_known2, prob_unknown2, 
		prob_known3, prob_unknown3,prob_known4, prob_unknown4, prob_known5, prob_unknown5, return_pred=TRUE)
	pred$status = test$status
	AUPRC = pr.curve(scores.class0 = pred[pred$status == "CRC", 1], scores.class1 = pred[pred$status == "control", 1], curve = T)$auc.integral
	return(AUPRC)
}
