library(readr)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(glmnet)

####################################README##################################
####################require species abundance table#########################
##require prediction probabilties for calculating LOSO model stacking AUC###
############################################################################

#################################within-dataset function####################################
within_model_stacking <- function(ds, prob_ds){
	for (i in 1:nrow(prob_ds)){
		known_coef = auc(ds$status[rownames(ds)%in%rownames(prob_ds)[-i]], prob_ds$known_CRC[-i])
		unknown_coef = auc(ds$status[rownames(ds)%in%rownames(prob_ds)[-i]], prob_ds$unknown_CRC[-i])
		known_coef = ifelse(known_coef < 0.5, 0.5, known_coef)
		unknown_coef = ifelse(unknown_coef < 0.5, 0.5, unknown_coef)
		prob_ds$avg_control[i]= ((known_coef-0.5)*prob_ds[i,1] + (unknown_coef-0.5)*prob_ds[i,3])/(known_coef+unknown_coef-1)
		prob_ds$avg_CRC[i]= ((known_coef-0.5)*prob_ds[i,2] + (unknown_coef-0.5)*prob_ds[i,4])/(known_coef+unknown_coef-1)
	}
	AUC = auc(ds$status[rownames(ds)%in%rownames(prob_ds)], prob_ds$avg_CRC)
	return(AUC)
}

#################################cross-dataset function####################################
cross_model_stacking <- function(ds, prob_known, prob_unknown){
	AUC_df = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(prob_known)))
	colnames(AUC_df) = c("known_AUC","unknown_AUC")
	rownames(AUC_df) = prob_known[,1]
	for (i in 1:nrow(AUC_df)){
		AUC_df[i,1] = auc(ds$status[rownames(ds)%in%rownames(AUC_df)[-i]], prob_known$CRC[prob_known[,1]%in%rownames(AUC_df)[-i]])
		AUC_df[i,2] = auc(ds$status[rownames(ds)%in%rownames(AUC_df)[-i]], prob_unknown$CRC[prob_unknown[,1]%in%rownames(AUC_df)[-i]])
	}
	AUC_df[AUC_df < 0.5] = 0.5
	prob_df = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(prob_known)))
	colnames(prob_df) = c("avg_control","avg_CRC")
	rownames(prob_df) = prob_known[,1]
	prob_df$avg_control = ((AUC_df[,1]-0.5)*prob_known[,2] + (AUC_df[,2]-0.5)*prob_unknown[,2])/(AUC_df[,1]+AUC_df[,2]-1)
	prob_df$avg_CRC = ((AUC_df[,1]-0.5)*prob_known[,3] + (AUC_df[,2]-0.5)*prob_unknown[,3])/(AUC_df[,1]+AUC_df[,2]-1)

	AUC = auc(ds$status, prob_df$avg_CRC)
	return(AUC)
}

#################################LODO function####################################
###function for create AUC table
AUC_table <- function(test, prob_known, prob_unknown) {
	AUC_df = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(prob_ds)))
	colnames(AUC_df) = c("known_AUC","unknown_AUC")
	rownames(AUC_df) = [,1]
	for (i in 1:nrow(AUC_df)){
		AUC_df[i,1] = auc(test$status[rownames(test)%in%rownames(AUC_df)[-i]], prob_known$CRC[prob_known[,1]%in%rownames(AUC_df)[-i]])
		AUC_df[i,2] = auc(test$status[rownames(test)%in%rownames(AUC_df)[-i]], prob_unknown$CRC[prob_unknown[,1]%in%rownames(AUC_df)[-i]])
	}
	AUC_df[AUC_df < 0.5] = 0.5
	return(AUC_df)

}

LODO_model_stacking <- function(test, prob_known1, prob_unknown1, prob_known2, prob_unknown2, 
	prob_known3, prob_unknown3,prob_known4, prob_unknown4, prob_known5, prob_unknown5) {
	AUC1 <- AUC_table(test, prob_known1, prob_unknown1)
	AUC2 <- AUC_table(test, prob_known2, prob_unknown2)
	AUC3 <- AUC_table(test, prob_known3, prob_unknown3)
	AUC4 <- AUC_table(test, prob_known4, prob_unknown4)
	AUC5 <- AUC_table(test, prob_known5, prob_unknown5)

	pred_known = AUC1[,1]-0.5)*prob_known1[,3] + (AUC2[,1]-0.5)*prob_known2[,3] + (AUC3[,1]-0.5)*prob_known3[,3] + (AUC4[,1]-0.5)*prob_known4[,3] + (AUC5[,1]-0.5)*prob_known5[,3]
	pred_unknown = AUC1[,2]-0.5)*prob_unknown1[,3] + (AUC2[,2]-0.5)*prob_unknown2[,3] + (AUC3[,2]-0.5)*prob_unknown3[,3] + (AUC4[,2]-0.5)*prob_unknown4[,3] + (AUC5[,2]-0.5)*prob_unknown5[,3]
	pred = (pred_known + pred_unknown) / (AUC1[,1] + AUC2[,1] + AUC3[,1] + AUC4[,1] + AUC5[,1] + AUC1[,2] + AUC2[,2] + AUC3[,2] + AUC4[,2] + AUC5[,2]- 5)

	AUC = auc(test$status, pred)
	return(AUC)
}

















