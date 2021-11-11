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
  dataset[,3] = round(dataset[,3], digits = 2)
  dataset[,2] = 1-dataset[,3]
  for (i in 1:nrow(dataset)) {
    if (dataset[i,1]%in%rownames(ds)) {
      dataset$status[i] = ds$status[rownames(ds) == dataset[i,1]]
    }
  }
  dataset$status = ifelse(dataset$status == "control", 0, 1)
  pr.curve(scores.class0 = dataset[,3], weights.class0 = dataset$status, curve = T)
  
}


#################################training:yu############################
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_feng.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_feng.csv"))
feng_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/combined_pred_test_feng.csv"))
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_hannigan.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_hannigan.csv"))
hannigan_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/combined_pred_test_hannigan.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_vogtmann.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_vogtmann.csv"))
vogtmann_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/combined_pred_test_vogtmann.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_zeller.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_zeller.csv"))
zeller_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/combined_pred_test_zeller.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_thomas.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_thomas.csv"))
thomas_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/combined_pred_test_thomas.csv"))

hannigan_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_hannigan.csv"))
feng_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_feng.csv"))
vogtmann_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_vogtmann.csv"))
zeller_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_zeller.csv"))
thomas_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_thomas.csv"))


AUPRC_training_yu = as.data.frame(matrix(NA,nrow = 4, ncol = 5))
colnames(AUPRC_training_yu) = c("hannigan", "feng", "vogtmann","zeller", "thomas")
rownames(AUPRC_training_yu) = c("known", "unknown", "LOSO", "combined")

################AUPRC for known speicies
AUPRC_training_yu[1,1] = auprc(hannigan_known,known_hannigan)$auc.integral
AUPRC_training_yu[1,2] = auprc(feng_known,known_feng)$auc.integral
AUPRC_training_yu[1,3] = auprc(vogtmann_known,known_vogtmann)$auc.integral
AUPRC_training_yu[1,4] = auprc(zeller_known,known_zeller)$auc.integral
AUPRC_training_yu[1,5] = auprc(thomas_known,known_thomas)$auc.integral

###############AUPRC for unknown speicies
AUPRC_training_yu[2,1] = auprc(hannigan_unknown,known_hannigan)$auc.integral
AUPRC_training_yu[2,2] = auprc(feng_unknown,known_feng)$auc.integral
AUPRC_training_yu[2,3] = auprc(vogtmann_unknown,known_vogtmann)$auc.integral
AUPRC_training_yu[2,4] = auprc(zeller_unknown,known_zeller)$auc.integral
AUPRC_training_yu[2,5] = auprc(thomas_unknown,known_thomas)$auc.integral

###############AUPRC for combined speicies
AUPRC_training_yu[4,1] = auprc(hannigan_combined,known_hannigan)$auc.integral
AUPRC_training_yu[4,2] = auprc(feng_combined,known_feng)$auc.integral
AUPRC_training_yu[4,3] = auprc(vogtmann_combined,known_vogtmann)$auc.integral
AUPRC_training_yu[4,4] = auprc(zeller_combined,known_zeller)$auc.integral
AUPRC_training_yu[4,5] = auprc(thomas_combined,known_thomas)$auc.integral


##############AUPRC for LOSO
AUPRC_training_yu[3,1] = auprc(hannigan_LOSO,known_hannigan)$auc.integral
AUPRC_training_yu[3,2] = auprc(feng_LOSO,known_feng)$auc.integral
AUPRC_training_yu[3,3] = auprc(vogtmann_LOSO,known_vogtmann)$auc.integral
AUPRC_training_yu[3,4] = auprc(zeller_LOSO,known_zeller)$auc.integral
AUPRC_training_yu[3,5] = auprc(thomas_LOSO,known_thomas)$auc.integral

AUPRC_training_yu = round(AUPRC_training_yu, digits = 2)
write.csv(AUPRC_training_yu, "/Users/lynngao/Desktop/results/AUPRC/cross_dataset/AUPRC_training_yu1.csv")
##########################################################################################

################################training:hannigan########################################
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_feng.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_feng.csv"))
feng_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/combined_pred_test_feng.csv"))
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_yu.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_yu.csv"))
yu_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/combined_pred_test_yu.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_vogtmann.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_vogtmann.csv"))
vogtmann_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/combined_pred_test_vogtmann.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_zeller.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_zeller.csv"))
zeller_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/combined_pred_test_zeller.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/known_pred_test_thomas.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/unknown_pred_test_thomas.csv"))
thomas_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/combined_pred_test_thomas.csv"))

yu_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_yu.csv"))
feng_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_feng.csv"))
vogtmann_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_vogtmann.csv"))
zeller_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_zeller.csv"))
thomas_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_thomas.csv"))

AUPRC_training_hannigan = as.data.frame(matrix(NA,nrow = 4, ncol = 5))
colnames(AUPRC_training_hannigan) = c("yu", "feng", "vogtmann","zeller", "thomas")
rownames(AUPRC_training_hannigan) = c("known", "unknown", "LOSO", "combined")

################AUPRC for known speicies
AUPRC_training_hannigan[1,1] = auprc(yu_known,known_yu)$auc.integral
AUPRC_training_hannigan[1,2] = auprc(feng_known,known_feng)$auc.integral
AUPRC_training_hannigan[1,3] = auprc(vogtmann_known,known_vogtmann)$auc.integral
AUPRC_training_hannigan[1,4] = auprc(zeller_known,known_zeller)$auc.integral
AUPRC_training_hannigan[1,5] = auprc(thomas_known,known_thomas)$auc.integral

###############AUPRC for unknown speicies
AUPRC_training_hannigan[2,1] = auprc(yu_unknown,known_yu)$auc.integral
AUPRC_training_hannigan[2,2] = auprc(feng_unknown,known_feng)$auc.integral
AUPRC_training_hannigan[2,3] = auprc(vogtmann_unknown,known_vogtmann)$auc.integral
AUPRC_training_hannigan[2,4] = auprc(zeller_unknown,known_zeller)$auc.integral
AUPRC_training_hannigan[2,5] = auprc(thomas_unknown,known_thomas)$auc.integral

###############AUPRC for combined speicies
AUPRC_training_hannigan[4,1] = auprc(yu_combined,known_yu)$auc.integral
AUPRC_training_hannigan[4,2] = auprc(feng_combined,known_feng)$auc.integral
AUPRC_training_hannigan[4,3] = auprc(vogtmann_combined,known_vogtmann)$auc.integral
AUPRC_training_hannigan[4,4] = auprc(zeller_combined,known_zeller)$auc.integral
AUPRC_training_hannigan[4,5] = auprc(thomas_combined,known_thomas)$auc.integral


##############AUPRC for LOSO
AUPRC_training_hannigan[3,1] = auprc(yu_LOSO,known_yu)$auc.integral
AUPRC_training_hannigan[3,2] = auprc(feng_LOSO,known_feng)$auc.integral
AUPRC_training_hannigan[3,3] = auprc(vogtmann_LOSO,known_vogtmann)$auc.integral
AUPRC_training_hannigan[3,4] = auprc(zeller_LOSO,known_zeller)$auc.integral
AUPRC_training_hannigan[3,5] = auprc(thomas_LOSO,known_thomas)$auc.integral

AUPRC_training_hannigan = round(AUPRC_training_hannigan, digits = 2)
write.csv(AUPRC_training_hannigan, "/Users/lynngao/Desktop/results/AUPRC/cross_dataset/AUPRC_training_hannigan1.csv")
##################################################################################################


##########training:feng################
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_hannigan.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_hannigan.csv"))
hannigan_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/combined_pred_test_hannigan.csv"))
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_yu.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_yu.csv"))
yu_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/combined_pred_test_yu.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_vogtmann.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_vogtmann.csv"))
vogtmann_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/combined_pred_test_vogtmann.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_zeller.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_zeller.csv"))
zeller_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/combined_pred_test_zeller.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/known_pred_test_thomas.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/unknown_pred_test_thomas.csv"))
thomas_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/combined_pred_test_thomas.csv"))

yu_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_yu.csv"))
hannigan_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_hannigan.csv"))
vogtmann_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_vogtmann.csv"))
zeller_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_zeller.csv"))
thomas_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_thomas.csv"))

AUPRC_training_feng = as.data.frame(matrix(NA,nrow = 4, ncol = 5))
colnames(AUPRC_training_feng) = c("yu", "hannigan", "vogtmann","zeller", "thomas")
rownames(AUPRC_training_feng) = c("known", "unknown", "LOSO", "combined")

################AUPRC for known speicies
AUPRC_training_feng[1,1] = auprc(yu_known,known_yu)$auc.integral
AUPRC_training_feng[1,2] = auprc(hannigan_known,known_hannigan)$auc.integral
AUPRC_training_feng[1,3] = auprc(vogtmann_known,known_vogtmann)$auc.integral
AUPRC_training_feng[1,4] = auprc(zeller_known,known_zeller)$auc.integral
AUPRC_training_feng[1,5] = auprc(thomas_known,known_thomas)$auc.integral

###############AUPRC for unknown speicies
AUPRC_training_feng[2,1] = auprc(yu_unknown,known_yu)$auc.integral
AUPRC_training_feng[2,2] = auprc(hannigan_unknown,known_hannigan)$auc.integral
AUPRC_training_feng[2,3] = auprc(vogtmann_unknown,known_vogtmann)$auc.integral
AUPRC_training_feng[2,4] = auprc(zeller_unknown,known_zeller)$auc.integral
AUPRC_training_feng[2,5] = auprc(thomas_unknown,known_thomas)$auc.integral

###############AUPRC for combined speicies
AUPRC_training_feng[4,1] = auprc(yu_combined,known_yu)$auc.integral
AUPRC_training_feng[4,2] = auprc(hannigan_combined,known_hannigan)$auc.integral
AUPRC_training_feng[4,3] = auprc(vogtmann_combined,known_vogtmann)$auc.integral
AUPRC_training_feng[4,4] = auprc(zeller_combined,known_zeller)$auc.integral
AUPRC_training_feng[4,5] = auprc(thomas_combined,known_thomas)$auc.integral

##############AUPRC for LOSO
AUPRC_training_feng[3,1] = auprc(yu_LOSO,known_yu)$auc.integral
AUPRC_training_feng[3,2] = auprc(hannigan_LOSO,known_hannigan)$auc.integral
AUPRC_training_feng[3,3] = auprc(vogtmann_LOSO,known_vogtmann)$auc.integral
AUPRC_training_feng[3,4] = auprc(zeller_LOSO,known_zeller)$auc.integral
AUPRC_training_feng[3,5] = auprc(thomas_LOSO,known_thomas)$auc.integral

AUPRC_training_feng = round(AUPRC_training_feng, digits = 2)
write.csv(AUPRC_training_feng, "/Users/lynngao/Desktop/results/AUPRC/cross_dataset/AUPRC_training_feng1.csv")
##################################################################################################


##########training:vogtmann#################
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_hannigan.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_hannigan.csv"))
hannigan_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/combined_pred_test_hannigan.csv"))
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_yu.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_yu.csv"))
yu_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/combined_pred_test_yu.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_feng.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_feng.csv"))
feng_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/combined_pred_test_feng.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_zeller.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_zeller.csv"))
zeller_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/combined_pred_test_zeller.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/known_pred_test_thomas.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/unknown_pred_test_thomas.csv"))
thomas_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/combined_pred_test_thomas.csv"))

yu_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_yu.csv"))
hannigan_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_hannigan.csv"))
feng_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_feng.csv"))
zeller_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_zeller.csv"))
thomas_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_thomas.csv"))

AUPRC_training_vogtmann = as.data.frame(matrix(NA,nrow = 4, ncol = 5))
colnames(AUPRC_training_vogtmann) = c("yu", "hannigan", "feng","zeller", "thomas")
rownames(AUPRC_training_vogtmann) = c("known", "unknown", "LOSO", "combined")

################AUPRC for known speicies
AUPRC_training_vogtmann[1,1] = auprc(yu_known,known_yu)$auc.integral
AUPRC_training_vogtmann[1,2] = auprc(hannigan_known,known_hannigan)$auc.integral
AUPRC_training_vogtmann[1,3] = auprc(feng_known,known_feng)$auc.integral
AUPRC_training_vogtmann[1,4] = auprc(zeller_known,known_zeller)$auc.integral
AUPRC_training_vogtmann[1,5] = auprc(thomas_known,known_thomas)$auc.integral

###############AUPRC for unknown speicies
AUPRC_training_vogtmann[2,1] = auprc(yu_unknown,known_yu)$auc.integral
AUPRC_training_vogtmann[2,2] = auprc(hannigan_unknown,known_hannigan)$auc.integral
AUPRC_training_vogtmann[2,3] = auprc(feng_unknown,known_feng)$auc.integral
AUPRC_training_vogtmann[2,4] = auprc(zeller_unknown,known_zeller)$auc.integral
AUPRC_training_vogtmann[2,5] = auprc(thomas_unknown,known_thomas)$auc.integral

###############AUPRC for unknown speicies
AUPRC_training_vogtmann[4,1] = auprc(yu_combined,known_yu)$auc.integral
AUPRC_training_vogtmann[4,2] = auprc(hannigan_combined,known_hannigan)$auc.integral
AUPRC_training_vogtmann[4,3] = auprc(feng_combined,known_feng)$auc.integral
AUPRC_training_vogtmann[4,4] = auprc(zeller_combined,known_zeller)$auc.integral
AUPRC_training_vogtmann[4,5] = auprc(thomas_combined,known_thomas)$auc.integral

##############AUPRC for LOSO
AUPRC_training_vogtmann[3,1] = auprc(yu_LOSO,known_yu)$auc.integral
AUPRC_training_vogtmann[3,2] = auprc(hannigan_LOSO,known_hannigan)$auc.integral
AUPRC_training_vogtmann[3,3] = auprc(feng_LOSO,known_feng)$auc.integral
AUPRC_training_vogtmann[3,4] = auprc(zeller_LOSO,known_zeller)$auc.integral
AUPRC_training_vogtmann[3,5] = auprc(thomas_LOSO,known_thomas)$auc.integral

AUPRC_training_vogtmann = round(AUPRC_training_vogtmann, digits = 2)
write.csv(AUPRC_training_vogtmann, "/Users/lynngao/Desktop/results/AUPRC/cross_dataset/AUPRC_training_vogtmann1.csv")
##################################################################################################



##########training:zeller######################
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_hannigan.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_hannigan.csv"))
hannigan_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/combined_pred_test_hannigan.csv"))
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_yu.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_yu.csv"))
yu_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/combined_pred_test_yu.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_feng.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_feng.csv"))
feng_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/combined_pred_test_feng.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_vogtmann.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_vogtmann.csv"))
vogtmann_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/combined_pred_test_vogtmann.csv"))
thomas_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/known_pred_test_thomas.csv"))
thomas_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/unknown_pred_test_thomas.csv"))
thomas_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/combined_pred_test_thomas.csv"))

yu_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_yu.csv"))
hannigan_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_hannigan.csv"))
feng_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_feng.csv"))
vogtmann_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_vogtmann.csv"))
thomas_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_thomas.csv"))

AUPRC_training_zeller = as.data.frame(matrix(NA,nrow = 4, ncol = 5))
colnames(AUPRC_training_zeller) = c("yu", "hannigan", "feng","vogtmann", "thomas")
rownames(AUPRC_training_zeller) = c("known", "unknown", "LOSO", "combined")

################AUPRC for known speicies
AUPRC_training_zeller[1,1] = auprc(yu_known,known_yu)$auc.integral
AUPRC_training_zeller[1,2] = auprc(hannigan_known,known_hannigan)$auc.integral
AUPRC_training_zeller[1,3] = auprc(feng_known,known_feng)$auc.integral
AUPRC_training_zeller[1,4] = auprc(vogtmann_known,known_vogtmann)$auc.integral
AUPRC_training_zeller[1,5] = auprc(thomas_known,known_thomas)$auc.integral

###############AUPRC for unknown speicies
AUPRC_training_zeller[2,1] = auprc(yu_unknown,known_yu)$auc.integral
AUPRC_training_zeller[2,2] = auprc(hannigan_unknown,known_hannigan)$auc.integral
AUPRC_training_zeller[2,3] = auprc(feng_unknown,known_feng)$auc.integral
AUPRC_training_zeller[2,4] = auprc(vogtmann_unknown,known_vogtmann)$auc.integral
AUPRC_training_zeller[2,5] = auprc(thomas_unknown,known_thomas)$auc.integral

###############AUPRC for combined speicies
AUPRC_training_zeller[4,1] = auprc(yu_combined,known_yu)$auc.integral
AUPRC_training_zeller[4,2] = auprc(hannigan_combined,known_hannigan)$auc.integral
AUPRC_training_zeller[4,3] = auprc(feng_combined,known_feng)$auc.integral
AUPRC_training_zeller[4,4] = auprc(vogtmann_combined,known_vogtmann)$auc.integral
AUPRC_training_zeller[4,5] = auprc(thomas_combined,known_thomas)$auc.integral

##############AUPRC for LOSO
AUPRC_training_zeller[3,1] = auprc(yu_LOSO,known_yu)$auc.integral
AUPRC_training_zeller[3,2] = auprc(hannigan_LOSO,known_hannigan)$auc.integral
AUPRC_training_zeller[3,3] = auprc(feng_LOSO,known_feng)$auc.integral
AUPRC_training_zeller[3,4] = auprc(vogtmann_LOSO,known_vogtmann)$auc.integral
AUPRC_training_zeller[3,5] = auprc(thomas_LOSO,known_thomas)$auc.integral

AUPRC_training_zeller = round(AUPRC_training_zeller, digits = 2)
write.csv(AUPRC_training_zeller, "/Users/lynngao/Desktop/results/AUPRC/cross_dataset/AUPRC_training_zeller1.csv")
##################################################################################################


##########training:thomas##########
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_hannigan.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_hannigan.csv"))
hannigan_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/combined_pred_test_hannigan.csv"))
yu_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_yu.csv"))
yu_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_yu.csv"))
yu_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/combined_pred_test_yu.csv"))
feng_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_feng.csv"))
feng_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_feng.csv"))
feng_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/combined_pred_test_feng.csv"))
vogtmann_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_vogtmann.csv"))
vogtmann_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_vogtmann.csv"))
vogtmann_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/combined_pred_test_vogtmann.csv"))
zeller_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/known_pred_test_zeller.csv"))
zeller_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/unknown_pred_test_zeller.csv"))
zeller_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/combined_pred_test_zeller.csv"))

yu_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_yu.csv"))
hannigan_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_hannigan.csv"))
feng_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_feng.csv"))
vogtmann_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_vogtmann.csv"))
zeller_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_zeller.csv"))

AUPRC_training_thomas = as.data.frame(matrix(NA,nrow = 4, ncol = 5))
colnames(AUPRC_training_thomas) = c("yu", "hannigan", "feng","vogtmann", "zeller")
rownames(AUPRC_training_thomas) = c("known", "unknown", "LOSO", "combined")

################AUPRC for known speicies
AUPRC_training_thomas[1,1] = auprc(yu_known,known_yu)$auc.integral
AUPRC_training_thomas[1,2] = auprc(hannigan_known,known_hannigan)$auc.integral
AUPRC_training_thomas[1,3] = auprc(feng_known,known_feng)$auc.integral
AUPRC_training_thomas[1,4] = auprc(vogtmann_known,known_vogtmann)$auc.integral
AUPRC_training_thomas[1,5] = auprc(zeller_known,known_zeller)$auc.integral

###############AUPRC for unknown speicies
AUPRC_training_thomas[2,1] = auprc(yu_unknown,known_yu)$auc.integral
AUPRC_training_thomas[2,2] = auprc(hannigan_unknown,known_hannigan)$auc.integral
AUPRC_training_thomas[2,3] = auprc(feng_unknown,known_feng)$auc.integral
AUPRC_training_thomas[2,4] = auprc(vogtmann_unknown,known_vogtmann)$auc.integral
AUPRC_training_thomas[2,5] = auprc(zeller_unknown,known_zeller)$auc.integral

###############AUPRC for combined seicies
AUPRC_training_thomas[4,1] = auprc(yu_combined,known_yu)$auc.integral
AUPRC_training_thomas[4,2] = auprc(hannigan_combined,known_hannigan)$auc.integral
AUPRC_training_thomas[4,3] = auprc(feng_combined,known_feng)$auc.integral
AUPRC_training_thomas[4,4] = auprc(vogtmann_combined,known_vogtmann)$auc.integral
AUPRC_training_thomas[4,5] = auprc(zeller_combined,known_zeller)$auc.integral

##############AUPRC for LOSO
AUPRC_training_thomas[3,1] = auprc(yu_LOSO,known_yu)$auc.integral
AUPRC_training_thomas[3,2] = auprc(hannigan_LOSO,known_hannigan)$auc.integral
AUPRC_training_thomas[3,3] = auprc(feng_LOSO,known_feng)$auc.integral
AUPRC_training_thomas[3,4] = auprc(vogtmann_LOSO,known_vogtmann)$auc.integral
AUPRC_training_thomas[3,5] = auprc(zeller_LOSO,known_zeller)$auc.integral

AUPRC_training_thomas = round(AUPRC_training_thomas, digits = 2)
write.csv(AUPRC_training_thomas, "/Users/lynngao/Desktop/results/AUPRC/cross_dataset/AUPRC_training_thomas1.csv")
##################################################################################################




