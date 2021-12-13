##load helper file
source("AUPRC.R")


#########################load microbial species abundance files########################################
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
#########################################################################################################


###########################example of within-dataset AUPRC###############################################
AUPRC_within = as.data.frame(matrix(NA,nrow = 4, ncol = 1))
colnames(AUPRC_within) = c("yu")
rownames(AUPRC_within) = c("known", "unknown", "LOSO", "combined")

res_yu_known = NA
res_yu_unknown = NA
res_yu_LOSO = NA
for (i in 1:30){
  ###load 30 prediction probability files###
  yu = as.data.frame(read_csv(paste0("pred_yu",i,".csv", sep="")))
  res_known_temp = unlist(auprc(known_yu, yu)[1]) #get AUPRC for known species
  res_unknown_temp = unlist(auprc(known_yu, yu)[2]) #get AUPRC for unknown species
  res_yu_known = c(res_yu_known, res_known_temp) 
  res_yu_unknown = c(res_yu_unknown, res_unknown_temp) 
  res_yu_temp_LOSO = auprc_LOSO_within(known_yu,yu) #calculate AUPRC for LOSO model stacking method
  res_yu_LOSO = c(res_yu_sample, res_yu_temp_LOSO)
}
AUPRC_within[1,1] = mean(res_yu_known[-1])
AUPRC_within[2,1]= mean(res_yu_unknown[-1])
AUPRC_within[3,1] = mean(res_yu_sample[-1])
#combined AUPRC
res_yu_combined = NA
for (i in 1:30){
  ###load 30 prediction probability files###
  yu = as.data.frame(read_csv(paste0("combined/pred_yu",i,".csv", sep="")))
  res_combined_temp = auprc_within_combined(known_yu,known_yu) #calculate AUPRC for combined species
  res_yu_combined = c(res_yu_combined, res_combined_temp) 
}
AUPRC_within[4,1] = mean(res_yu_combined[-1])
###########################################################################################################


###########################example of cross-dataset AUPRC##################################################
######load prediction probability files for "training:Yu, Testing:Hannigan" pair
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_hannigan.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_hannigan.csv"))
hannigan_combined = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/combined_pred_test_hannigan.csv"))
hannigan_LOSO = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_hannigan.csv"))

AUPRC_training_yu = as.data.frame(matrix(NA,nrow = 4, ncol = 1))
colnames(AUPRC_training_yu) = c("hannigan")
rownames(AUPRC_training_yu) = c("known", "unknown", "LOSO", "combined")

AUPRC_training_yu[1,1] = auprc(known_hannigan, hannigan_known) #AUPRC for known species
AUPRC_training_yu[2,1] = auprc(known_hannigan, hannigan_unknown) #AUPRC for unknown species
AUPRC_training_yu[3,1] = auprc(known_hannigan, hannigan_LOSO) #AUPRC for LOSO model stacking
AUPRC_training_yu[4,1] = auprc(known_hannigan, hannigan_combined) #AUPRC for combined species
###########################################################################################################


###########################example of LODO AUPRC###########################################################
AUPRC_LODO = as.data.frame(matrix(NA,nrow = 2, ncol = 1))
colnames(AUPRC_LODO) = c("yu")
rownames(AUPRC_LODO) = c("known", "LOSO")

######AUPRC for known species from LODO prediction probabilities
res_yu = NA
for (i in 1:10){
  yu = as.data.frame(read_csv(paste0("/Users/lynngao/Desktop/pred/LODO/pred_yu",i,".csv", sep="")))
  res_temp = auprc(yu,known_yu)
  res_yu = c(res_yu, res_temp)
}
AUPRC_LODO[1,1] = mean(res_yu[-1])

#######AUPRC for LOSO from cross-dataset prediction probabilities, testing on Yu dataset
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

AUPRC_LODO[1,2] = auprc_LODO(test, hannigan_known, hannigan_unknown, feng_known, feng_unknown, 
	vogtmann_known, vogtmann_unknown, zeller_known, zeller_unknown, thomas_known, thomas_unknown, return_pred=TRUE)
###########################################################################################################







