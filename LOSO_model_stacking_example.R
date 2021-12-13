##load helper file
source("LOSO_model_stacking.R")

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


###########################example of within-dataset LOSO model stacking#################################
model_stacking = as.data.frame(matrix(NA, ncol = 1, nrow = 30))
colnames(model_stacking) = c("yu")

for(i in 1:30) {
	###load 30 prediction probability files###
	pred_yu = as.data.frame(read_csv(paste0("within_dataset/pred_yu",i,".csv")))
	rownames(pred_yu) = pred_yu[,1]
	pred_yu[,1] = NULL
	model_stacking[i,1] = within_model_stacking(known_yu, pred_yu)
}
#########################################################################################################

###########################example of cross-dataset LOSO model stacking##################################
######load prediction probability files for "training:Yu, Testing:Hannigan" pair
hannigan_known = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/known_pred_test_hannigan.csv"))
hannigan_unknown = as.data.frame(read_csv("/Users/lynngao/Desktop/pred/cross_dataset/yu/unknown_pred_test_hannigan.csv"))

cross_model_stacking(known_hannigan, hannigan_known, hannigan_unknown)
cross_model_stacking(known_hannigan, hannigan_known, hannigan_unknown, TRUE) #save LOSO probabilities for AUPRC calculation
#########################################################################################################

###########################example of LODO LOSO model stacking###########################################
#####load files for cross-dataset prediction probabilities, testing on Yu dataset
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

LODO_model_stacking(test, hannigan_known, hannigan_unknown, feng_known, feng_unknown, 
	vogtmann_known, vogtmann_unknown, zeller_known, zeller_unknown, thomas_known, thomas_unknown)
###########################################################################################################

