##load helper file
source("ComBat.R")

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

################################ComBat croo-dataset example##############################################
##create dataframe to store results
auc_cross_cohort = as.data.frame(matrix(NA,ncol = 5, nrow = 32))
colnames(auc_cross_cohort) = c("training:yu, test:hannigan", "training:yu, test:feng", "training:yu, test:vogtmann","training:yu, test:zeller", "training:yu, test:thomas")
rownames(auc_cross_cohort) = c(1:30,"mean","sd")

##training on Yu dataset, the other five datasets are severed as testing in turn
training = known_yu[1:(ncol(known_yu)-1)]
training_status = known_yu$status
test1 = known_hannigan[1:(ncol(known_hannigan)-1)]
test1_status = known_hannigan$status
test2 = known_feng[1:(ncol(known_feng)-1)]
test2_status = known_feng$status
test3 = known_vogtmann[1:(ncol(known_vogtmann)-1)]
test3_status = known_vogtmann$status
test4 = known_zeller[1:(ncol(known_zeller)-1)]
test4_status = known_zeller$status
test5 = known_thomas[1:(ncol(known_thomas)-1)]
test5_status = known_thomas$status

##training RF with 1000 decision trees and 30 reptitions
auc_cross_cohort <- rf_combat(training, training_status, test1, test1_status, auc_cross_cohort, 1000, 30, col=1)
auc_cross_cohort <- rf_combat(training, training_status, test2, test2_status, auc_cross_cohort, 1000, 30, col=2)
auc_cross_cohort <- rf_combat(training, training_status, test3, test3_status, auc_cross_cohort, 1000, 30, col=3)
auc_cross_cohort <- rf_combat(training, training_status, test4, test4_status, auc_cross_cohort, 1000, 30, col=4)
auc_cross_cohort <- rf_combat(training, training_status, test5, test5_status, auc_cross_cohort, 1000, 30, col=5)


