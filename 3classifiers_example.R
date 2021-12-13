##load helper file
source("3classifiers.R")

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


###########################example of within-dataset classifiers#########################################
##creating dataframe to save results
auc_df = as.data.frame(matrix(NA,ncol = 3, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = c("RF","LASSO","SVM")

##random forests classifier with 1000 decision trees, 30 reptitions, training on Yu dataset
auc_df = within_rf(known_yu, 1000, 30, auc_df)

##LASSO classifier with 30 reptitions, training on Yu dataset
auc_df = within_lasso(known_yu, 30, auc_df, col=2)

##SVM classifier with 30 reptitions, training on Yu dataset
auc_df = within_svm(known_yu, 30, auc_df, col=3)
#########################################################################################################

###########################example of cross-dataset classifiers##########################################
##creating dataframe to save results
auc_df = as.data.frame(matrix(NA,ncol = 3, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = c("RF","LASSO","SVM")

##training on Yu dataset, testing on Hannigan dataset
training = known_yu[1:(ncol(known_yu)-1)]
training_status = known_yu$status
test = known_hannigan[1:(ncol(known_hannigan)-1)]
test_status = known_hannigan$status

##random forests classifier with 1000 decision trees, 30 reptitions
auc_df = cross_rf(training, training_status, test, test_status, 1000, 30, auc_df)

##LASSO classifier with 30 reptitions
auc_df = cross_lasso(training, training_status, test, test_status, 30, auc_df, col=2)

##SVM classifier with 30 reptitions
auc_df = cross_svm(training, training_status, test, test_status, 30, auc_df, col=3)
#########################################################################################################

####################################example of LODO classifiers##########################################
##creating dataframe to save results
auc_df = as.data.frame(matrix(NA,ncol = 3, nrow = 32))
rownames(auc_known) <- c(1:30, "mean", "sd")
colnames(auc_known) = c("RF","LASSO","SVM")

##testing on Yu dataset, the other five datasets are combined as training

##random forests classifier with 1000 decision trees, 10 reptitions
auc_df = LODO_rf(known_hannigan, known_feng, known_vogtmann, known_zeller, known_thomas, test=known_yu, 1000, 10, auc_df)

##LASSO classifier with 10 reptitions
auc_df = LODO_lasso(known_hannigan, known_feng, known_vogtmann, known_zeller, known_thomas, test=known_yu, 10, auc_df, col=2)

##SVM classifier with 10 reptitions
auc_df = LODO_svm(known_hannigan, known_feng, known_vogtmann, known_zeller, known_thomas, test=known_yu, 10, auc_df, col=3)
#########################################################################################################









