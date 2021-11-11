library(readr)
library(ggplot2)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)

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

##########training:yu
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

##test:hannigan
test_hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(test_hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(test_hannigan_AUC)){
	test_hannigan_AUC[i,1] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(test_hannigan_AUC)[-i]])
	test_hannigan_AUC[i,2] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(test_hannigan_AUC)[-i]])
}
test_hannigan_AUC[test_hannigan_AUC < 0.5] = 0.5
test_hannigan = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan) = c("avg_control","avg_CRC")
rownames(test_hannigan) = hannigan_known[,1]
test_hannigan$avg_control = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,2] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,2])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)
test_hannigan$avg_CRC = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,3])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)


##test:feng
test_feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng_AUC) = c("known_AUC","unknown_AUC")
rownames(test_feng_AUC) = feng_known[,1]
for (i in 1:nrow(test_feng_AUC)){
	test_feng_AUC[i,1] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(test_feng_AUC)[-i]])
	test_feng_AUC[i,2] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(test_feng_AUC)[-i]])
}
test_feng_AUC[test_feng_AUC < 0.5] = 0.5
test_feng = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng) = c("avg_control","avg_CRC")
rownames(test_feng) = feng_known[,1]
test_feng$avg_control = ((test_feng_AUC[,1]-0.5)*feng_known[,2] + (test_feng_AUC[,2]-0.5)*feng_unknown[,2])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)
test_feng$avg_CRC = ((test_feng_AUC[,1]-0.5)*feng_known[,3] + (test_feng_AUC[,2]-0.5)*feng_unknown[,3])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)


##test:vogtmann
test_vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(test_vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(test_vogtmann_AUC)){
	test_vogtmann_AUC[i,1] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(test_vogtmann_AUC)[-i]])
	test_vogtmann_AUC[i,2] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(test_vogtmann_AUC)[-i]])
}
test_vogtmann_AUC[test_vogtmann_AUC < 0.5] = 0.5
test_vogtmann = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann) = c("avg_control","avg_CRC")
rownames(test_vogtmann) = vogtmann_known[,1]
test_vogtmann$avg_control = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,2] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,2])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)
test_vogtmann$avg_CRC = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)


##test:zeller
test_zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(test_zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(test_zeller_AUC)){
	test_zeller_AUC[i,1] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(test_zeller_AUC)[-i]])
	test_zeller_AUC[i,2] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(test_zeller_AUC)[-i]])
}
test_zeller_AUC[test_zeller_AUC < 0.5] = 0.5
test_zeller = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller) = c("avg_control","avg_CRC")
rownames(test_zeller) = zeller_known[,1]
test_zeller$avg_control = ((test_zeller_AUC[,1]-0.5)*zeller_known[,2] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,2])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)
test_zeller$avg_CRC = ((test_zeller_AUC[,1]-0.5)*zeller_known[,3] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,3])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)


##test:thomas
test_thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(test_thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(test_thomas_AUC)){
	test_thomas_AUC[i,1] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(test_thomas_AUC)[-i]])
	test_thomas_AUC[i,2] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(test_thomas_AUC)[-i]])
}
test_thomas_AUC[test_thomas_AUC < 0.5] = 0.5
test_thomas = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas) = c("avg_control","avg_CRC")
rownames(test_thomas) = thomas_known[,1]
test_thomas$avg_control = ((test_thomas_AUC[,1]-0.5)*thomas_known[,2] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,2])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)
test_thomas$avg_CRC = ((test_thomas_AUC[,1]-0.5)*thomas_known[,3] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,3])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)

write.csv(test_hannigan, "/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_hannigan.csv")
write.csv(test_feng, "/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_feng.csv")
write.csv(test_vogtmann, "/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_vogtmann.csv")
write.csv(test_zeller, "/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_zeller.csv")
write.csv(test_thomas, "/Users/lynngao/Desktop/pred/cross_dataset/yu/LOSO_prob_test_thomas.csv")

auc(known_hannigan$status, test_hannigan$avg_CRC) #0.668
auc(known_feng$status, test_feng$avg_CRC) #0.8982
auc(known_vogtmann$status, test_vogtmann$avg_CRC) #0.7629
auc(known_zeller$status, test_zeller$avg_CRC) #0.8002
auc(known_thomas$status, test_thomas$avg_CRC) #0.7374



##########training:hannigan   
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


##test:yu
test_yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu_AUC) = c("known_AUC","unknown_AUC")
rownames(test_yu_AUC) = yu_known[,1]
for (i in 1:nrow(test_yu_AUC)){
	test_yu_AUC[i,1] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(test_yu_AUC)[-i]])
	test_yu_AUC[i,2] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(test_yu_AUC)[-i]])
}
test_yu_AUC[test_yu_AUC < 0.5] = 0.5
test_yu = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu) = c("avg_control","avg_CRC")
rownames(test_yu) = yu_known[,1]
test_yu$avg_control = ((test_yu_AUC[,1]-0.5)*yu_known[,2] + (test_yu_AUC[,2]-0.5)*yu_unknown[,2])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)
test_yu$avg_CRC = ((test_yu_AUC[,1]-0.5)*yu_known[,3] + (test_yu_AUC[,2]-0.5)*yu_unknown[,3])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)


##test:feng
test_feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng_AUC) = c("known_AUC","unknown_AUC")
rownames(test_feng_AUC) = feng_known[,1]
for (i in 1:nrow(test_feng_AUC)){
	test_feng_AUC[i,1] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(test_feng_AUC)[-i]])
	test_feng_AUC[i,2] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(test_feng_AUC)[-i]])
}
test_feng_AUC[test_feng_AUC < 0.5] = 0.5
test_feng = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng) = c("avg_control","avg_CRC")
rownames(test_feng) = feng_known[,1]
test_feng$avg_control = ((test_feng_AUC[,1]-0.5)*feng_known[,2] + (test_feng_AUC[,2]-0.5)*feng_unknown[,2])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)
test_feng$avg_CRC = ((test_feng_AUC[,1]-0.5)*feng_known[,3] + (test_feng_AUC[,2]-0.5)*feng_unknown[,3])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)


##test:vogtmann
test_vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(test_vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(test_vogtmann_AUC)){
	test_vogtmann_AUC[i,1] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(test_vogtmann_AUC)[-i]])
	test_vogtmann_AUC[i,2] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(test_vogtmann_AUC)[-i]])
}
test_vogtmann_AUC[test_vogtmann_AUC < 0.5] = 0.5
test_vogtmann = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann) = c("avg_control","avg_CRC")
rownames(test_vogtmann) = vogtmann_known[,1]
test_vogtmann$avg_control = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,2] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,2])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)
test_vogtmann$avg_CRC = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)


##test:zeller
test_zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(test_zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(test_zeller_AUC)){
	test_zeller_AUC[i,1] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(test_zeller_AUC)[-i]])
	test_zeller_AUC[i,2] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(test_zeller_AUC)[-i]])
}
test_zeller_AUC[test_zeller_AUC < 0.5] = 0.5
test_zeller = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller) = c("avg_control","avg_CRC")
rownames(test_zeller) = zeller_known[,1]
test_zeller$avg_control = ((test_zeller_AUC[,1]-0.5)*zeller_known[,2] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,2])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)
test_zeller$avg_CRC = ((test_zeller_AUC[,1]-0.5)*zeller_known[,3] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,3])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)


##test:thomas
test_thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(test_thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(test_thomas_AUC)){
	test_thomas_AUC[i,1] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(test_thomas_AUC)[-i]])
	test_thomas_AUC[i,2] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(test_thomas_AUC)[-i]])
}
test_thomas_AUC[test_thomas_AUC < 0.5] = 0.5
test_thomas = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas) = c("avg_control","avg_CRC")
rownames(test_thomas) = thomas_known[,1]
test_thomas$avg_control = ((test_thomas_AUC[,1]-0.5)*thomas_known[,2] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,2])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)
test_thomas$avg_CRC = ((test_thomas_AUC[,1]-0.5)*thomas_known[,3] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,3])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)

write.csv(test_yu, "/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_yu.csv")
write.csv(test_feng, "/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_feng.csv")
write.csv(test_vogtmann, "/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_vogtmann.csv")
write.csv(test_zeller, "/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_zeller.csv")
write.csv(test_thomas, "/Users/lynngao/Desktop/pred/cross_dataset/hannigan/LOSO_prob_test_thomas.csv")

auc(known_yu$status, test_yu$avg_CRC)
auc(known_feng$status, test_feng$avg_CRC)
auc(known_vogtmann$status, test_vogtmann$avg_CRC)
auc(known_zeller$status, test_zeller$avg_CRC)
auc(known_thomas$status, test_thomas$avg_CRC)


##########training:feng
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

##test:yu
test_yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu_AUC) = c("known_AUC","unknown_AUC")
rownames(test_yu_AUC) = yu_known[,1]
for (i in 1:nrow(test_yu_AUC)){
	test_yu_AUC[i,1] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(test_yu_AUC)[-i]])
	test_yu_AUC[i,2] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(test_yu_AUC)[-i]])
}
test_yu_AUC[test_yu_AUC < 0.5] = 0.5
test_yu = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu) = c("avg_control","avg_CRC")
rownames(test_yu) = yu_known[,1]
test_yu$avg_control = ((test_yu_AUC[,1]-0.5)*yu_known[,2] + (test_yu_AUC[,2]-0.5)*yu_unknown[,2])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)
test_yu$avg_CRC = ((test_yu_AUC[,1]-0.5)*yu_known[,3] + (test_yu_AUC[,2]-0.5)*yu_unknown[,3])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)


##test:hannigan
test_hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(test_hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(test_hannigan_AUC)){
	test_hannigan_AUC[i,1] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(test_hannigan_AUC)[-i]])
	test_hannigan_AUC[i,2] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(test_hannigan_AUC)[-i]])
}
test_hannigan_AUC[test_hannigan_AUC < 0.5] = 0.5
test_hannigan = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan) = c("avg_control","avg_CRC")
rownames(test_hannigan) = hannigan_known[,1]
test_hannigan$avg_control = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,2] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,2])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)
test_hannigan$avg_CRC = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,3])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)


##test:vogtmann
test_vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(test_vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(test_vogtmann_AUC)){
	test_vogtmann_AUC[i,1] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(test_vogtmann_AUC)[-i]])
	test_vogtmann_AUC[i,2] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(test_vogtmann_AUC)[-i]])
}
test_vogtmann_AUC[test_vogtmann_AUC < 0.5] = 0.5
test_vogtmann = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann) = c("avg_control","avg_CRC")
rownames(test_vogtmann) = vogtmann_known[,1]
test_vogtmann$avg_control = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,2] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,2])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)
test_vogtmann$avg_CRC = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)


##test:zeller
test_zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(test_zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(test_zeller_AUC)){
	test_zeller_AUC[i,1] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(test_zeller_AUC)[-i]])
	test_zeller_AUC[i,2] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(test_zeller_AUC)[-i]])
}
test_zeller_AUC[test_zeller_AUC < 0.5] = 0.5
test_zeller = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller) = c("avg_control","avg_CRC")
rownames(test_zeller) = zeller_known[,1]
test_zeller$avg_control = ((test_zeller_AUC[,1]-0.5)*zeller_known[,2] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,2])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)
test_zeller$avg_CRC = ((test_zeller_AUC[,1]-0.5)*zeller_known[,3] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,3])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)


##test:thomas
test_thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(test_thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(test_thomas_AUC)){
	test_thomas_AUC[i,1] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(test_thomas_AUC)[-i]])
	test_thomas_AUC[i,2] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(test_thomas_AUC)[-i]])
}
test_thomas_AUC[test_thomas_AUC < 0.5] = 0.5
test_thomas = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas) = c("avg_control","avg_CRC")
rownames(test_thomas) = thomas_known[,1]
test_thomas$avg_control = ((test_thomas_AUC[,1]-0.5)*thomas_known[,2] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,2])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)
test_thomas$avg_CRC = ((test_thomas_AUC[,1]-0.5)*thomas_known[,3] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,3])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)

write.csv(test_yu, "/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_yu.csv")
write.csv(test_hannigan, "/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_hannigan.csv")
write.csv(test_vogtmann, "/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_vogtmann.csv")
write.csv(test_zeller, "/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_zeller.csv")
write.csv(test_thomas, "/Users/lynngao/Desktop/pred/cross_dataset/feng/LOSO_prob_test_thomas.csv")

auc(known_yu$status, test_yu$avg_CRC)
auc(known_hannigan$status, test_hannigan$avg_CRC)
auc(known_vogtmann$status, test_vogtmann$avg_CRC)
auc(known_zeller$status, test_zeller$avg_CRC)
auc(known_thomas$status, test_thomas$avg_CRC)



##########training:vogtmann
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

##test:yu
test_yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu_AUC) = c("known_AUC","unknown_AUC")
rownames(test_yu_AUC) = yu_known[,1]
for (i in 1:nrow(test_yu_AUC)){
	test_yu_AUC[i,1] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(test_yu_AUC)[-i]])
	test_yu_AUC[i,2] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(test_yu_AUC)[-i]])
}
test_yu_AUC[test_yu_AUC < 0.5] = 0.5
test_yu = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu) = c("avg_control","avg_CRC")
rownames(test_yu) = yu_known[,1]
test_yu$avg_control = ((test_yu_AUC[,1]-0.5)*yu_known[,2] + (test_yu_AUC[,2]-0.5)*yu_unknown[,2])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)
test_yu$avg_CRC = ((test_yu_AUC[,1]-0.5)*yu_known[,3] + (test_yu_AUC[,2]-0.5)*yu_unknown[,3])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)


##test:hannigan
test_hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(test_hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(test_hannigan_AUC)){
	test_hannigan_AUC[i,1] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(test_hannigan_AUC)[-i]])
	test_hannigan_AUC[i,2] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(test_hannigan_AUC)[-i]])
}
test_hannigan_AUC[test_hannigan_AUC < 0.5] = 0.5
test_hannigan = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan) = c("avg_control","avg_CRC")
rownames(test_hannigan) = hannigan_known[,1]
test_hannigan$avg_control = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,2] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,2])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)
test_hannigan$avg_CRC = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,3])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)


##test:feng
test_feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng_AUC) = c("known_AUC","unknown_AUC")
rownames(test_feng_AUC) = feng_known[,1]
for (i in 1:nrow(test_feng_AUC)){
	test_feng_AUC[i,1] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(test_feng_AUC)[-i]])
	test_feng_AUC[i,2] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(test_feng_AUC)[-i]])
}
test_feng_AUC[test_feng_AUC < 0.5] = 0.5
test_feng = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng) = c("avg_control","avg_CRC")
rownames(test_feng) = feng_known[,1]
test_feng$avg_control = ((test_feng_AUC[,1]-0.5)*feng_known[,2] + (test_feng_AUC[,2]-0.5)*feng_unknown[,2])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)
test_feng$avg_CRC = ((test_feng_AUC[,1]-0.5)*feng_known[,3] + (test_feng_AUC[,2]-0.5)*feng_unknown[,3])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)


##test:zeller
test_zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(test_zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(test_zeller_AUC)){
	test_zeller_AUC[i,1] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(test_zeller_AUC)[-i]])
	test_zeller_AUC[i,2] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(test_zeller_AUC)[-i]])
}
test_zeller_AUC[test_zeller_AUC < 0.5] = 0.5
test_zeller = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller) = c("avg_control","avg_CRC")
rownames(test_zeller) = zeller_known[,1]
test_zeller$avg_control = ((test_zeller_AUC[,1]-0.5)*zeller_known[,2] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,2])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)
test_zeller$avg_CRC = ((test_zeller_AUC[,1]-0.5)*zeller_known[,3] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,3])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)


##test:thomas
test_thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(test_thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(test_thomas_AUC)){
	test_thomas_AUC[i,1] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(test_thomas_AUC)[-i]])
	test_thomas_AUC[i,2] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(test_thomas_AUC)[-i]])
}
test_thomas_AUC[test_thomas_AUC < 0.5] = 0.5
test_thomas = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas) = c("avg_control","avg_CRC")
rownames(test_thomas) = thomas_known[,1]
test_thomas$avg_control = ((test_thomas_AUC[,1]-0.5)*thomas_known[,2] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,2])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)
test_thomas$avg_CRC = ((test_thomas_AUC[,1]-0.5)*thomas_known[,3] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,3])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)

write.csv(test_yu, "/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_yu.csv")
write.csv(test_hannigan, "/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_hannigan.csv")
write.csv(test_feng, "/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_feng.csv")
write.csv(test_zeller, "/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_zeller.csv")
write.csv(test_thomas, "/Users/lynngao/Desktop/pred/cross_dataset/vogtmann/LOSO_prob_test_thomas.csv")


auc(known_yu$status, test_yu$avg_CRC)
auc(known_hannigan$status, test_hannigan$avg_CRC)
auc(known_feng$status, test_feng$avg_CRC)
auc(known_zeller$status, test_zeller$avg_CRC)
auc(known_thomas$status, test_thomas$avg_CRC)


##########training:zeller
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

##test:yu
test_yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu_AUC) = c("known_AUC","unknown_AUC")
rownames(test_yu_AUC) = yu_known[,1]
for (i in 1:nrow(test_yu_AUC)){
	test_yu_AUC[i,1] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(test_yu_AUC)[-i]])
	test_yu_AUC[i,2] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(test_yu_AUC)[-i]])
}
test_yu_AUC[test_yu_AUC < 0.5] = 0.5
test_yu = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu) = c("avg_control","avg_CRC")
rownames(test_yu) = yu_known[,1]
test_yu$avg_control = ((test_yu_AUC[,1]-0.5)*yu_known[,2] + (test_yu_AUC[,2]-0.5)*yu_unknown[,2])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)
test_yu$avg_CRC = ((test_yu_AUC[,1]-0.5)*yu_known[,3] + (test_yu_AUC[,2]-0.5)*yu_unknown[,3])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)


##test:hannigan
test_hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(test_hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(test_hannigan_AUC)){
	test_hannigan_AUC[i,1] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(test_hannigan_AUC)[-i]])
	test_hannigan_AUC[i,2] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(test_hannigan_AUC)[-i]])
}
test_hannigan_AUC[test_hannigan_AUC < 0.5] = 0.5
test_hannigan = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan) = c("avg_control","avg_CRC")
rownames(test_hannigan) = hannigan_known[,1]
test_hannigan$avg_control = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,2] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,2])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)
test_hannigan$avg_CRC = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,3])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)


##test:feng
test_feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng_AUC) = c("known_AUC","unknown_AUC")
rownames(test_feng_AUC) = feng_known[,1]
for (i in 1:nrow(test_feng_AUC)){
	test_feng_AUC[i,1] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(test_feng_AUC)[-i]])
	test_feng_AUC[i,2] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(test_feng_AUC)[-i]])
}
test_feng_AUC[test_feng_AUC < 0.5] = 0.5
test_feng = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng) = c("avg_control","avg_CRC")
rownames(test_feng) = feng_known[,1]
test_feng$avg_control = ((test_feng_AUC[,1]-0.5)*feng_known[,2] + (test_feng_AUC[,2]-0.5)*feng_unknown[,2])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)
test_feng$avg_CRC = ((test_feng_AUC[,1]-0.5)*feng_known[,3] + (test_feng_AUC[,2]-0.5)*feng_unknown[,3])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)


##test:vogtmann
test_vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(test_vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(test_vogtmann_AUC)){
	test_vogtmann_AUC[i,1] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(test_vogtmann_AUC)[-i]])
	test_vogtmann_AUC[i,2] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(test_vogtmann_AUC)[-i]])
}
test_vogtmann_AUC[test_vogtmann_AUC < 0.5] = 0.5
test_vogtmann = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann) = c("avg_control","avg_CRC")
rownames(test_vogtmann) = vogtmann_known[,1]
test_vogtmann$avg_control = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,2] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,2])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)
test_vogtmann$avg_CRC = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)


##test:thomas
test_thomas_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas_AUC) = c("known_AUC","unknown_AUC")
rownames(test_thomas_AUC) = thomas_known[,1]
for (i in 1:nrow(test_thomas_AUC)){
	test_thomas_AUC[i,1] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_known$CRC[thomas_known[,1]%in%rownames(test_thomas_AUC)[-i]])
	test_thomas_AUC[i,2] = auc(known_thomas$status[rownames(known_thomas)%in%rownames(test_thomas_AUC)[-i]], thomas_unknown$CRC[thomas_unknown[,1]%in%rownames(test_thomas_AUC)[-i]])
}
test_thomas_AUC[test_thomas_AUC < 0.5] = 0.5
test_thomas = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(thomas_known)))
colnames(test_thomas) = c("avg_control","avg_CRC")
rownames(test_thomas) = thomas_known[,1]
test_thomas$avg_control = ((test_thomas_AUC[,1]-0.5)*thomas_known[,2] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,2])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)
test_thomas$avg_CRC = ((test_thomas_AUC[,1]-0.5)*thomas_known[,3] + (test_thomas_AUC[,2]-0.5)*thomas_unknown[,3])/(test_thomas_AUC[,1]+test_thomas_AUC[,2]-1)

write.csv(test_yu, "/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_yu.csv")
write.csv(test_hannigan, "/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_hannigan.csv")
write.csv(test_feng, "/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_feng.csv")
write.csv(test_vogtmann, "/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_vogtmann.csv")
write.csv(test_thomas, "/Users/lynngao/Desktop/pred/cross_dataset/zeller/LOSO_prob_test_thomas.csv")


auc(known_yu$status, test_yu$avg_CRC)
auc(known_hannigan$status, test_hannigan$avg_CRC)
auc(known_feng$status, test_feng$avg_CRC)
auc(known_vogtmann$status, test_vogtmann$avg_CRC)
auc(known_thomas$status, test_thomas$avg_CRC)



##########training:thomas
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

##test:yu
test_yu_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu_AUC) = c("known_AUC","unknown_AUC")
rownames(test_yu_AUC) = yu_known[,1]
for (i in 1:nrow(test_yu_AUC)){
	test_yu_AUC[i,1] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_known$CRC[yu_known[,1]%in%rownames(test_yu_AUC)[-i]])
	test_yu_AUC[i,2] = auc(known_yu$status[rownames(known_yu)%in%rownames(test_yu_AUC)[-i]], yu_unknown$CRC[yu_unknown[,1]%in%rownames(test_yu_AUC)[-i]])
}
test_yu_AUC[test_yu_AUC < 0.5] = 0.5
test_yu = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(yu_known)))
colnames(test_yu) = c("avg_control","avg_CRC")
rownames(test_yu) = yu_known[,1]
test_yu$avg_control = ((test_yu_AUC[,1]-0.5)*yu_known[,2] + (test_yu_AUC[,2]-0.5)*yu_unknown[,2])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)
test_yu$avg_CRC = ((test_yu_AUC[,1]-0.5)*yu_known[,3] + (test_yu_AUC[,2]-0.5)*yu_unknown[,3])/(test_yu_AUC[,1]+test_yu_AUC[,2]-1)


##test:hannigan
test_hannigan_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan_AUC) = c("known_AUC","unknown_AUC")
rownames(test_hannigan_AUC) = hannigan_known[,1]
for (i in 1:nrow(test_hannigan_AUC)){
	test_hannigan_AUC[i,1] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_known$CRC[hannigan_known[,1]%in%rownames(test_hannigan_AUC)[-i]])
	test_hannigan_AUC[i,2] = auc(known_hannigan$status[rownames(known_hannigan)%in%rownames(test_hannigan_AUC)[-i]], hannigan_unknown$CRC[hannigan_unknown[,1]%in%rownames(test_hannigan_AUC)[-i]])
}
test_hannigan_AUC[test_hannigan_AUC < 0.5] = 0.5
test_hannigan = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(hannigan_known)))
colnames(test_hannigan) = c("avg_control","avg_CRC")
rownames(test_hannigan) = hannigan_known[,1]
test_hannigan$avg_control = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,2] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,2])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)
test_hannigan$avg_CRC = ((test_hannigan_AUC[,1]-0.5)*hannigan_known[,3] + (test_hannigan_AUC[,2]-0.5)*hannigan_unknown[,3])/(test_hannigan_AUC[,1]+test_hannigan_AUC[,2]-1)


##test:feng
test_feng_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng_AUC) = c("known_AUC","unknown_AUC")
rownames(test_feng_AUC) = feng_known[,1]
for (i in 1:nrow(test_feng_AUC)){
	test_feng_AUC[i,1] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_known$CRC[feng_known[,1]%in%rownames(test_feng_AUC)[-i]])
	test_feng_AUC[i,2] = auc(known_feng$status[rownames(known_feng)%in%rownames(test_feng_AUC)[-i]], feng_unknown$CRC[feng_unknown[,1]%in%rownames(test_feng_AUC)[-i]])
}
test_feng_AUC[test_feng_AUC < 0.5] = 0.5
test_feng = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(feng_known)))
colnames(test_feng) = c("avg_control","avg_CRC")
rownames(test_feng) = feng_known[,1]
test_feng$avg_control = ((test_feng_AUC[,1]-0.5)*feng_known[,2] + (test_feng_AUC[,2]-0.5)*feng_unknown[,2])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)
test_feng$avg_CRC = ((test_feng_AUC[,1]-0.5)*feng_known[,3] + (test_feng_AUC[,2]-0.5)*feng_unknown[,3])/(test_feng_AUC[,1]+test_feng_AUC[,2]-1)


##test:vogtmann
test_vogtmann_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann_AUC) = c("known_AUC","unknown_AUC")
rownames(test_vogtmann_AUC) = vogtmann_known[,1]
for (i in 1:nrow(test_vogtmann_AUC)){
	test_vogtmann_AUC[i,1] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_known$CRC[vogtmann_known[,1]%in%rownames(test_vogtmann_AUC)[-i]])
	test_vogtmann_AUC[i,2] = auc(known_vogtmann$status[rownames(known_vogtmann)%in%rownames(test_vogtmann_AUC)[-i]], vogtmann_unknown$CRC[vogtmann_unknown[,1]%in%rownames(test_vogtmann_AUC)[-i]])
}
test_vogtmann_AUC[test_vogtmann_AUC < 0.5] = 0.5
test_vogtmann = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(vogtmann_known)))
colnames(test_vogtmann) = c("avg_control","avg_CRC")
rownames(test_vogtmann) = vogtmann_known[,1]
test_vogtmann$avg_control = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,2] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,2])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)
test_vogtmann$avg_CRC = ((test_vogtmann_AUC[,1]-0.5)*vogtmann_known[,3] + (test_vogtmann_AUC[,2]-0.5)*vogtmann_unknown[,3])/(test_vogtmann_AUC[,1]+test_vogtmann_AUC[,2]-1)


##test:zeller
test_zeller_AUC = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller_AUC) = c("known_AUC","unknown_AUC")
rownames(test_zeller_AUC) = zeller_known[,1]
for (i in 1:nrow(test_zeller_AUC)){
	test_zeller_AUC[i,1] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_known$CRC[zeller_known[,1]%in%rownames(test_zeller_AUC)[-i]])
	test_zeller_AUC[i,2] = auc(known_zeller$status[rownames(known_zeller)%in%rownames(test_zeller_AUC)[-i]], zeller_unknown$CRC[zeller_unknown[,1]%in%rownames(test_zeller_AUC)[-i]])
}
test_zeller_AUC[test_zeller_AUC < 0.5] = 0.5
test_zeller = as.data.frame(matrix(NA, ncol = 2, nrow = nrow(zeller_known)))
colnames(test_zeller) = c("avg_control","avg_CRC")
rownames(test_zeller) = zeller_known[,1]
test_zeller$avg_control = ((test_zeller_AUC[,1]-0.5)*zeller_known[,2] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,2])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)
test_zeller$avg_CRC = ((test_zeller_AUC[,1]-0.5)*zeller_known[,3] + (test_zeller_AUC[,2]-0.5)*zeller_unknown[,3])/(test_zeller_AUC[,1]+test_zeller_AUC[,2]-1)

write.csv(test_yu, "/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_yu.csv")
write.csv(test_hannigan, "/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_hannigan.csv")
write.csv(test_feng, "/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_feng.csv")
write.csv(test_vogtmann, "/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_vogtmann.csv")
write.csv(test_zeller, "/Users/lynngao/Desktop/pred/cross_dataset/thomas/LOSO_prob_test_zeller.csv")

auc(known_yu$status, test_yu$avg_CRC)
auc(known_hannigan$status, test_hannigan$avg_CRC)
auc(known_feng$status, test_feng$avg_CRC)
auc(known_vogtmann$status, test_vogtmann$avg_CRC)
auc(known_zeller$status, test_zeller$avg_CRC)







