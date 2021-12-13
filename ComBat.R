library(readr)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(sva)

##########################README#########################
##require species abundance table for calculating AUC####
#########################################################

##############ComBat normalization is only available for cross-dataset analysi####################
rf_combat = function(training, training_status, test, test_status, num_trees, num_rep, dataframe, col=1){
	training = training/rowSums(training) #renormalize
	training = training[,order(colnames(training))]

	for (j in 1:ncol(training)) {
		if(colnames(training)[j] %!in% colnames(test)){
			test = cbind(test,0)
			colnames(test)[ncol(test)] = colnames(training)[j]
		}
	}
	test = test/rowSums(test) #renormalize
	test = test[,colnames(test)%in%colnames(training)]
	test = test[,order(colnames(test))]
	training$batch = 0 #add 0 for training batch
	test$batch = 1 #add 1 for testing batch

	data = rbind(training, test)
	batch = data$batch
	data = t(data[,1:(ncol(data)-1)])
	data_norm = as.data.frame(t(ComBat(data, batch)))
	training_norm = data_norm[1:nrow(training),]
	training_norm$status = training_status
	test_norm = data_norm[-(1:nrow(training)),]
	test_norm$status = test_status

	for (i in 1:num_rep){
		rf <- train(
		  status ~ ., data = training_norm, method = "rf", ntree = num_trees, metric = "ROC",
		  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
		)
		pred <- predict(rf, test_norm, type="prob")
		dataframe[i,col] <- auc(test_norm$status, pred[,1])
	}
	return(dataframe)
}
