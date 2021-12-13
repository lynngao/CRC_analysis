library(readr)
library(metafor)
library(caret)
library(randomForest)
library(pROC)
library(operators)
library(glmnet)

##########################README#########################
##require species abundance table for calculating AUC####
#########################################################

#################################within-dataset functions####################################
within_rf <- function(dataset, num_trees, num_rep, df, col=1){
	#renormalization
	dataset[,1:(ncol(dataset)-1)] = dataset[,1:(ncol(dataset)-1)]/rowSums(dataset[,1:(ncol(dataset)-1)]) 
	
	for (i in 1:num_rep){
		sample_crc = sample.int(n=nrow(dataset[dataset$status == "CRC",]),size=floor(0.8*nrow(dataset[dataset$status == "CRC",])),replace=F)
  		sample_control = sample.int(n=nrow(dataset[dataset$status == "control",]),size=floor(0.8*nrow(dataset[dataset$status == "control",])),replace=F)
  		training_crc = dataset[dataset$status == "CRC",][sample_crc,]
  		training_control = dataset[dataset$status == "control",][sample_control,]
  		training = rbind(training_crc,training_control)
  		test = dataset[rownames(dataset)%!in%rownames(training),]
  
  		rf <- train(
    		status ~ ., data = training, method = "rf", ntree = num_trees, metric = "ROC",
    		trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  		)
  		pred <- predict(rf, test, type="prob")
  		df[i,col] <- auc(test$status, pred[,1])
  	}

  	return(df)
}

within_lasso <- function(dataset, num_rep, df, col=1){
	#renormalization
	dataset[,1:(ncol(dataset)-1)] = dataset[,1:(ncol(dataset)-1)]/rowSums(dataset[,1:(ncol(dataset)-1)]) 
	#log10 transformation
	dataset[,1:(ncol(dataset)-1)] = log10(dataset[,1:(ncol(dataset)-1)]+1e-100)

	parameters <- seq(0,0.5,0.001)

	for (i in 1:num_rep){
		sample_crc = sample.int(n=nrow(dataset[dataset$status == "CRC",]),size=floor(0.8*nrow(dataset[dataset$status == "CRC",])),replace=F)
  		sample_control = sample.int(n=nrow(dataset[dataset$status == "control",]),size=floor(0.8*nrow(dataset[dataset$status == "control",])),replace=F)
  		training_crc = dataset[dataset$status == "CRC",][sample_crc,]
  		training_control = dataset[dataset$status == "control",][sample_control,]
  		training = rbind(training_crc,training_control)
  		test = dataset[rownames(dataset)%!in%rownames(training),]
  
  		lassofit<- train(
    	status ~ ., data = training, method = "glmnet", metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = parameters),
    	trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  		)
  		pred <- predict(lassofit, test, type="prob")
  		df[i,col] <- auc(test$status, pred[,1])
  	}

  	return(df)
}

within_svm <- function(dataset, num_rep, df, col=1){
	#renormalization
	dataset[,1:(ncol(dataset)-1)] = dataset[,1:(ncol(dataset)-1)]/rowSums(dataset[,1:(ncol(dataset)-1)]) 
	#log10 transformation
	dataset[,1:(ncol(dataset)-1)] = log10(dataset[,1:(ncol(dataset)-1)]+1e-100)

	for (i in 1:num_rep){
		sample_crc = sample.int(n=nrow(dataset[dataset$status == "CRC",]),size=floor(0.8*nrow(dataset[dataset$status == "CRC",])),replace=F)
  		sample_control = sample.int(n=nrow(dataset[dataset$status == "control",]),size=floor(0.8*nrow(dataset[dataset$status == "control",])),replace=F)
  		training_crc = dataset[dataset$status == "CRC",][sample_crc,]
  		training_control = dataset[dataset$status == "control",][sample_control,]
  		training = rbind(training_crc,training_control)
  		test = dataset[rownames(dataset)%!in%rownames(training),]
  
  		rf <- train(
    		status ~ ., data = known_training, method = "svmRadial", metric = "ROC",
    		trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  		)
  		pred <- predict(rf, test, type="prob")
  		df[i,col] <- auc(test$status, pred[,1])
  	}

  	return(df)
}

#################################cross-dataset functions####################################
cross_rf <- function(training, training_status, test, test_status, num_trees, num_rep, df, col=1){
	training = training/rowSums(training) #renormalize
	for (j in 1:ncol(training)) {
  		if(colnames(training)[j] %!in% colnames(test)){
    	test = cbind(test,0)
    	colnames(test)[ncol(test)] = colnames(training)[j]
  		}
	}
	test = test/rowSums(test)

	training$status = training_status
	test$status = test_status

	for (i in 1:num_rep){
		rf <- train(
		  status ~ ., data = training, method = "rf", ntree = num_trees, metric = "ROC",
		  trControl = trainControl(method = "repeatedcv", number = 10,search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
		)
		pred <- predict(rf, test, type="prob")
		df[i,col] <- auc(test$status, pred[,1])
	}
	return(df)
}

cross_lasso <- function(training, training_status, test, test_status, num_rep, df, col=1){
	training = training/rowSums(training) #renormalize
	#log10 transformation
	training = log10(training+1e-100)

	for (j in 1:ncol(training)) {
		if(colnames(training)[j] %!in% colnames(test)){
	    	test = cbind(test,0)
	    	colnames(test)[ncol(test)] = colnames(training)[j]
	  	}
	}
	test = test[,colnames(test)%in%colnames(training)]
	test = test/rowSums(test)
	#log10 transformation
	test = log10(test+1e-100)

	training$status = training_status
	test$status = test_status

	parameters <- seq(0,0.5,0.001)

	for (i in 1:num_rep){
		lassofit<- train(
    	status ~ ., data = training, method = "glmnet", metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = parameters),
    	trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  		)
		pred <- predict(lassofit, test, type="prob")
		dataframe[i,col] <- auc(test$status, pred[,1])
	}
	return(df)
}

cross_svm <- function(training, training_status, test, test_status, num_rep, df, col=1){
	training = training/rowSums(training) #renormalize
	#log10 transformation
	training = log10(training+1e-100)
	for (j in 1:ncol(training)) {
		if(colnames(training)[j] %!in% colnames(test)){
	    	test = cbind(test,0)
	    	colnames(test)[ncol(test)] = colnames(training)[j]
	  	}
	}
	test = test[,colnames(test)%in%colnames(training)]
	test = test/rowSums(test)
	#log10 transformation
	test = log10(test+1e-100)

	training$status = training_status
	test$status = test_status

	for (i in 1:30){
		svm <- train(status ~ ., data = training, method = "svmRadial", metric = "ROC",
	    trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE))
		pred <- predict(svm, test, type="prob")
		df[i,col] <- auc(test$status, pred[,1])
	}
	return(df)
}


#################################LODO functions####################################
###########function for selecting only training species
feature_selection = function(training){
for (i in 1:length(features)) {
  if(features[i] %!in% colnames(training)){
    training = cbind(training,0)
    colnames(training)[ncol(training)] = features[i]
  }
}
training = training[,order(colnames(training))]
return(training)
}


LODO_rf <- function(ds1, ds2, ds3, ds4, ds5, test, num_trees, num_rep, df, col=1){
	features = unique(c(colnames(ds1),colnames(ds2),colnames(ds3),colnames(ds4),colnames(ds5)))
	known_1 = feature_selection(ds1)
	known_2 = feature_selection(ds2)
	known_3 = feature_selection(ds3)
	known_4 = feature_selection(ds4)
	known_5 = feature_selection(ds5)

	training = rbind(known_1,known_2,known_3,known_4,known_5)
	training[,1:(ncol(training)-1)] = training[,1:(ncol(training)-1)]/rowSums(training[,1:(ncol(training)-1)]) #renormalize
	for (i in 1:ncol(training)) {
	  	if(colnames(training)[i] %!in% colnames(test)){
	    	test = cbind(test,0)
	    	colnames(test)[ncol(test)] = colnames(training)[i]
	  	}
	}
	test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
	for (i in 1:num_rep){
		rf <- train(status ~ ., data = training, method = "rf", ntree = num_trees, metric = "ROC",
	  				trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE))
	pred <- predict(rf, test, type="prob")
	df[i,col] <- auc(test$status, pred[,1])
	}
}


LODO_lasso <- function(ds1, ds2, ds3, ds4, ds5, test, num_rep, df, lambda="lambda.min", col=1){
	features = unique(c(colnames(ds1),colnames(ds2),colnames(ds3),colnames(ds4),colnames(ds5)))
	known_1 = feature_selection(ds1)
	known_2 = feature_selection(ds2)
	known_3 = feature_selection(ds3)
	known_4 = feature_selection(ds4)
	known_5 = feature_selection(ds5)

	training = rbind(known_1,known_2,known_3,known_4,known_5)
	training[,1:(ncol(training)-1)] = training[,1:(ncol(training)-1)]/rowSums(training[,1:(ncol(training)-1)]) #renormalize
	#log10 transformation
	training[,1:(ncol(training)-1)] = log10(training[,1:(ncol(training)-1)]+1e-100)
	test = known_vogtmann
	for (i in 1:ncol(training)) {
	  	if(colnames(training)[i] %!in% colnames(test)){
	    	test = cbind(test,0)
	    	colnames(test)[ncol(test)] = colnames(training)[i]
	  	}
	}
	test = test[,colnames(test)%in%colnames(training)]
	test = test[,order(colnames(test))]
	test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
	#log10 transformation
	test[,-which(colnames(test)%in%c("status"))] = log10(test[,-which(colnames(test)%in%c("status"))]+1e-100)
	
	parameters <- seq(0,0.5,0.001)

	for (i in 1:num_rep){
		lassofit<- train(
    	status ~ ., data = training, method = "glmnet", metric = "ROC", tuneGrid = expand.grid(alpha = 1, lambda = parameters),
    	trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
  		)
		pred <- predict(lassofit, test, type="prob")
		dataframe[i,col] <- auc(test$status, pred[,1])
	}
}


LODO_svm <- function(ds1, ds2, ds3, ds4, ds5, test, num_rep, df, col=1){
	features = unique(c(colnames(ds1),colnames(ds2),colnames(ds3),colnames(ds4),colnames(ds5)))
	known_1 = feature_selection(ds1)
	known_2 = feature_selection(ds2)
	known_3 = feature_selection(ds3)
	known_4 = feature_selection(ds4)
	known_5 = feature_selection(ds5)

	training = rbind(known_1,known_2,known_3,known_4,known_5)
	training[,1:(ncol(training)-1)] = training[,1:(ncol(training)-1)]/rowSums(training[,1:(ncol(training)-1)]) #renormalize
	#log10 transformation
	training[,1:(ncol(training)-1)] = log10(training[,1:(ncol(training)-1)]+1e-100)
	test = known_vogtmann
	for (i in 1:ncol(training)) {
	  	if(colnames(training)[i] %!in% colnames(test)){
	    	test = cbind(test,0)
	    	colnames(test)[ncol(test)] = colnames(training)[i]
	  	}
	}
	test = test[,colnames(test)%in%colnames(training)]
	test = test[,order(colnames(test))]
	test[,-which(colnames(test)%in%c("status"))] = test[,-which(colnames(test)%in%c("status"))]/rowSums(test[,-which(colnames(test)%in%c("status"))])
	#log10 transformation
	test[,-which(colnames(test)%in%c("status"))] = log10(test[,-which(colnames(test)%in%c("status"))]+1e-100)
	for (i in 1:num_rep){
		svm <- train(
	    	status ~ ., data = training, method = "svmRadial", metric = "ROC",
	    	trControl = trainControl(method = "repeatedcv", number = 10, search = "grid", summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
		)
	pred <- predict(svm, test, type="prob")
	df[i,col] <- auc(test$status, pred[,1])
	}
}


