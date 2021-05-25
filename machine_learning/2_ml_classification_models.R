##run different classification machine learning methods for Sandy 
set.seed(1234)

library(caret)
library(randomForest)
library(glmnet)
library(e1071)
library(gbm)

setwd("results")
inputDir = "results/"

##function to combine the results in comb and res into one data frame, returning the combined dataframe
combineResults <- function(comb, res) {
  if(!all(is.na(comb))) {
    comb = rbind(comb, res)
  } else {
    comb = res
  }
  return(comb)
}

##function to apply the given machine learning function fun on each quantification, and returns the combined results from the 4 models
##fun = machine learning function to apply
##funName = function name, to be used to identify the results
##any = predictors for any quantification
##pct = predictors for percent quantification
##hlf = predictors for half quantification
##mixed = predictors for mix quantification
callFunction <- function(fun, funName, any, pct, hlf, mixed) {
  print(funName)
  comb = NA
  if(!all(is.na(any))) {
    res = fun(any, paste0(funName, "_train_any"), train)
    comb = combineResults(comb, res)
  }
  if(!all(is.na(pct))) {
    res = fun(pct, paste0(funName, "_train_pct"), train)
    comb = combineResults(comb, res)
  }
  if(!all(is.na(hlf))) {
    res = fun(hlf, paste0(funName, "_train_half"), train)
    comb = combineResults(comb, res)
  }
  if(!all(is.na(mixed))) {
    res = fun(mixed, paste0(funName, "_train_mix"), train)
    comb = combineResults(comb, res)
  }
  return(comb)
}

##for the given model, return the statistics of the model performance on the given data
##model = model used to get the performance statistics
##preds = the set of predictors used to generate the model
##name = the name of the algorithm/set/group, used for identification in returned table
##data = the data to be used for the predictions and stats
##dataName = the name of the data set, used for identification in returned table
##method = the method for prediction (different predict functions are needed for different models)
getStats <- function(model, preds, name, data, dataName, method = "class") { 
  ## get prediction for the data
  if(method == "class") {
    pred = factor(predict(model, data[,preds], type="class"), ordered = T, levels = c("outbreak", "cluster"))
  } else if(method == "response") {
    pred = factor(ifelse(predict(model, data[,preds], type="response") < 0.5, "cluster", "outbreak"), 
                  ordered = T, levels = c("outbreak", "cluster"))
  } else if(method == "partyTree") {
    pred = factor(predict(model, data[,preds], type="response"), levels = c("outbreak", "cluster"), ordered=T)
  } else if(method == "rfParty") {
    pred = factor(predict(model, newdata = data[,preds], type="response", OOB=T), 
                  levels = c("outbreak", "cluster"), ordered = T)
  } else if(method == "glmnet") {
    x <- model.matrix(outcome ~., data[,c(preds, "outcome")])[,-1]
    pred = factor(ifelse(predict(model, x, type="response") < 0.5, "cluster", "outbreak"), 
                  levels = c("outbreak", "cluster"), ordered = T)
  } else if(method == "elasticNet") {
    x <- model.matrix(outcome ~., data[,c(preds, "outcome")])[,-1]
    pred = factor(predict(model, x, type="raw"), ordered = T, levels = c("outbreak", "cluster"))
  } else if(method == "gbm") {
    pred = factor(ifelse(predict(model, data[,preds], n.trees = n.gbm.trees, type="response") < 0.5, "cluster", "outbreak"),
                  levels = c("outbreak", "cluster"), ordered = T)
  } else if(method == "xgboost") {
    tmp = data[,preds]
    for(p in preds){#xgboost doesn't do factors
      tmp[,p] = as.numeric(as.character(tmp[,p]))
    }
    pred = factor(ifelse(predict(model, as.matrix(tmp)) < 0.5, "cluster", "outbreak"),
                  levels = c("outbreak", "cluster"), ordered = T)
  } else if(method == "ada") {
    pred = factor(predict(model, data[,preds]),
                  levels = c("outbreak", "cluster"), ordered = T)
  } else if(method == "empty") { #return empty result
    res = data.frame(alg_set_group = name,
                     accuracy = NA,
                     balancedAccuracy = NA,
                     sensitivity = NA,
                     specificity = NA,
                     PPV = NA,
                     NPV = NA)
    names(res)[-1] = paste0(names(res)[-1], ".", dataName)
    return(res)
  } else if(method == "NN") {
    pred = factor(ifelse(predict(nn, test, type="response")[,1] < 0.5, "cluster", "outbreak"), 
                  levels = c("outbreak", "cluster"))
  } else {
    stop("Bad stats method", method)
  }
  
  ##get performance statistics
  cm = confusionMatrix(data=pred, reference = data$outcome, positive = "outbreak")
  sens = cm$byClass[1]
  spec = cm$byClass[2]
  
  res = data.frame(alg_set_group = name,
                   accuracy = cm$overall[1],
                   balancedAccuracy = cm$byClass[11],
                   sensitivity = sens,
                   specificity = spec,
                   PPV = cm$byClass[3],
                   NPV = cm$byClass[4],
                   youden = sens + spec - 1)
  names(res)[-1] = paste0(names(res)[-1], ".", dataName)
  return(res)
}

##wrapper function to get results for test set
getStatsWrapper <- function(model, preds, name, method="class") {
  testRes = getStats(model, preds, name, test, "test", method)
  return(testRes)
}

###functions to generate each ML model
##each requires preds (the set of predictors), name (name to write results to), and train (the training set to use)

##random forest
rfTree <- function(preds, name, train) {
  data = train[,c(preds, "outcome")]
  data$outcome = factor(data$outcome)
  rf = randomForest(outcome ~ ., data = data)
  return(getStatsWrapper(rf, preds, name))
}

##lasso penalized regression model
lassoRegression <- function(preds, name, train) {
  x = model.matrix(outcome ~., train[,c(preds, "outcome")])[,-1] #-1 takes out intercept
  y = ifelse(train$outcome=="outbreak", 1, 0)
  ##Find the best lambda using cross-validation
  cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
  ##Fit the final model on the training data
  model <- glmnet(x, y, alpha = 1, family = "binomial",
                  lambda = cv.lasso$lambda.min)
  return(getStatsWrapper(model, preds, name, method = "glmnet"))
}

##ridge regression model
ridgeRegression <- function(preds, name, train) {
  x = model.matrix(outcome ~., train[,c(preds, "outcome")])[,-1] #-1 takes out intercept
  y = ifelse(train$outcome=="outbreak", 1, 0)
  ##Find the best lambda using cross-validation
  cv.lasso <- cv.glmnet(x, y, alpha = 0, family = "binomial")
  ##Fit the final model on the training data
  model <- glmnet(x, y, alpha = 0, family = "binomial",
                  lambda = cv.lasso$lambda.min)
  return(getStatsWrapper(model, preds, name, method = "glmnet"))
}

##elastic net model
elasticNet <- function(preds, name, train) {
  x = model.matrix(outcome ~., train[,c(preds, "outcome")])[,-1] #-1 takes out intercept
  y = factor(train$outcome)
  model <- train(x, y, 
                 method = "glmnet",
                 trControl = trainControl("cv", number = 10),
                 tuneLength = 10)
  return(getStatsWrapper(model, preds, name, method = "elasticNet"))
}

##SVM
##SVM with linear kernel
linearSupportVectorMachine <- function(preds, name, train) {
  data = train[,c(preds, "outcome")]
  data$outcome = as.factor(data$outcome)
  model = svm(outcome ~ ., data=data, kernel = "linear")
  return(getStatsWrapper(model, preds, name))
}
##SVM with polynomial kernel
polySupportVectorMachine <- function(preds, name, train) {
  data = train[,c(preds, "outcome")]
  data$outcome = as.factor(data$outcome)
  model = svm(outcome ~ ., data=data, kernel="polynomial")
  return(getStatsWrapper(model, preds, name))
}
##SVM with radial kernel
radialSupportVectorMachine <- function(preds, name, train) {
  data = train[,c(preds, "outcome")]
  data$outcome = as.factor(data$outcome)
  model = svm(outcome ~ ., data=data, kernel="radial")
  return(getStatsWrapper(model, preds, name))
}
##SVM with sigmoid kernel
sigmoidSupportVectorMachine <- function(preds, name, train) {
  data = train[,c(preds, "outcome")]
  data$outcome = as.factor(data$outcome)
  model = svm(outcome ~ ., data=data, kernel="sigmoid")
  return(getStatsWrapper(model, preds, name))
}

##gradient boosting machine
n.gbm.trees = 100 #default
gradientBoosting <- function(preds, name, train) {
  data = train[,c(preds, "outcome")]
  data$outcome = ifelse(data$outcome=="cluster", 1, 0)
  model = gbm(outcome ~ ., data=data, distribution = "bernoulli", n.trees = n.gbm.trees)
  return(getStatsWrapper(model, preds, name, method = "gbm"))
}


####generate all models
flag = "first" #use the first flag only
cohorts = c("1Q", "2Q", "4Q", "9Q") #time frames to iterate over
yr = paste0(1, "yrfu") #follow up period number of years for outcome
fs = "noFS" #no feature selection
repeats = 1:5 #number of repeats in cross validation
##4 timeframes
for(co in cohorts) {
  ##5 CV repeats
  for(r in repeats) {
    ex = paste0("trainingSet", r, "C_", flag, "_", yr, "_", co)  # name of this iteration
    print(paste(ex, fs))
    
    ##get train and test set
    train = read.csv(paste0(inputDir, ex, ".csv"), stringsAsFactors = F, header = T)
    test = read.csv(paste0(inputDir, "testSet", r, "C_", flag, "_", yr, "_", co, ".csv"), stringsAsFactors = F, header = T)
    train$outcome = factor(train$outcome, ordered = T, levels = c("outbreak", "cluster"))
    test$outcome = factor(test$outcome, ordered = T, levels = c("outbreak", "cluster"))

    ###set up predictors
    vfname = paste0(inputDir, "variablesList_comb_", flag, "_", co)
    ##mixed variable set
    fname = paste0(vfname, "_all.txt")
    if(file.exists(fname)) {
      mixed = read.table(fname, header = F, stringsAsFactors = F)
      mixed = as.character(mixed$V1)
    } else {
      mixed = NA
      print(paste("no mixed for", ex, fs))
    }
    
    ##any variable set
    fname = paste0(vfname, "_any.txt")
    if(file.exists(fname)) {
      any = read.table(fname, header = F, stringsAsFactors = F)
      any = as.character(any$V1)
      any = gsub("1", "", any) #for factor, 1 gets added on end
    } else {
      any = NA
      print(paste("no any for", ex, fs))
    }
    
    ##percent variable set
    fname = paste0(vfname, "_pct.txt")
    if(file.exists(fname)) {
      pct = read.table(fname, header = F, stringsAsFactors = F)
      pct = as.character(pct$V1)
    } else {
      pct = NA
      print(paste("no pct for", ex, fs))
    }
    
    ##half variable set
    fname = paste0(vfname, "_half.txt")
    if(file.exists(fname)) {
      hlf = read.table(fname, header = F, stringsAsFactors = F)
      hlf = as.character(hlf$V1)
      hlf = gsub("1", "", hlf) #for factor, 1 gets added on end
    } else {
      hlf = NA
      print(paste("no half for", ex, fs))
    }
    
    if(length(mixed) < 2 & length(any) < 2 & length(pct) < 2 & length(hlf) < 2) {
      print(paste(ex, "has no variable combos that passed feature selection"))
      next()
    }
    
    ##fix names
    names(train) = gsub("prev_", "prevtb_", names(train)) 
    names(test) = gsub("prev_", "prevtb_", names(test)) 
    
    ##list of continuous and factored variables
    facVars = c(names(train)[c(grep("half", names(train)), grep("any", names(train)))], "initflag", "incident")
    contVars = c(names(train)[grep("pct", names(train))], "casecnt", "timediff", "maxaccumxs", "nbhflag", "qtrsbtwflags", "nbhbase", "nbhcount", "prebasecases")
    facVars = facVars[facVars %in% names(train)]
    contVars = contVars[contVars %in% names(train)]
    
    ##fix classes (factor or make numeric)
    for(c in facVars) {
      train[,c] = factor(train[,c], levels=c(1,2))
      test[,c] = factor(test[,c], levels=c(1,2))
    }
    for(c in contVars) { #sometimes continuous variables are read as integers
      train[,c] = as.numeric(train[,c])
      test[,c] = as.numeric(test[,c])
    }
    
    ex = paste0("C_", yr, "_", co, "_", flag, "_", fs, "_Set", r)
    
    ###run each model
    summary = callFunction(rfTree, paste(ex, "randomForest", sep="_"), any, pct, hlf, mixed)
    summary = rbind(summary, callFunction(lassoRegression, paste(ex, "lassoRegression", sep="_"), any, pct, hlf, mixed))
    summary = rbind(summary, callFunction(ridgeRegression, paste(ex, "ridgeRegression", sep="_"), any, pct, hlf, mixed))
    summary = rbind(summary, callFunction(elasticNet, paste(ex, "elasticNet", sep="_"), any, pct, hlf, mixed))
    summary = rbind(summary, callFunction(linearSupportVectorMachine, paste(ex, "svmlinear", sep="_"), any, pct, hlf, mixed))
    summary = rbind(summary, callFunction(polySupportVectorMachine, paste(ex, "svmpoly", sep="_"), any, pct, hlf, mixed))
    summary = rbind(summary, callFunction(radialSupportVectorMachine, paste(ex, "svmradial", sep="_"), any, pct, hlf, mixed))
    summary = rbind(summary, callFunction(sigmoidSupportVectorMachine, paste(ex, "svmsigmoid", sep="_"), any, pct, hlf, mixed))
    summary = rbind(summary, callFunction(gradientBoosting, paste(ex, "gbm", sep="_"), any, pct, hlf, mixed))
    
    write.csv(summary, paste0("concatenatedModelStats_", ex, ".csv"), row.names = F, quote = F, na = "")
  }
}