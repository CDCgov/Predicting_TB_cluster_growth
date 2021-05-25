###generate final model using best model from cross validation (random forest on 2Q with half quantification) and test on validation set
set.seed(1234)

library(randomForest)
library(caret)

setwd("results")
inputDir = "../data/"

##set up training table
tab = read.csv(paste0(inputDir, "Fulldata2q_072920.csv"), header = T, stringsAsFactors = F)
names(tab) = gsub("prev_", "prevtb_", names(tab)) #these are the same variables but have different starts
tab$initflag = ifelse(tab$initflag=="Y", 2, 1)
tab$incident = ifelse(tab$clustyp=="Incident", 2, 1)
tab$prebasecases[is.na(tab$prebasecases)] = 0
##set up outcome based on number of years of follow-up
tab$maxaccumxs = tab$accum1yr
##classification
tab$outcome = ifelse(tab$maxaccumxs <= 0, "cluster", "outbreak")

train = tab[tab$initflag==2,]

##set up predictors
preds = read.table("variablesList_comb_first_2Q_half.txt", header = F, stringsAsFactors = F)
preds = as.character(preds$V1)

##fix training data classes
facVars = c(names(tab)[c(grep("half", names(tab)), grep("any", names(tab)))], "initflag", "incident")
contVars = c(names(tab)[grep("pct", names(tab))], "casecnt", "timediff", "maxaccumxs", "nbhflag", "qtrsbtwflags", "nbhbase", "nbhcount", "prebasecases")
for(c in facVars) {
  tab[,c] = factor(tab[,c], levels=c(1,2))
  train[,c] = factor(train[,c], levels=c(1,2))
}
for(c in contVars) {
  tab[,c] = as.numeric(tab[,c])
  train[,c] = as.numeric(train[,c])
}

##generate and save model (function copied from )
data = train[,c(preds, "outcome")]
data$outcome = factor(data$outcome)
model = randomForest(outcome ~ ., data = data)
save(model,file="modelRF_classification.rda")

##function to get performance statistics (and saves confusion matrix)
##model = model used to get the performance statistics
##preds = the set of predictors used to generate the model
##name = the name of the algorithm/set/group, used for identification in returned table
##data = the data to be used for the predictions and stats
getStats <- function(model, preds, name, data) { 
  pred = predict(model, data[,preds], type="class")

  predClass = factor(pred, levels=c("outbreak", "cluster"))
  actClass = factor(data$outcome, levels=c("outbreak", "cluster"))
  cm = confusionMatrix(data=predClass, reference = actClass, positive = "outbreak")
  sink(paste0("StatsConfusionMatrix_", name, ".txt"))
  print(cm)
  sink()
  
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
  
  comb = data.frame(nbh_no = data$nbh_no,
                    clus_no = data$clus_no,
                    actualClass = data$outcome,
                    predictedClass = pred)
  write.csv(comb, paste0("prediction_", name, ".csv"), row.names = F)
  
  return(res)
}

##performance on second and third flags
res1 = getStats(model, preds, "addlFlags", tab[tab$initflag==1,])

##set up validation set
valid = read.csv(paste0(inputDir, "VFulldata2Q_072920.csv"), header = T, stringsAsFactors = F)
vout = read.csv(paste0(inputDir, "ValidAccum_081420.csv"), header = T, stringsAsFactors = F)
names(valid)[names(valid)=="vaccum1yr"] = "old_vaccum1yr"
valid = merge(valid, vout[,c("nbh_no", "vaccum1yr")], by="nbh_no")

valid$outcome = ifelse(valid$vaccum1yr <= 0, "cluster", "outbreak")
valid$incident = ifelse(valid$clustyp=="Incident", 2, 1)
names(valid) = gsub("prev_", "prevtb_", names(valid))
valid$initflag = ifelse(valid$initflag=="Y", 2, 1)
valid$maxaccumxs = valid$vaccum1yr
for(c in facVars) {
  valid[,c] = factor(valid[,c], levels=c(1,2))
}
for(c in contVars) {
  valid[,c] = as.numeric(valid[,c])
}

##performance on validation set (all flags and separating the flags)
res2 = getStats(model, preds, "validationSet", valid)
res3 = getStats(model, preds, "validationSetFirstFlag", valid[valid$initflag==2,])
res4 = getStats(model, preds, "validationSetAddlFlag", valid[valid$initflag==1,])
write.csv(rbind(res1, res2, res3, res4), "classification_prediction_stats.csv", row.names = F)
write.csv(rbind(res3), "classification_prediction_stats_validfirst.csv", row.names = F)
