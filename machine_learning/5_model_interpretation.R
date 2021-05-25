##generate ALE and importance plots
library(iml)
library(randomForest)
library(ggplot2)
library(gridExtra)
library(xlsx)

set.seed(1234)

setwd("C:\\Users\\nrf1\\Desktop\\github_repos\\Predicting_TB_cluster_growth\\machine_learning\\results")
load("modelRF_classification.rda")

##set up training data
tab = read.csv("../data/Fulldata2Q_072920.csv", header = T, stringsAsFactors = F)
names(tab) = gsub("prev_", "prevtb_", names(tab)) #these are the same variables but have different starts
tab$initflag = ifelse(tab$initflag=="Y", 2, 1)
tab$incident = ifelse(tab$clustyp=="Incident", 2, 1)
tab$prebasecases[is.na(tab$prebasecases)] = 0
##set up outcome based on number of years of follow-up
tab$maxaccumxs = tab$accum1yr
##classification
tab$outcome = ifelse(tab$maxaccumxs <= 0, "cluster", "outbreak")
###factor
facVars = c(names(tab)[c(grep("half", names(tab)), grep("any", names(tab)))], "initflag", "incident")
contVars = c(names(tab)[grep("pct", names(tab))], "casecnt", "timediff", "maxaccumxs", "nbhflag", "qtrsbtwflags", "nbhbase", "nbhcount", "prebasecases")
for(c in facVars) {
  tab[,c] = factor(tab[,c], levels=c(1,2))
}
for(c in contVars) {
  tab[,c] = as.numeric(tab[,c])
}

train = tab[tab$initflag==2,]

##get predictors
preds = read.table("variablesList_comb_first_2Q_half.txt")
preds = as.character(preds$V1)

##fix variable names
labels = read.xlsx("../data/variable_labels.xlsx", sheetIndex = 1)
labels = labels[!is.na(labels$variable.name),]

##fix labels from any to half
labels$variable.name = as.character(labels$variable.name)
labels$label = as.character(labels$label)
labels$variable.name = gsub("_any", "_half", labels$variable.name)
labels$label = gsub("any ", "at least half ", labels$label)

##function that returns the label for the given variable var
getLabel <- function(var) {
  var = gsub("2", "", var)
  if(var %in% labels$variable.name) {
    return(labels$label[labels$variable.name==var])
  } else {
    print(var)
    return(var)
  }
}

##create model with updated labels
data = train[,c(preds, "outcome")]
names(data) = sapply(names(data), getLabel)
imp = model$importance
preds = row.names(imp)[order(imp, decreasing = T)] #get in order of importance
preds = sapply(preds, getLabel)
data = data[,c(preds, "outcome")]
names(data) = gsub("-| ", "_", names(data))
preds = gsub("-| ", "_", preds)
data$outcome = factor(data$outcome)
model = randomForest(outcome ~ ., data = data)

predictor = Predictor$new(model, data = data[,preds], y = data$outcome, class="outbreak") #without the class option, will print both outbreak and cluster

###feature effects
##ale
effs = FeatureEffects$new(predictor, method="ale")
df = effs$results
# jpeg("alePlots_all_iml.jpg", height=1200, width=900)
# plot(effs, ncols=3)
# dev.off()

##default has missing bars in the bar plot so plot each predictor separately with fixed limits
plist = vector('list', length(preds))
for(i in 1:length(preds)){
  e = FeatureEffect$new(predictor = predictor, feature = preds[i])
  p <- plot(e, ylim=c(-0.4, 0.4)) +
    labs(y="ALE")
  if(names(preds[i]) %in% facVars) {
    res = e$results
    p <- p +
      geom_text(label=format(res$.value, digits=2), vjust=ifelse(res$.value < 0, 1, 0))
  }
  plist[[i]] <- p
}
# jpeg("alePlots_all_iml_sep.jpg", height=1200, width=1000) -> used in paper so increase res
jpeg("alePlots_all_iml_sep.jpg", height=4000, width=4000, res = 300)
do.call("grid.arrange", c(plist, ncol=3, top="Accumulated Local Effects (ALE) for Excess Growth"))
dev.off()

##set up for other options
preds = read.table("variablesList_comb_first_2Q_half.txt")
preds = as.character(preds$V1)
load("modelRF_classification.rda")
predictor = Predictor$new(model, data = train[,preds], y = train$outcome) 
effs = FeatureEffects$new(predictor, method="ale")
df = effs$results
# jpeg("alePlots_all_iml_default.jpg", height=1000, width=1000)
# plot(effs)
# dev.off()

###feature importance for RF (do not use default so can add variable type to plot)
imp = model$importance #importance(model, type=2)
# jpeg("importancePlot.jpg", height=800, width=800)
# varImpPlot(model,type=2, labels = sapply(row.names(imp)[order(imp)], getLabel))
# dev.off()

ip = data.frame(variable = sapply(row.names(imp), getLabel),
                gini = imp)
ip = ip[order(ip$MeanDecreaseGini, decreasing = T),]
ip$variable = factor(ip$variable, levels = rev(ip$variable))
clust.char = c("timediff", "nbhbase", "prebasecases", "nbhcount", "incident")
ip$variableType = ifelse(row.names(ip) %in% clust.char,
                         "cluster characteristic", "patient characteristic")

jpeg("importancePlot_category.jpg", height=1500, width=2000, res = 300) #used in paper
ggplot(ip, aes(x=MeanDecreaseGini, y=variable, fill = variableType)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_manual(values = c("black", "white")) +
  theme_bw() +
  labs(fill = "",
       x = "Mean Decrease in Gini Index",
       y = "") +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(face = ifelse(rev(row.names(ip)) %in% clust.char, "bold", "plain")))
dev.off()

