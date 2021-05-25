##combine model results, average and organize for easier sorting and analysis, and make clean table for publication

setwd("results")

library(xlsx)

flag = "first" #use the first flag only
cohorts = c("1Q", "2Q", "4Q", "9Q") #time frames to iterate over
yr = paste0(1, "yrfu") #follow up period number of years for outcome
fs = "noFS" #no feature selection
repeats = 1:5 #number of repeats in cross validation
a = "C_" #analytic prefix

##read all results into res
res = NA 
for(co in cohorts) {
  for(r in repeats) {
    group = paste0(a, yr, "_", co, "_", flag, "_", fs, "_Set", r)
    print(group)
    fname = paste0("concatenatedModelStats_", group, ".csv")
    if(file.exists(fname)) {
      tab = read.csv(fname, header = T, stringsAsFactors = F)
      if(all(is.na(res))) {
        res = tab
      } else {
        res = rbind(res, tab)
      }
    } else {
      print(paste("missing", fname))
    }    
  }
}

##split alg_set_group for easier sorting
res$method = NA
res$featureSelection = NA
res$variableSet = NA
res$timeframe = NA
res$set = NA
res$flag = NA
for(r in 1:nrow(res)) {
  ex = regexpr("[1-9]Q", res$alg_set_group[r])
  res$timeframe[r] = substr(res$alg_set_group[r], ex[1], ex[1] + attr(ex, "match.length")-1)
  
  ex = regexpr("_[no]*FS_", res$alg_set_group[r])
  res$featureSelection[r] = substr(res$alg_set_group[r], ex[1]+1, ex[1] + attr(ex, "match.length")-2)
  
  res$flag[r] = flag
  
  ex = regexpr("Set[1-5]", res$alg_set_group[r])
  index = ex[1] + attr(ex, "match.length")-1
  res$set[r] = substr(res$alg_set_group[r], index, index)
  
  ex = regexpr("_[a-z]*$", res$alg_set_group[r])
  res$variableSet[r] = substr(res$alg_set_group[r], ex[1]+1, ex[1] + attr(ex, "match.length")-1)
  
  ex = regexpr("Set[0-9]_[a-zA-Z]*_", res$alg_set_group[r])
  ml = substr(res$alg_set_group[r], ex[1]+5, ex[1] + attr(ex, "match.length")-2)
  if(ml=="nn") {
    ex = regexpr("FS_[a-z]*_[1-9]*_", res$alg_set_group[r])
    ml = substr(res$alg_set_group[r], ex[1]+3, ex[1] + attr(ex, "match.length")-2)
    ml = gsub("nn", "neural_network", ml)
  }
  res$method[r] = ml
}
write.csv(res, paste0("combinedClassificationModelResults_", a, "_", yr, "_", flag, ".csv"), row.names = F, na = "")

##calculate average of the 5 sets
len = length(unique(res$timeframe)) * length(unique(res$variableSet)) * length(unique(res$method)) * length(unique(res$featureSelection))

ave = data.frame(cohort_set_method = rep(NA, len),
                 featureSelection = rep(NA, len),
                 timeframe = rep(NA, len),
                 flag = rep(NA, len),
                 variableSet = rep(NA, len),
                 method = rep(NA, len),
                 numRuns = rep(NA, len))
tmp = as.data.frame(matrix(NA, nrow = len, ncol = 2*sum(grepl("test", names(res)))))
metrics = names(res)[grepl("test", names(res))]
names(tmp) = c(metrics, paste0(metrics, ".sd"))
ave = cbind(ave, tmp)
variableSets = sort(unique(res$variableSet))
methods = sort(unique(res$method))
index = 1
for(c in cohorts) {
  for(vs in variableSets) {
    for(m in methods) {
      sub = res[res$timeframe == c & res$flag==flag & res$variableSet==vs & res$method==m & res$featureSelection==fs ,]
      sub[is.na(sub)] = 0 #sometimes one is missing
      if(nrow(sub) > 0) {
        ave$cohort_set_method[index] = paste(c, flag, vs, m, sep="_")
        ave$featureSelection[index] = fs
        ave$timeframe[index] = c
        ave$flag[index] = flag
        ave$variableSet[index] = vs
        ave$method[index] = m
        ave$numRuns[index] = nrow(sub)
        for(met in metrics) {
          ave[index, met] = mean(sub[,met])
          ave[index, paste0(met, ".sd")] = sd(sub[,met])
        }
        index = index + 1
      }
    }
  }
}

write.csv(ave, paste0("averageClassificationModelResults_", a, "_", yr, "_", flag, ".csv"), row.names = F, na = "")

####make table for publication
##function that for the given method abbreviation used in code return name to be used in table
methodNames <- function(method) {
  if(method=="cforest") {
    return("random forest (party package)")
  } else if(method=="ctree") {
    return("CART (party package)")
  } else if(method=="elasticNet") {
    return("elastic net")
  } else if(method=="gbm") {
    return("GBM") 
  } else if(method=="lassoRegression") {
    return("Lasso regression")
  } else if(method=="randomForest") {
    return("random forest")
  } else if(method=="ridgeRegression") {
    return("Ridge regression")
  } else if(method=="rpart") {
    return("CART (party package)")
  } else if(method=="svmlinear") {
    return("SVM with linear kernel")
  } else if(method=="svmpoly") {
    return("SVM with polynomial kernel")
  } else if(method=="svmsigmoid") {
    return("SVM with sigmoid kernel")
  } else if(method=="svmradial") {
    return("SVM with radial kernel")
  } else if(method=="wsvmlinear") {
    return("weighted SVM with linear kernel")
  } else if(method=="wsvmpoly") {
    return("weighted SVM with polynomial kernel")
  } else if(method=="wsvmradial") {
    return("weighted SVM with radial kernel")
  } else if(method=="wsvmsigmoid") {
    return("weighted SVM with sigmoid kernel")
  } else if(method=="xgboost") {
    return("XGBoost")
  } else {
    warning(paste("Missing method:", method))
    return(method)
  }
}

##write Excel table
writeExcelTable <- function(fileName, df, sheetName = "Sheet1", wrapHeader=T) {
  workbook = createWorkbook()
  sheet = createSheet(workbook, sheetName = sheetName)
  rows = createRow(sheet, rowIndex = 1:(nrow(df)+1))
  cells = createCell(rows, colIndex = 1:ncol(df))
  
  ##header
  headstyle = CellStyle(workbook) + Fill(foregroundColor = rgb(red = 192, green = 192, blue = 192, maxColorValue = 255))
  if(wrapHeader) {
    headstyle = headstyle + Alignment(wrapText = T)
  }
  for(c in 1:ncol(df)) {
    setCellValue(cells[[1,c]], names(df)[c])
    setCellStyle(cells[[1,c]], headstyle)
  }
  ##write data
  for(r in 1:nrow(df)) {
    for(c in 1:ncol(df)) {
      if(!is.na(df[r,c])) {
        if(as.character(df[r,c])!="") {
          if(all(grepl("[0-9.]", strsplit(as.character(df[r,c]), "")[[1]]))) {
            setCellValue(cells[[r+1,c]], as.numeric(as.character(df[r,c])))
          } else {
            setCellValue(cells[[r+1,c]], df[r,c])
          }
        }
      }
    }
  }
  ##set column widths
  autoSizeColumn(sheet, 1:ncol(df))
  saveWorkbook(workbook, fileName)
}

##clean averaged classification table
cleanAveClassification <- function(inFile, outFile, metrics) {
  class = read.csv(inFile, header = T, stringsAsFactors = F)
  
  ##clean model names
  class$method = sapply(class$method, methodNames, USE.NAMES = F)
  ##feature selection
  class = class[class$featureSelection=="noFS",]
  ##remove extra columns
  class = class[,!names(class) %in% c("cohort_set_method", "numRuns", "AUC.test", "flag")]
  
  ##sort by Youden index
  class = class[order(class$youden.test, decreasing = T),]
  
  ##round the numbers
  for(c in names(class)[grepl(".test", names(class)) | grepl("PPV", names(class))]) {
    class[,c] = round(class[,c], digits = 3)
  }
  ##clean header
  keep = c("timeframe", "variableSet", "method")
  for(m in metrics) {
    keep = c(keep, paste0(m, ".test"), paste0(m, ".test.sd"))
  }
  keep = gsub("balanced.accuracy", "balancedAccuracy", keep)
  class = class[,keep]
  names(class) = gsub(".test.sd", " standard deviation", names(class))
  names(class) = gsub(".test", " mean", names(class))
  names(class) = gsub("balancedAccuracy", "balanced accuracy", names(class))
  names(class) = gsub("youden", "Youden index", names(class))
  names(class)[names(class)=="variableSet"] = "quantification"
  ##write
  if(file.exists(outFile)) {
    file.remove(outFile)
  }
  writeExcelTable(df = class, fileName = outFile)
}


metrics = c("youden", "accuracy", "sensitivity", "PPV", "specificity", "NPV")
cleanAveClassification(inFile = "averageClassificationModelResults_C__1yrfu_first.csv",
                       outFile = "AverageClassificationResults_casedate_1yrfu_firstflag.xlsx",
                       metrics = metrics)

###make a nice table (not sortable) for paper
cleanForPaper<-function(inFile, outFile, metrics) {
  class = read.xlsx(inFile, sheetIndex = 1)
  
  ##remove low sensitivity
  class = class[class$sensitivity.mean > 0.25,]
  
  ##combine mean and sd to mean(SD)
  for(m in metrics) {
    mean = format(class[,paste0(m, ".mean")], nsmall=2)
    class$tmp = paste0(mean, "(", class[,paste0(m, ".standard.deviation")], ")")
    class = class[,!grepl(paste0("^", m), names(class))]
    names(class)[names(class)=="tmp"] = m
  }
  
  ##remove feature selection
  class = class[,!names(class)=="feature.selection"]
  
  
  ##clean header
  names(class) = gsub(".", " ", names(class), fixed = T)
  
  ##write
  if(file.exists(outFile)) {
    file.remove(outFile)
  }
  writeExcelTable(df = class, fileName = outFile)
}
cleanForPaper(inFile = "AverageClassificationResults_casedate_1yrfu_firstflag.xlsx",
              outFile = "AverageClassificationResults_casedate_1yrfu_firstflag_forPaper.xlsx",
              metrics = c("Youden.index", "accuracy", "sensitivity", "PPV",  "specificity", "NPV"))

##subset given table to just the best results
cleanForPaperTop<-function(inFile, outFile) {
  class = read.xlsx(inFile, sheetIndex = 1)
  
  ##subset to first for each method
  class = class[!duplicated(class$method),]
  
  ##write
  if(file.exists(outFile)) {
    file.remove(outFile)
  }
  writeExcelTable(df = class, fileName = outFile)
}
cleanForPaperTop(inFile = "AverageClassificationResults_casedate_1yrfu_firstflag_forPaper.xlsx",
                 outFile = "AverageClassificationResults_casedate_1yrfu_firstflag_forPaper_top.xlsx")
