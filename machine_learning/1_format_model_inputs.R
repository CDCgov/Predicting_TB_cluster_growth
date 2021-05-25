##split test/train dataset for cross-validation, generating the input tables and lists of predictors for the models
set.seed(1234)

##function to generate the list of variables to include in an analysis
##vars = initial list of variables to include (some may be removed if there are few unique values)
##tab = table of data
##name = suffix used in output files
fixVars <- function(vars, tab, name) {
  print(name)
  ##remove variables with only one value (not informative and can mess up some algorithms)
  for(v in vars) {
    if(!(grepl("pct", v) | grepl("cnt", v)) & !v %in% c("timediff", "casecnt", "nbhflag", "qtrsbtwflags", "nbhbase", "nbhcount", "prebasecases")) { #keep continuous variables if have more than 5 unique values
      counts = table(tab[,v])
      if(length(unique(tab[,v])) < 2 | any(counts < 5)) {
        print(paste("removed", v))
        vars = vars[vars!=v]
      }
    } else {
      if(length(unique(tab[,v])) < 2) {
        print(paste("removed", v))
        vars = vars[vars!=v]
      }
    }
  }
  any = vars[grepl("_any", vars)]
  pct = vars[grepl("_pct", vars)]
  cnt = vars[grepl("_cnt", vars)]
  vars = vars[!vars %in% cnt] #remove count
  hlf = vars[grepl("_half", vars)]
  oth = c("timediff", "nbhcount", "nbhbase", "incident", "prebasecases")
  print(vars[!vars %in% c(any,pct,cnt,hlf,oth)]) #should only be empty
  
  write.table(data.frame(k = c(any, oth)), paste0("variablesList_", name, "_", flag, "_", co, "_any.txt"), row.names = F, col.names = F, quote = F)
  write.table(data.frame(k = c(pct, oth)), paste0("variablesList_", name, "_", flag, "_", co, "_pct.txt"), row.names = F, col.names = F, quote = F)
  write.table(data.frame(k = c(hlf, oth)), paste0("variablesList_", name, "_", flag, "_", co, "_half.txt"), row.names = F, col.names = F, quote = F)
  write.table(data.frame(k = vars), paste0("variablesList_", name, "_", flag, "_", co, "_all.txt"), row.names = F, col.names = F, quote = F)
  
  return(vars)
}

###split each cohort into cross-validation datasets
flag = "first" #include first flag only
cohorts = c("1Q", "2Q", "4Q", "9Q")
yr = paste0(1, "yrfu") #follow up period number of years for outcome (years)
for(co in cohorts) {
  tab = read.csv(paste0("../data/Fulldata", tolower(co), "_072920.csv"), header = T, stringsAsFactors = F)
  names(tab) = gsub("prev_", "prevtb_", names(tab)) #these are the same variables but have different starts
  allflags = tab
  
  ##set up outcome based on number of years of follow-up
  tab$maxaccumxs = tab$accum1yr
  ##classification
  tab$outcome = ifelse(tab$maxaccumxs <= 0, "cluster", "outbreak")
  
  ##set up incident/prevalent
  tab$incident = ifelse(tab$clustyp=="Incident", 2, 1)
  tab$prebasecases[is.na(tab$prebasecases)] = 0
  
  name = paste(flag, yr, co, sep="_")
  print(name)
  ###split data for cross validation
  tab = tab[tab$initflag=="Y",]
  n = nrow(tab)/5
  sp5 = 1:nrow(tab) #5th split will be the remainder of the others
  sp1 = sample(x = sp5, size=n+1, replace = F) #remainder of 3 when divide by 5, distribute
  sp5 = sp5[!sp5 %in% sp1]
  sp2 = sample(x = sp5, size=n+1, replace = F)
  sp5 = sp5[!sp5 %in% sp2]
  sp3 = sample(x = sp5, size=n, replace = F)
  sp5 = sp5[!sp5 %in% sp3]
  sp4 = sample(x = sp5, size=n, replace = F)
  sp5 = sp5[!sp5 %in% sp4]
  
  write.csv(tab[c(sp2, sp3, sp4, sp5),], paste0("trainingSet1C_", name, ".csv"), row.names = F)
  write.csv(tab[sp1,], paste0("testSet1C_", name, ".csv"), row.names = F)
  
  write.csv(tab[c(sp1, sp3, sp4, sp5),], paste0("trainingSet2C_", name, ".csv"), row.names = F)
  write.csv(tab[sp2,], paste0("testSet2C_", name, ".csv"), row.names = F)
  
  write.csv(tab[c(sp1, sp2, sp4, sp5),], paste0("trainingSet3C_", name, ".csv"), row.names = F)
  write.csv(tab[sp3,], paste0("testSet3C_", name, ".csv"), row.names = F)
  
  write.csv(tab[c(sp1, sp2, sp3, sp5),], paste0("trainingSet4C_", name, ".csv"), row.names = F)
  write.csv(tab[sp4,], paste0("testSet4C_", name, ".csv"), row.names = F)
  
  write.csv(tab[c(sp1, sp2, sp3, sp4),], paste0("trainingSet5C_", name, ".csv"), row.names = F)
  write.csv(tab[sp5,], paste0("testSet5C_", name, ".csv"), row.names = F)
  
  ####check splits
  for(r in 1:5) {
    train = read.csv(paste0("trainingSet", r, "C_", name, ".csv"), header = T, stringsAsFactors = F)
    test = read.csv(paste0("testSet", r, "C_", name, ".csv"), header = T, stringsAsFactors = F)
    if(nrow(train) + nrow(test) != nrow(tab)) {
      stop("Incorrect number of rows")
    }
    clust = c(train$nbh_no, test$nbh_no)
    if(any(duplicated(clust))) {
      stop("Duplicated IDs")
    }
    if(!all(tab$nbh_no %in% clust)) {
      stop("Missing IDs")
    }
  }
  
  ####set up list of variables
  vars = c("nbhbase", "nbhcount", "timediff", "incident", "prebasecases",
           names(tab)[which(names(tab)=="hl_pct"):which(names(tab)=="margin_half")])
  ##remove drug, substance, and margin
  vars = vars[!grepl("subs", vars) & !grepl("margin", vars) & !grepl("drug", vars)]
  vc = fixVars(vars, tab, "comb")
}
