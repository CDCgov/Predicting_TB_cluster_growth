## script to do moving hurdle fits on quarterly time series 

library(zoo) # for moving averages zoo::rollapply and zoo::as.yearqtr
library(pscl) # for ZIP and ZINB, hurdle etc, pscl::zeroinfl and pscl::hurdle
library(MASS)
## Data Prep
prv.cl <- read.csv("Elig4Enrol_n3.csv", header=TRUE, sep=",") 
names(prv.cl)

prv.cl$gimsdt <- as.yearqtr(prv.cl$gimsdt, format = "%m/%d/%Y") # so R reads as quarterly dates
prv.cl$clus_no <- as.factor(prv.cl$clus_no) # so R reads as categorical
str(prv.cl) # check 
head(prv.cl) 
tail(prv.cl)

# Data fossick
# date range
range(prv.cl$gimsdt) 
# corresponding number of quarters 
( nq <- length(unique(prv.cl$gimsdt)) ) 
# number of clusters
( nc <- nlevels(prv.cl$clus_no) ) 
# number of rows should be the product of these 2 
nrow(prv.cl) == nq*nc

# save orginal dataset and discard clusters of fewer than 3 cases as these won't be flagged 
prv.cl.original <- prv.cl 

# sum cases within each cluster 
sum_by_clus <- lapply( with(prv.cl, split(count, clus_no)) , sum) 
# store cluster IDs which have more than 2 cases in total 
keep_clus <- levels(prv.cl$clus_no)[which(unlist(sum_by_clus) > 2)] 
# number of clusters retained 
keep_clus_n <- length(keep_clus)
# store cluster IDs which have fewer than 3 cases in total 
drop_clus <- levels(prv.cl$clus_no)[which(unlist(sum_by_clus) <= 2)] 

# subset the dataframe 
prv.cl <- prv.cl[prv.cl$clus_no %in% keep_clus, ]
# drop empty levels for which < 3 cases in total
prv.cl$clus_no <- droplevels(prv.cl$clus_no) 

## call the ZINB and hurdle MLE functions 
source("cluster_func.R") # specify path if this file is in diferent directory 

## Calcs  
mv.w <- 8  # width of moving window in quarters (usually a 2 year window)
timeline <- unique(prv.cl$gimsdt)[-(1:mv.w)] # apply moving fits starting after the first full window

start.time <- Sys.time()
# rolling fits 
# initialize outputs   
flag.qtr <- flag.case <- flag.up <- flag.ma <- list()    
# qtr is the timepoint, case is number of cases, up is upper band, ma is moving average 
for (j in 1:nlevels(prv.cl$clus_no)) { # j indexes distinct clusters  
  cases <- prv.cl$count[prv.cl$clus_no == levels(prv.cl$clus_no)[j] ] # extract timeseries for each cluster ID  
  mle <- rollapply(cases, width=mv.w, hrdlnb.mle) 
  mv.up <- apply(mle, MARGIN = 1, FUN = hrdlnb.u95)
  ma <- rollapply(cases, width=mv.w, mean)
  x <- which(cases[-(1:mv.w)] > mv.up[-length(mv.up)]) 
  # record positions at which cases post baseline period onwards exceed upper band
  # if no instances, x will be zero length, returned as integer(0) 
  if ( length(x) > 0 ) {    
    flag.qtr[[j]] <- timeline[x] # timepoint corresponding to position x 
    flag.case[[j]] <- cases[-(1:mv.w)][x] 
    flag.up[[j]] <- mv.up[x] 
    flag.ma[[j]] <- ma[x] 
  }
  txtProgressBar(1, nlevels(prv.cl$clus_no), j, style=3) # progress indicator 
} # may give warnings  
end.time <- Sys.time()
(time.taken <- end.time - start.time)

# put it all together 
# store which clusters were flagged and which weren't 
clusters_flagged <- levels(prv.cl$clus_no)[which(lengths(flag.case) > 0)] 
keep_clus_n_flagged <- length(clusters_flagged) 
clusters_not_flagged <- levels(prv.cl$clus_no)[-which(lengths(flag.case) > 0)]
keep_clus_n_not_flagged <- length(clusters_not_flagged) 

# remove clusters for which there was no flag (NULL entries of list) 
flag.qtr <- Filter(length, flag.qtr)  
flag.case <- Filter(length, flag.case) 
flag.up <- Filter(length, flag.up) 
flag.ma <- Filter(length, flag.ma) 

# results
hrdlnb.results <- data.frame(rep(clusters_flagged, lengths(flag.case)), 
                             as.yearqtr(unlist(flag.qtr)), 
                             unlist(flag.case), 
                             unlist(flag.up), 
                             unlist(flag.ma)
                             )    
names(hrdlnb.results) <- c("cluster", "quarter", "count", "upper", "MA")  

cat("There are", nc, "clusters", 
    "of which", length(drop_clus), "have fewer than 3 cases in total.","\n",
    "Of the remaining", keep_clus_n, "clusters,", keep_clus_n_not_flagged, "were not flagged and", 
    keep_clus_n_flagged, "were flagged.")  


# export to CSV
write.csv(hrdlnb.results, "n3clusters_011118_hrdlNB.csv", row.names = F)  
# save the R data set if desired 
save(hrdlnb.results, file = "en3clusters_011118_hrdlNB.Rda") 


