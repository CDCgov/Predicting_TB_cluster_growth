## functions formoving hurdle fits on quarterly time series 

library(pscl) # for ZIP and ZINB, hurdle etc, pscl::zeroinfl and pscl::hurdle
            
# functions for ZINB and NB hurdle  

# zero inflated negative binomial (ZINB) 
# MLE 
zinb.mle <- function(x) {
  n0 = sum(x==0)  
  if (n0 == length(x)) {  # all zeros
    mle = c(1,0,0)
  } else if (n0 > 0) {  # some zeros
    zinb.fit = zeroinfl(x~1, data=as.data.frame(x), dist = "negbin")
    p.hat = plogis(zinb.fit$coefficients$zero)
    mu.hat = exp(zinb.fit$coefficients$count)
    size.hat = zinb.fit$theta 
    mle = c(p.hat, mu.hat, size.hat)  
  } else {  # no zeros 
    nb.fit <- glm.nb(x~1, data=as.data.frame(x)) # straight Neg Bin fit if no zeros 
    mu.hat = exp(nb.fit$coef)
    size.hat = nb.fit$theta 
    mle = c(0, mu.hat, size.hat) 
  }
  return(mle)  
}

# 95th percentile 
zinb.u95 <- function(x) { 
  if(x[1] >= .95) { 
    upper=3 # arbitrary, unlikely to happen in the dataset
  } else { upper = qnbinom((.95-x[1])/(1-x[1]), mu = x[2], size = x[3]) 
  }  
  return(upper) 
}

# hurdle neg bin 
# MLE 
hrdlnb.mle <- function(x) { 
  n0 = sum(x==0)  
  if (n0 == length(x)) {  # all zeros
    mle = c(1,0,0)
  } else if (n0 > 0 & n0 < length(x)) {  # some zeros 
    hrdlnb.fit = hurdle(x~1, data=as.data.frame(x), dist = "negbin")
    p.hat = plogis(hrdlnb.fit$coefficients$zero) # back transform from logit scale 
    mu.hat = exp(hrdlnb.fit$coefficients$count) # back transform from log scale 
    size.hat = hrdlnb.fit$theta 
    mle = c(1-p.hat, mu.hat, size.hat)  # hurdle p estimates are backwards (1-p) in pscl...  
  } else {  # no zeros
    nb.fit <- glm.nb(x~1, data=as.data.frame(x)) # straight Neg Bin fit if no zeros 
    mu.hat = exp(nb.fit$coef)
    size.hat = nb.fit$theta 
    mle = c(0, mu.hat, size.hat) 
  }
  return(mle)  
}

# 95th percentile 
hrdlnb.u95 <- function(x) { 
  if(x[1] >= .95) { # arbitrary, unlikely to happen in the dataset
    upper=3 
  } else { upper = qnbinom((.95-x[1])/(1-x[1])*(1-dnbinom(0, mu = x[2], size = x[3])) + dnbinom(0, mu = x[2], size = x[3]), 
                           mu = x[2], size = x[3]) 
  }
  return(upper)
}





