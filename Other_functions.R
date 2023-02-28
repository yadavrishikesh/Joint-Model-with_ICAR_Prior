#rm(list=ls())
### pdf of extednded GPD1
library(evd) ## FOR WRITING THE EGPD DISTRIBUTION 
library(actuar) ## FOR BURR DISTRIBUTION 
dEGPD1<-function(x, k, xi, sigma, log=NULL){
  
  if(log==FALSE){
    dens<- k * (evd::dgpd(x=x, loc=0, scale=sigma, shape=xi)) * (evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE))^(k-1) 
  } else{
    dens<- log(k) + (evd::dgpd(x=x, loc=0, scale=sigma, shape=xi, log = TRUE)) + ((k-1)* log(evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE)))
  }
  return(dens)
  
}

# dEGPD1<-Vectorize(dEGPD1)
# 
# ### Check the function  for k=1, which is same as GPD
#  dEGPD1(x=1, k=1, xi=0.2, sigma = 2, log=TRUE)
# dgpd(x=1:10, loc=0, scale=1:10, shape=0.2, log = TRUE) 
# dEGPD1(x=1:10, k=1, xi=0.2, sigma = 1:10, log=FALSE)
# dgpd(x=1:10, loc=0, scale=1:10, shape=0.2, log = FALSE) 



### cdf of extednded GPD1
pEGPD1<-function(x, k, xi, sigma, log=NULL){
  if(log==FALSE){
    cdf<-  (evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE))^(k) 
  } else{
    cdf<- (k)* log(evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE))
  }
  return(cdf)
}
# pEGPD1(x=1:10, k=1, xi=0.2, sigma = 1:10, log=TRUE)
# log(pgpd(q=1:10, loc=0, scale=1:10, shape=0.2, lower.tail = TRUE)) 


### simulating random number fronm GPD1
rEGPD1<-function(n, k, xi, sigma){
  X<- (sigma/xi) * (((1-(runif(n=n))^(1/k))^(-xi))-1)
  return(X)
}

# set.seed(1)
# sim_egpd1<-rEGPD1(n=10000, k=1, xi=0.8, sigma = 2)
# set.seed(1)
# sim_gpd<-rgpd(n=10000, loc=0, scale = 2, shape=0.8)
# sim_egpd1; sim_gpd
# par(mfrow=c(1,2))
# hist(sim_gpd); hist(sim_egpd1)
# mean(sim_gpd);mean(sim_egpd1)
# median(sim_gpd);median(sim_egpd1)

###############################################################
############## FUNCTION FOR GENEARLIZED GAMMA DISTRIBUTION #############
###############################################################


# theta   : scale parameter
# kappa      : shape parameter
# delta   : shape parameter
# t       : positive argument   
# p       : probability (0,1)
# n       : number of simulations

# Probability Density Function
dggamma <- function(t, theta, kappa, delta, log = FALSE){
  val <- log(delta) - (kappa*log(theta)) - lgamma(kappa/delta) + ((kappa - 1)*log(t)) - ((t/theta)^delta)
  if(log) return(val) else return(exp(val))
}

# GG CDF
pggamma <- function(t, theta, kappa, delta, log.p = FALSE){
  val <- pgamma(t^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE) 
  if(log.p) return(val) else return(exp(val))
}

# GG Survival Function
sggamma <- function(t, theta, kappa, delta, log.p = FALSE){
  val <- pgamma( t^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE, lower.tail =  FALSE) 
  if(log.p) return(val) else return(exp(val))
}


# GG Hazard Function
hggamma <- function(t, theta, kappa, delta, log = FALSE){
  val <- dggamma(t, theta, kappa, delta, log = TRUE) - sggamma(t, theta, kappa, delta, log.p = TRUE)
  if(log) return(val) else return(exp(val))
}

# GG Cumulative Hazard Function
chggamma <- function(t, theta, kappa, delta){
  val <- -pgamma( t^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE, lower.tail =  FALSE) 
  return(val) 
}

# Quantile Function
qggamma <- function(p, theta, kappa, delta){
  out <- theta * (qgamma(p, shape = kappa/delta, rate = 1)^(1/delta))
  # out <- qgamma(u, shape = kappa/delta, scale = theta^delta)^(1/delta)
  return(out)
}

# Random number Generation Function
rggamma <- function(n, theta, kappa, delta){
  u <- runif(n)
 # out <- qgamma(u, shape = kappa/delta, scale = theta^delta)^(1/delta)
  out <- theta * (qgamma(u, shape = kappa/delta, rate = 1)^(1/delta))
  return(as.vector(out))
}


