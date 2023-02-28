######################################################################################
####### log-posterior for the hyperparameters in the Mark distributiobns##########################
######################################################################################
s.fixed<-1
### Imputing NA for the landslides counts
impute.NA.Y<-function(ind_NA_Y, eta, CV){
  n1<-length(eta)
  if(CV=="OSD"){
  imputed_NA_Y<-rpois(sum(ind_NA_Y), lambda = exp(eta)[ind_NA_Y])
  } else if(CV=="WSD"){
    imputed_NA_Y<- rpois(n1, lambda = exp(eta))
  }
  return(imputed_NA_Y)
}

### Imputing NA for landslides sizes
impute.NA.A<-function(ind_NA_A, mu, log.hyper.mu, family, CV){
  n2<-length(mu)
  ###### Gamma-gamma model
  if(family=="Gamma-gamma"){
    if(CV=="OSD"){
    gamma1<- exp(log.hyper.mu[1])
    gamma2<- exp(log.hyper.mu[2])
    alpha<- exp(mu)/qf(0.5, df1=gamma1, df2=gamma2)
    imputed_NA_A<- alpha[ind_NA_A]* rf(n=sum(ind_NA_A), df1=gamma1, df2=gamma2)
    } else if(CV=="WSD"){
      gamma1<- exp(log.hyper.mu[1])
      gamma2<- exp(log.hyper.mu[2])
      alpha<- exp(mu)/qf(0.5, df1=gamma1, df2=gamma2)
      imputed_NA_A<- alpha * rf(n=n2, df1=gamma1, df2=gamma2)
    }
    ###### extended GPD
  } else if (family=="extended-GPD"){
    if(CV=="OSD"){
    k<-exp(log.hyper.mu[1])
    xi<-exp(log.hyper.mu[2])
   # xi<- 2 * (exp(log.hyper.mu[2])-1) / (exp(log.hyper.mu[2])+1)
    # if(k < 0.2){
    #   k<-0.2
    # } else{
    #   k<- k
    # }
    # if(xi> 0.001){
    #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
    # } else{
    #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
    # }
    sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    imputed_NA_A<- rEGPD1(n=sum(ind_NA_A), k=k, xi=xi, sigma=sigma[ind_NA_A])
    } else if(CV=="WSD"){
      k<-exp(log.hyper.mu[1])
      xi<-exp(log.hyper.mu[2])
      sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
      imputed_NA_A<- rEGPD1(n=n2, k=k, xi=xi, sigma=sigma)
    }
      
   ##### GPD   
  } else if (family=="GPD"){
    if(CV=="OSD"){
    k<- 1
    xi<- exp(log.hyper.mu)
    #xi<- 2 * (exp(log.hyper.mu)-1) / (exp(log.hyper.mu)+1)
    # if(xi > 0.001){
    #   sigma<- exp(mu) * xi / ((0.5)^(-xi)-1)
    # } else{
    #   sigma<- - exp(mu)/ log(0.5)
    # }
    sigma<- exp(mu) / evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    #imputed_NA_A<- evd::rgpd(n=sum(ind_NA_A), loc=0, shape =xi, scale =sigma[ind_NA_A])
    imputed_NA_A<- rEGPD1(n=sum(ind_NA_A), k=k, xi=xi, sigma=sigma[ind_NA_A])
    } else if(CV=="WSD"){
      k<- 1
      xi<- exp(log.hyper.mu)
      sigma<- exp(mu) / evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
      #imputed_NA_A<- evd::rgpd(n=sum(ind_NA_A), loc=0, shape =xi, scale =sigma[ind_NA_A])
      imputed_NA_A<- rEGPD1(n=n2, k=k, xi=xi, sigma=sigma)
    }
    ##### Gamma
  } else if(family=="Gamma"){
    if(CV=="OSD"){
    # kappa_a<-exp(log.hyper.mu)
    # imputed_NA_A<-rgamma(n=sum(ind_NA_A), shape = s.fixed* kappa_a, rate = s.fixed * kappa_a* exp(-mu)[ind_NA_A])

    k<-exp(log.hyper.mu)
    c<-1
    #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) # mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    imputed_NA_A<- rggamma(n=sum(ind_NA_A), theta = sigma[ind_NA_A],  kappa=k, delta=c)
    } else if(CV=="WSD"){
      k<-exp(log.hyper.mu)
      c<-1
      #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) # mean
      sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
      imputed_NA_A<- rggamma(n=n2, theta = sigma,  kappa=k, delta=c)
    }
    
    ### Log-Gamma
  } else if (family=="log-Gamma"){
    if(CV=="OSD"){
    kappa_a<-exp(log.hyper.mu)
    #rate.a<-1/(1-((1+exp(mu))^(-1/kappa_a)))
    # if(kappa_a>1){
    #   rate.a<- ((kappa_a-1)/(log(1+exp(mu))))-1  #Covariates in mode
    # } else {
    #   rate.a<- 10^10
    # }
    rate.a<- qgamma(0.5, shape = kappa_a, rate = 1)/log(1+exp(mu))
    imputed_NA_A<-exp(rgamma(n=sum(ind_NA_A), shape = kappa_a, rate = rate.a[ind_NA_A]))-1
    } else if(CV=="WSD"){
      kappa_a<- exp(log.hyper.mu)
      rate.a<- qgamma(0.5, shape = kappa_a, rate = 1)/log(1+exp(mu))
      imputed_NA_A<- exp(rgamma(n=n2, shape = kappa_a, rate = rate.a))-1
    }
    ### Weibull 
    
  } else if(family=="Weibull"){
    if(CV=="OSD"){
    # kappa_a<-exp(log.hyper.mu)
    # lambda<-exp(mu)/gamma(1+kappa_a)
    # imputed_NA_A<-rweibull(n=sum(ind_NA_A), shape = kappa_a, scale = lambda[ind_NA_A])
    k<- exp(log.hyper.mu)
    c<- k
    #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) # mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    imputed_NA_A<- rggamma(n=sum(ind_NA_A), theta = sigma[ind_NA_A],  kappa=k, delta=c)
    } else if(CV=="WSD"){
      k<- exp(log.hyper.mu)
      c<- k
      #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) # mean
      sigma<-  exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
      imputed_NA_A<- rggamma(n=n2, theta = sigma,  kappa=k, delta=c)
    }
 #### log-Normal
  } else if(family=="log-normal"){
    if(CV=="OSD"){
    kappa_a<-exp(log.hyper.mu)
    imputed_NA_A<-exp(rnorm(n=sum(ind_NA_A), mean =mu[ind_NA_A], sd=sqrt(1/kappa_a)))
    } else if(CV=="WSD"){
      kappa_a<-exp(log.hyper.mu)
      imputed_NA_A<-exp(rnorm(n=n2, mean =mu, sd=sqrt(1/kappa_a)))
    }
  ##### Burr
  } else if(family=="Burr"){
    if(CV=="OSD"){
    k<-exp(log.hyper.mu[1])
    c<-exp(log.hyper.mu[2])
    #sigma<-exp(mu)/(((2^(1/k))-1)^(1/c))
    sigma<- exp(mu) /qburr(0.5, shape1 = k, shape2 = c, scale = 1)
    imputed_NA_A<-rburr(n=sum(ind_NA_A), shape1 = k, shape2 = c, scale = sigma[ind_NA_A])
    } else if(CV=="WSD"){
      k<-exp(log.hyper.mu[1])
      c<-exp(log.hyper.mu[2])
      #sigma<-exp(mu)/(((2^(1/k))-1)^(1/c))
      sigma<- exp(mu) /qburr(0.5, shape1 = k, shape2 = c, scale = 1)
      imputed_NA_A<-rburr(n=n2, shape1 = k, shape2 = c, scale = sigma)
    }
    #### generalized-gamma distribution
  }  else if(family=="generalized-Gamma"){
    if(CV=="OSD"){
    k<-exp(log.hyper.mu[1])
    c<-exp(log.hyper.mu[2])
    #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) # mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    imputed_NA_A<- rggamma(n=sum(ind_NA_A), theta = sigma[ind_NA_A],  kappa=k, delta=c)
    } else if(CV=="WSD"){
      k<-exp(log.hyper.mu[1])
      c<-exp(log.hyper.mu[2])
      #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) # mean
      sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
      imputed_NA_A<- rggamma(n=n2, theta = sigma,  kappa=k, delta=c)
    }
  }
  return(imputed_NA_A)
}
