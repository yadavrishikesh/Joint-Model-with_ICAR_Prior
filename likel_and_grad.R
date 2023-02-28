######################################################################################
####### log-posterior for the hyperparameters in the Mark distributiobns##########################
######################################################################################
log_post_hyper.mu<-function(A, mu, log.hyper.mu, hyper.mu_fixed, s.fixed, family){
  n2<-length(A)
  if(family=="Gamma-gamma"){
    gamma1<- exp(log.hyper.mu[1])
    gamma2<- exp(log.hyper.mu[2])
    alpha<- exp(mu) / qf(0.5, df1=gamma1, df2=gamma2)
    log_post <- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu) + sum(df(x=A/alpha, df1=gamma1, df2=gamma2, log = TRUE))-sum(log(alpha)) 
  } else if (family=="extended-GPD"){
    k<- exp(log.hyper.mu[1])
    xi<- exp(log.hyper.mu[2])
   #xi<- 2 * (exp(log.hyper.mu[2])) / (exp(log.hyper.mu[2])+1)
    # if(k < 0.2){
    #   k<-0.2
    # } else{
    #   k<- k
    # }
    # if(abs(xi)> 0.001){
    #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
    # } else{
    #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
    # }
    sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
   # log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) +  sum(log.hyper.mu) + sum(dEGPD1(x=A, k=k, xi=xi, sigma = sigma, log = TRUE)) 
     log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) +  sum(log.hyper.mu) + 
                 n2 * log(k) - sum(log(sigma)) + (k-1) * sum(log(evd::pgpd(A/sigma, loc = 0, scale = 1, shape = xi))) + sum(evd::dgpd(x=A/sigma, loc = 0, scale = 1, shape = xi, log = TRUE)) 
    
  } else if (family=="GPD"){
    k<- 1
    xi<- exp(log.hyper.mu)
    #xi<- 2 * (exp(log.hyper.mu)) / (exp(log.hyper.mu)+1)
    # if(abs(xi)> 0.001){
    #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
    # } else{
    #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
    # }
    sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    #log_post <- sum(dgamma(xi, shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu)  + sum(dEGPD1(x=A, k=k, xi=xi, sigma = sigma, log = TRUE)) 
    log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) +  sum(log.hyper.mu) + 
      n2 * log(k) - sum(log(sigma)) + (k-1) * sum(log(evd::pgpd(A/sigma, loc = 0, scale = 1, shape = xi))) + sum(evd::dgpd(x=A/sigma, loc = 0, scale = 1, shape = xi, log = TRUE)) 
    
    } else if(family=="Gamma"){
    # kappa_a<-exp(log.hyper.mu)
    # log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu) + sum(dgamma(A, shape=s.fixed*kappa_a, rate = s.fixed * kappa_a * exp(-mu), log=TRUE)) 
    k<- exp(log.hyper.mu)
    c<- 1
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu) + sum(dggamma(t=A, theta = sigma, kappa = k, delta = c, log = TRUE))  
     
  } else if (family=="log-Gamma"){
    kappa_a<- exp(log.hyper.mu)
    #rate.a<- 1/(1-((1+exp(mu))^(-1/kappa_a)))  #Covariates in mean, bt mean does not exit for rate <1
    # if(kappa_a>1){
    # rate.a<- ((kappa_a-1)/(log(1+exp(mu))))-1  #Covariates in mode
    # } else {
    #   rate.a<- 10^10
    # }
    sigma<- log(1+exp(mu)) / qgamma(0.5, shape = kappa_a, scale = 1) ## in median
    log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) +  sum(log.hyper.mu) + sum(dgamma(log(A+1), shape = kappa_a, scale = sigma, log=TRUE)) - sum(log(A+1)) 
    
  } else if(family=="Weibull"){
    # kappa_a<- exp(log.hyper.mu)
    # lambda<- exp(mu)/(log(2)^(1/kappa_a))  # in median
    # #lambda<- exp(mu)/gamma(1+(1/kappa_a))  ## in mean
    # log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu) + sum(dweibull(A, shape=kappa_a, scale = lambda, log=TRUE)) 
    
    k<- exp(log.hyper.mu)
    c<- k
    # sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu) + sum(dggamma(t=A, theta = sigma, kappa = k, delta = c, log = TRUE))  
  
    
    }  else if(family=="Burr"){
    k<- exp(log.hyper.mu[1])
    c<- exp(log.hyper.mu[2])
    #sigma<- exp(mu)/(((2^(1/k))-1)^(1/c))
    sigma<- exp(mu) / qburr(0.5, shape1 = k, shape2 = c, scale = 1)
    log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu) + sum(dburr(x=A, shape1 = k, shape2 = c, scale = sigma, log = TRUE)) 

  } else if(family=="generalized-Gamma"){
    k<-exp(log.hyper.mu[1])
    c<-exp(log.hyper.mu[2])
    # sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    log_post<- sum(dgamma(exp(log.hyper.mu), shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)) + sum(log.hyper.mu) + sum(dggamma(t=A, theta = sigma, kappa = k, delta = c, log = TRUE))  
  }
  return(log_post)
}
########################################################################################################
################ log-posetries that depends on mu
########################################################################################################
log_post_mu<-function(A, mu, intercept2, W1, W2, beta, beta2,  kappa_mu, log.hyper.mu, Z2, A2, s.fixed, family){
  n2<-length(A)
  log.post.process_mu<- -0.5*kappa_mu*sum(mu^2) + kappa_mu * sum(mu *(intercept2+Z2%*%beta2+beta*(A2%*%W1)+A2%*%W2))
  if(family=="Gamma-gamma"){
    gamma1<- exp(log.hyper.mu[1])
    gamma2<- exp(log.hyper.mu[2])
    alpha<- exp(mu)/qf(0.5, df1=gamma1, df2=gamma2)
    log.post_A<- sum(df(A/alpha, df1=gamma1, df2=gamma2, log = TRUE))- sum(log(alpha)) 
  } else if(family=="extended-GPD"){
    k<- exp(log.hyper.mu[1])
    xi<- exp(log.hyper.mu[2])
    #xi<-2 * (exp(log.hyper.mu[2])) / (exp(log.hyper.mu[2])+1)
    # if(k < 0.2){
    #   k=0.2
    # } else{
    #   k<- k
    # }
    # if (xi > 0.001) {
    #   sigma <- exp(mu) * xi / (((1 - ((0.5) ^ (1 / k))) ^ (-xi)) - 1)
    # } else{
    #   sigma <- -exp(mu) / (log(1 - ((0.5) ^ (1 / k))))
    # }
    sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    #log.post_A<-sum(dEGPD1(x=A, k=k, xi=xi, sigma = sigma, log=TRUE)) 
    log.post_A<- n2 * log(k) - sum(log(sigma)) + (k-1) * sum(log(evd::pgpd(A/sigma, loc = 0, scale = 1, shape = xi))) + sum(evd::dgpd(x=A/sigma, loc = 0, scale = 1, shape = xi, log = TRUE)) 
    
    
  } else if(family=="GPD"){
    k<-1
    xi<- exp(log.hyper.mu)
    #xi<- 2 * (exp(log.hyper.mu)) / (exp(log.hyper.mu)+1)
    # if(abs(xi)> 0.001){
    #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
    # } else{
    #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
    # }
    sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    #log.post_A<-  sum(dEGPD1(x=A, k=k, xi=xi, sigma = sigma, log = TRUE)) 
    log.post_A<- n2 * log(k) - sum(log(sigma)) + (k-1) * sum(log(evd::pgpd(A/sigma, loc = 0, scale = 1, shape = xi))) + sum(evd::dgpd(x=A/sigma, loc = 0, scale = 1, shape = xi, log = TRUE)) 
    
  } else if(family=="Gamma"){
    # kappa_a<- exp(log.hyper.mu)
    # log.post_A<- sum(dgamma(A, shape=s.fixed*kappa_a, rate = s.fixed * kappa_a * exp(-mu), log = TRUE)) 
    
    k<- exp(log.hyper.mu)
    c<- 1
    #sigma<- exp(mu) * gamma(k/c) / gamma((k+1)/c)
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    log.post_A<- sum(dggamma(t=A, theta = sigma, kappa = k, delta = c, log = TRUE)) 
    
  } else if(family=="log-Gamma"){
    kappa_a<- exp(log.hyper.mu)
   # rate.a<- 1/(1-((1+exp(mu))^(-1/kappa_a)))
    # if(kappa_a>1){
    #   rate.a<- ((kappa_a-1)/(log(1+exp(mu))))-1  #Covariates in mode
    # } else {
    #   rate.a<- 10^10
    # }
    sigma<- log(1+exp(mu)) / qgamma(0.5, shape = kappa_a, scale = 1) ## in median
    log.post_A<- sum(dgamma(log(A+1), shape=kappa_a, scale  = sigma, log = TRUE)) - sum(log(A+1)) 
  }  else if(family=="Weibull"){
    # kappa_a<- exp(log.hyper.mu)
    # #lambda<-exp(mu)/gamma(1+(1/kappa_a))
    # lambda<- exp(mu)/(log(2)^(1/kappa_a))  # in median
    # log.post_A<- sum(dweibull(A, shape=kappa_a, scale =lambda, log = TRUE)) 
    
    k<- exp(log.hyper.mu)
    c<- k
    #sigma<- exp(mu) * gamma(k/c) / gamma((k+1)/c)
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    log.post_A<- sum(dggamma(t=A, theta = sigma, kappa = k, delta = c, log = TRUE))
    
  } else if(family=="Burr"){
    k<- exp(log.hyper.mu[1])
    c<- exp(log.hyper.mu[2])
    #sigma<- exp(mu)/(((2^(1/k))-1)^(1/c))
    sigma<- exp(mu) / qburr(0.5, shape1 = k, shape2 = c, scale = 1)
    log.post_A<- sum(dburr(x=A, shape1=k, shape2 =c, scale = sigma, log=TRUE)) 
  } else if(family=="generalized-Gamma"){
    k<- exp(log.hyper.mu[1])
    c<- exp(log.hyper.mu[2])
    #sigma<- exp(mu) * gamma(k/c) / gamma((k+1)/c)
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    log.post_A<- sum(dggamma(t=A, theta = sigma, kappa = k, delta = c, log = TRUE)) 
  }
  return(log.post_A+log.post.process_mu)
}

########################################################################################################
################ Gradients of the log-likelihood wrt mu
########################################################################################################
grad_mu<-function(A, mu, intercept2, W1, W2, beta, beta2,  kappa_mu, log.hyper.mu, Z2, A2, s.fixed,  family){
  n2<-length(A)
  grad.at.process_mu<- -kappa_mu*mu  + kappa_mu * (intercept2+Z2%*%beta2+beta*(A2%*%W1)+ A2%*%W2)
  if(family=="Gamma-gamma"){
    gamma1<- exp(log.hyper.mu[1])
    gamma2<- exp(log.hyper.mu[2])
    alpha<- exp(mu)/qf(0.5, df1=gamma1, df2=gamma2)
    grad_A<- -((0.5*gamma1)) + (0.5 * (gamma1+gamma2) * gamma1 * A / (gamma1 * A + gamma2 * alpha)) 
  } else if(family=="extended-GPD"){
    k<- exp(log.hyper.mu[1])
    xi<- exp(log.hyper.mu[2])
   # xi<- 2 * (exp(log.hyper.mu[2])) / (exp(log.hyper.mu[2])+1)
    # if(k < 0.2){
    #   k<-0.2
    # } else{
    #   k<-k
    # }
    # if(abs(xi)> 0.001){
    #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
    # } else{
    #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
    # }
    sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    #grad_A<- -1 + (1+xi) * A /(sigma + xi * A) -  (((k-1) * A * (sigma + xi * A)^((-1/xi)-1)) / ((sigma^(-1/xi)) - (sigma + xi * A)^(-1/xi)))
    grad_A<- ((A-sigma) /(sigma + xi * A)) -  ((k-1) * A * evd::dgpd(x=A/sigma, loc=0, scale = 1, shape = xi) / (sigma* evd::pgpd(q=A/sigma, loc=0, scale = 1, shape = xi)))
    if(mean(is.na(grad_A))!=0){ ## to avoid numerical instability
      ind.na<-is.na(grad_A)
      grad_A[ind.na]<- 10^(-5)
    }
    } else if (family=="GPD"){
      k<-1
      xi<- exp(log.hyper.mu)
      #xi<- 2 * (exp(log.hyper.mu)) / (exp(log.hyper.mu)+1)
      # if(abs(xi)> 0.001){
      #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
      # } else{
      #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
      # }
      #grad_A<- -1 + ((1+xi)* A /(sigma + xi * A)) 
      sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
     # grad_A<- -1 + (1+xi) * A /(sigma + xi * A) -  (((k-1) * A * (sigma + xi * A)^((-1/xi)-1)) / ((sigma^(-1/xi)) - (sigma + xi * A)^(-1/xi)))
      grad_A<- ((A-sigma) /(sigma + xi * A)) -  ((k-1) * A * evd::dgpd(x=A/sigma, loc=0, scale = 1, shape = xi) / (sigma* evd::pgpd(q=A/sigma, loc=0, scale = 1, shape = xi)))
      
  } else if (family=="Gamma"){
    # kappa_a<- exp(log.hyper.mu)
    # grad_A<- - s.fixed * kappa_a * (1 - A * exp(-mu)) 
    k<- exp(log.hyper.mu)
    c<- 1
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    #sigma<- exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    grad_A<- -k + c * ((A/sigma)^c) 
  } else if (family=="log-Gamma"){
    kappa_a<- exp(log.hyper.mu)
    #rate.a<- 1/(1-((1+exp(mu))^(-1/kappa_a)))
    # if(kappa_a>1){
    #   rate.a<- ((kappa_a-1)/(log(1+exp(mu))))-1  #Covariates in mode
    # } else {
    #   rate.a<- 10^10
    # }
    # rate.a<- qgamma(0.5, shape = kappa_a, rate = 1)/log(1+exp(mu)) ## in median
    # grad_A<- -((kappa_a/rate.a)-(log(A+1))) * ((rate.a^2)/kappa_a) * exp(mu) *(1+exp(mu))^(-1-(1/kappa_a)) 
    # 
    sigma<- qgamma(0.5, shape = kappa_a, rate = 1)/log(1+exp(mu)) ## in median
    rate.a<- 1/sigma
    grad_A<- -((kappa_a/rate.a)-(log(A+1))) * (1/(1+exp(-mu))) * qgamma(0.5, shape = kappa_a, rate = 1) / ((log(1+exp(mu)))^2)

  } else if (family=="Weibull"){
    # kappa_a<- exp(log.hyper.mu)
    # lambda<- exp(mu)/(log(2)^(1/kappa_a))  # in median
    # #lambda<- exp(mu)/gamma(1+kappa_a)
    # grad_A<- - kappa_a * (1 - (A/lambda)^(kappa_a)) 
    
    k<- exp(log.hyper.mu)
    c<- k
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    #sigma<- exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    grad_A<- - k + c * ((A/sigma)^c) 

  } 
  else if (family=="Burr"){
    k<- exp(log.hyper.mu[1])
    c<- exp(log.hyper.mu[2])
   # sigma<- exp(mu)/(((2^(1/k))-1)^(1/c))
    sigma<- exp(mu) / qburr(0.5, shape1 = k, shape2 = c, scale = 1)
    grad_A<-  -c + c*(k+1)/(1+(A/sigma)^(-c)) 
  } else if(family=="generalized-Gamma"){
    k<- exp(log.hyper.mu[1])
    c<- exp(log.hyper.mu[2])
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    #sigma<- exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    grad_A<- -k + c * ((A/sigma)^c) 
  }
  return(grad_A+grad.at.process_mu)
}




