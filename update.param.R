#####################################################################################
################################## function to simulate the kappa_eta ##########################
######################################################################################
kappa_eta_sim<-function(eta, intercept1, beta1, W1, Z1, A1, model, fixed.kappa_eta){
  n1<-length(eta)
  if(model=="Model1" | model=="Model3"){
  sim_kappa_eta<-rgamma(1, shape=hyper_fixed[1]+0.5*n1, 
                        rate= hyper_fixed[2] + 0.5*(sum((eta-intercept1-Z1%*%beta1-A1%*%W1)^2)))
  } else if(model=="Model2" | model=="Model4"){
    sim_kappa_eta<-fixed.kappa_eta
  }
  return(sim_kappa_eta)
}

######################################################################################
################################## function to simulate the kappa_w1 ##############
######################################################################################
kappa_w1_sim<-function(W1, node1, node2){
  N<-length(W1)
  sim_kappa_w1<-rgamma(1, shape=hyper_fixed[3]+0.5*(N-1), 
                       rate= hyper_fixed[4] + 0.5*sum((W1[node1]-W1[node2])^2))
  return(sim_kappa_w1)
}


######################################################################################
################################## function to simulate the kappa_w2 ##############
######################################################################################
kappa_w2_sim<-function(W2, node1, node2){
  N<-length(W2)
  sim_kappa_w2<-rgamma(1, shape=hyper_fixed[5]+0.5*(N-1), 
                       rate= hyper_fixed[6] + 0.5*sum((W2[node1]-W2[node2])^2))
  return(sim_kappa_w2)
}


######################################################################################
################################## function to simulate the kappa_mu ##########################
######################################################################################
kappa_mu_sim<-function(mu, intercept2, beta2, beta, W1, W2, Z2, A2, model, fixed.kappa_mu){
  n2<-length(mu)
  if(model=="Model1" | model=="Model4"){
  sim_kappa_mu<-rgamma(1, shape=hyper_fixed[7] + 0.5 * n2, 
                       rate= hyper_fixed[8] + 0.5 * (sum((mu-intercept2-Z2%*%beta2-beta*(A2%*%W1)-A2%*%W2)^2)))
  } else if(model=="Model2" | model=="Model3"){
    sim_kappa_mu<-fixed.kappa_mu
  }
  return(sim_kappa_mu)
}


######################################################################################
################################## function to simulate the log.hyper.mu ##########################
######################################################################################
prop_log.hyper.mu<- function(A, mu, log.hyper.mu, sigma2.hyper.mu, hyper.mu_fixed, conjugate.hyper.mu, family){
if(conjugate.hyper.mu=="FALSE"){
  mean<-log.hyper.mu
  var<-rep(sigma2.hyper.mu, length(log.hyper.mu))
  prop<-rnorm(n=length(log.hyper.mu), mean=mean, sd=sqrt(var))
} else if (conjugate.hyper.mu=="TRUE" & family=="log-normal"){
  n2<-length(A)
 # prop<-log(rgamma(1, shape=hyper.mu_fixed[1]+0.5*n2, rate= hyper.mu_fixed[2] + 0.5*(sum((log(A)-mu)^2))))  ## in log scale
  prop<- rgamma(1, shape=hyper.mu_fixed[1]+0.5*n2, rate= hyper.mu_fixed[2] + 0.5*(sum((log(A)-mu)^2)))  ## in log scale
  
 }
return(prop)
}


######################################################################################
################################## function to simulate the beta ##############
######################################################################################
t.A2.A2<-t(A2)%*%A2
beta_sim<-function(mu, intercept2, kappa_mu, beta2, W1, W2, Z2, A2, t.A2.A2){
  prec_beta<-kappa_mu * (t(W1) %*% t.A2.A2 %*% W1) + hyper_fixed[9]
  mean_beta<-kappa_mu*(t(W1)%*% t(A2)%*% (mu-intercept2-Z2%*%beta2- A2%*%W2))/prec_beta  ## may be we can simplify it later
  sim_beta<-rnorm(n=1, mean=mean_beta, sd=sqrt(1/prec_beta))
  return(sim_beta)
}

######################################################################################
####################### function to simulate the intercept1  ##############
######################################################################################
intercept1_sim<-function(eta, beta1, W1, kappa_eta, Z1, A1){
  n1<-length(eta)
  prec.intercept1<- hyper_fixed[10]+n1*kappa_eta
  mean.intercept1<- kappa_eta*sum((eta-Z1%*%beta1-A1%*%W1))/prec.intercept1
  sim_intercept1<-rnorm(n=1, mean=mean.intercept1, sd=sqrt(1/prec.intercept1))
  return(sim_intercept1)
}

######################################################################################
####################### function to simulate the intercept2  ##############
######################################################################################
intercept2_sim<-function(mu, beta2, beta, W1, W2, kappa_mu, Z2, A2){
  n2<-length(mu)
  prec.intercept2<-hyper_fixed[11]+n2*kappa_mu
  mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2-beta*(A2%*%W1)-A2%*%W2)/prec.intercept2
  sim_intercept2<-rnorm(n=1, mean=mean.intercept2, sd=sqrt(1/prec.intercept2))
  return(sim_intercept2)
}

############################################################################
################################## function to simulate beta1 ##############
############################################################################
beta1_sim<-function(eta, intercept1, W1, kappa_eta, Z1, A1){
  p<-ncol(Z1)
  latent.cov.inv<- hyper_fixed[12] * diag(1,p) + kappa_eta * (t(Z1)%*%Z1)
  latent.mean.part<- kappa_eta * (t(Z1)%*%(eta-intercept1-A1%*%W1))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(p)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
}

############################################################################
################################## function to simulate beta2 ##############
############################################################################
beta2_sim<-function(mu, intercept2, beta, kappa_mu, W1, W2, Z2, A2){
  q<-ncol(Z2)
  latent.cov.inv<- hyper_fixed[13] * diag(1,q) + kappa_mu * (t(Z2)%*%Z2)
  latent.mean.part<- kappa_mu * (t(Z2)%*%(mu-intercept2-beta*(A2%*%W1)-A2%*%W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
}
######################################################################################
################################## Simulating W1 ##############
######################################################################################
t.A1.A1<-t(A1)%*%A1
t.one.one<-rep(1, dim(Q)[1])%*% t(rep(1, dim(Q)[1]))
W1_sim<-function(eta, mu, intercept1, intercept2, beta1, beta2, W2, kappa_w1, kappa_eta, kappa_mu, beta, Q, Z1, Z2, A1, A2, t.A1.A1, t.A2.A2, t.one.one, tau.w1_sum.const){
  N<-dim(Q)[1]
  # latent.cov.inv<- kappa_w1 * Q +  kappa_eta * t.A1.A1 + (beta^2) * kappa_mu * t.A2.A2 +  
  #   (1/(N*tau.w1_sum.const))* t.one.one  #### adding prior in  W1~Normal(0, tau.w1_sum.const*N) to make it close to zero
  latent.cov.inv<- kappa_w1 * Q +  kappa_eta * t.A1.A1 + (beta^2) * kappa_mu * t.A2.A2 ### Without any prior but at every iterations we need to replave W1 by W1-mean(W1)
  latent.mean.part<- beta * kappa_mu * (t(A2)%*%(mu-intercept2-Z2%*%beta2-A2%*%W2)) + kappa_eta * (t(A1)%*%(eta-intercept1-Z1%*%beta1))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(N)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
}

######################################################################################
################################## Simulating W2 ##############
######################################################################################
W2_sim<-function(mu, intercept2, beta2, W1, kappa_w2, kappa_mu, beta, Q, Z2, A2, t.A2.A2, t.one.one, tau.w2_sum.const){
  N<-dim(Q)[1]
  # latent.cov.inv<- kappa_w2 * Q +  kappa_mu * t.A2.A2 + 
  #   (1/(N*tau.w2_sum.const)) * t.one.one #### adding prior in W2~Normal(0, tau.w2_sum.const*N) to make it close to zero
  latent.cov.inv<- kappa_w2 * Q +  kappa_mu * t.A2.A2 ### Without any prior but at every iterations we need to replave W1 by W2-mean(W2)
  latent.mean.part<- kappa_mu * (t(A2)%*%(mu-intercept2-Z2%*%beta2-beta*(A2%*%W1)))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(N)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
}
######################################################################################
########################### Simulating mu using MALA and reverse proposal##############
######################################################################################

prop_mu<-function(A, mu, intercept2, W1, W2, beta, beta2,  kappa_mu, log.hyper.mu, Z2, A2, sigma2.mu, conjugate.mu,  s.fixed, family){
  n2<-length(A)
  if(conjugate.mu=="FALSE"){
  mean<-mu + 0.5*sigma2.mu*grad_mu(A, mu, intercept2, W1, W2, beta, beta2,  kappa_mu, log.hyper.mu, Z2, A2,  s.fixed, family)
  var<-rep(sigma2.mu, n2)
  proposals<-rnorm(n=n2, mean=mean, sd=sqrt(var))
  } else if (conjugate.mu==TRUE & family=="log-normal"){
    #prec.mu<-(exp(log.hyper.mu) + kappa_mu)
    #mean.mu<- (kappa_mu * (intercept2+Z2%*%beta2+beta*(A2%*%W1)+A2%*%W2) + exp(log.hyper.mu)*log(A))/prec.mu
    prec.mu<-(log.hyper.mu + kappa_mu)
    mean.mu<- (kappa_mu * (intercept2+Z2%*%beta2+beta*(A2%*%W1)+A2%*%W2) + log.hyper.mu*log(A))/prec.mu
    proposals <- rnorm(n=n2, mean =mean.mu, sd=sqrt(1/prec.mu))
  }
  return(proposals)
}

# ## MALA reverse density
MALA_dens_mu<-function(A, mu, mu.star, intercept2, W1, W2, beta, beta2,  kappa_mu, log.hyper.mu, Z2, A2, sigma2.mu,  s.fixed, family){
  n2<-length(A)
  mean<-mu + 0.5*sigma2.mu*grad_mu(A, mu, intercept2, W1, W2, beta, beta2,  kappa_mu, log.hyper.mu, Z2, A2,  s.fixed, family)
  var<-rep(sigma2.mu, n2)
  dens<-sum(dnorm(x=mu.star, mean=mean, sd=sqrt(var), log=TRUE))
  return(dens)
}


################################################################
################################## Simulating eta ##############
#################################################################
log.post_eta<-function(Y, eta, intercept1, beta1, kappa_eta, W1, Z1, A1){
  log.post<-sum(Y*eta)-sum(exp(eta))-0.5 * kappa_eta * sum(eta^2) +
    kappa_eta * t(eta) %*% (intercept1+Z1%*%beta1+ A1%*%W1)
  return(log.post)
}

grad_eta<-function(Y, eta, intercept1, beta1, kappa_eta, W1, Z1, A1){
  grad<- Y-exp(eta)-kappa_eta * eta + kappa_eta * (intercept1+Z1%*%beta1+ A1%*%W1)
  return(grad)
}

## MALA proposals
MALA_prop_eta<-function(Y, eta, intercept1, beta1, kappa_eta, W1, Z1, A1, sigma2_eta){
  n<-length(Y)
  mean<-eta+0.5*sigma2_eta*grad_eta(Y=Y, eta=eta, intercept1=intercept1, beta1=beta1, kappa_eta=kappa_eta, W1=W1, Z1=Z1, A1=A1)
  var<-rep(sigma2_eta, n)
  prop<-rnorm(n=n, mean=mean, sd=sqrt(var))
  return(prop)
}
# ## MALA reverse density
MALA_dens_eta<-function(Y, eta, eta.star, intercept1, beta1, kappa_eta, W1, Z1, A1, sigma2_eta){
  n<-length(Y)
  mean<-eta+0.5*sigma2_eta*grad_eta(Y=Y, eta=eta, intercept1=intercept1, beta1=beta1, kappa_eta=kappa_eta, W1=W1, Z1=Z1, A1=A1)
  var<-rep(sigma2_eta, n)
  dens<-sum(dnorm(x=eta.star, mean=mean, sd=sqrt(var), log=TRUE))
  return(dens)
}



