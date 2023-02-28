#rm(list=ls())
###############################################################################
######### Marked point processe #############
###############################################################################
#' @Model:
# Y given \eta independent Poisson( exp(eta_i)),
# eta=\beta_1 Z_1+W_1, 
# A given \mu independent Log-normal (\mu_i, 1/kappa_a)
# mu=\beta_2 Z_2+ \beta W_1+ W_2, where W_2 independent of W_1 

#' @Prior:
# W_1:- ICAR priors: multivariate Gaussian of dimension (n-1) with mean zero and precision matrix \kappa^{w_1} Q 
# W_2:- ICAR priors: multivariate Gaussian of dimension (n-1) with mean zero and precision matrix \kappa^{w_2} Q, such that W_1 \indp W_2
# \beta_1: multivarite Gaussian of dimension (p) with mean zero and precision \kappa_{\beta_1}=0.001
# \beta_2: multivarite Gaussian of dimension q with mean zero and precision \kappa_{\beta_2}=0.001
# \beta: Gaussian of dimension 1 with mean zero and precision \kappa_{\beta}=0.001
# \kappa_{w1}: Gamma with shape 0.5 and rate 0.5
# \kappa_{w2}: Gamma with shape 0.5 and rate 0.5
# \kappa_{a}: Gamma with shape 0.5 and rate 0.5

#' @Data:
# Y: data of dimension nx1
# A: the marked data of dimension
# Z_1: Set of first covariates of dimesion pxn
# Z_2: Set of first covariates of dimesion qxn

####################################################################################################
## Extracting the adjacency structure that is the same as the slope unit of the landslides data ####
####################################################################################################
############################################################
####################### Defining the projection matrices
############################################################

fun_proj_matrix<-function(Q){
  N<-dim(Q)[1]
  set.seed(1)
  #sample.ind<-sample(1:3, size = N, replace = TRUE)
  sample.ind<-sample(1, size = N, replace = TRUE)
  n.count<-sum(sample.ind)
  A1<-matrix(0, nrow = n.count, ncol = N)
  jj<-1
  for (i in 1:N) {
    initial.p<-jj
    last.p<-jj+sample.ind[i]-1
    A1[initial.p:last.p,i]<-1
    jj<-jj+sample.ind[i]
  }
  ## Projection matrix for Size 
  set.seed(2)
  #sample.ind<-sample(1:3, size = N, replace = TRUE)
  sample.ind<-sample(1, size = N, replace = TRUE)
  
  n.size<-sum(sample.ind)
  n.size
  A2<-matrix(0, nrow = n.size, ncol = N)
  jj<-1
  for (i in 1:N) {
    initial.p<-jj
    last.p<-jj+sample.ind[i]-1
    A2[initial.p:last.p,i]<-1
    jj<-jj+sample.ind[i]
  }
  return(list("A1"=A1,"A2"=A2, "n1"=n.count, "n2"=n.size))
  
}

############################################################
######### Simulation of Marked point processe  #############
############################################################
# sim_mark_function: function to simulate the Marked point procees defined above
## Input
#' @Q: the precision matrix in CAR priors
#' @hyper.prec = c(taw_w1, taw_w2, taw_a): Hyperparameters
#' @covar.param =c(beta1, beta2): parameters corresponding to the two covariates
#' @Z1, @Z2: Coavariates
## Ouput
#' @Y: Data 
sim_mark_function<-function(Q,  other.hyper, beta1, beta2, Z1, Z2, A1, A2, family){
  N<-dim(Q)[1]
  n1<-nrow(A1)
  n2<-nrow(A2)
  p<-ncol(Z1)
  q<-ncol(Z2)
  log.hyper.mu<-fun.hyper.mu(family = family)
  kappa_w1<-other.hyper[1]
  kappa_w2<-other.hyper[2]
  kappa_eta<-other.hyper[3]
  kappa_mu<- other.hyper[4]
  intercept1<-other.hyper[5]
  intercept2<-other.hyper[6]
  beta<-other.hyper[7]
  beta1<-beta1
  beta2<-beta2
  W1_s<-rMVNormP_eigen(n=1, mu=rep(0,N), Sigma = kappa_w1*Q)
  W2_s<-rMVNormP_eigen(n=1, mu=rep(0,N), Sigma = kappa_w2*Q)
  eta<-intercept1+Z1%*%beta1+ (A1 %*% W1_s) + rnorm(n1, mean=0, sd=sqrt(1/kappa_eta))
  Y=rpois(n=n1, lambda = exp(eta))
  mu<-intercept2+Z2%*% beta2+ beta* (A2%*%W1_s) + (A2 %*% W2_s)+ rnorm(n2, mean=0, sd=sqrt(1/kappa_mu))
  ### simulating it for different mark distributions 
  if(family=="Gamma-gamma"){
    gamma1<- exp(log.hyper.mu[1])
    gamma2<- exp(log.hyper.mu[2])
    alpha<- exp(mu)/qf(0.5, df1=gamma1, df2=gamma2)
    A<- alpha*rf(n=n2, df1 = gamma1, df2=gamma2)
    } else if(family=="extended-GPD"){
      k<- exp(log.hyper.mu[1])
      xi<-  exp(log.hyper.mu[2])
      #xi<- 2 * (exp(log.hyper.mu[2])) / (exp(log.hyper.mu[2])+1)
      # if(xi> 0.001){
      #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
      # } else{
      #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
      # }
      sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
      A<- rEGPD1(n=n2, k=k, xi=xi, sigma=sigma)
  } else if(family=="GPD"){
    k<-1
    xi<- exp(log.hyper.mu)
    #xi<- 2 * (exp(log.hyper.mu)) / (exp(log.hyper.mu)+1)
    # if(xi> 0.001){
    #   sigma<- exp(mu) * xi / (((1-((0.5)^(1/k)))^(-xi))-1)
    # } else{
    #   sigma<- -exp(mu)/(log(1-((0.5)^(1/k))))
    # }
    sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    A<-rEGPD1(n=n2, k=k, xi=xi, sigma=sigma)
  } else if(family=="generalized-Gamma"){
    k<-exp(log.hyper.mu[1])
    c<-exp(log.hyper.mu[2])
    #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    A<- rggamma(n=n2, theta = sigma,  kappa=k, delta=c)
  } else if(family=="Gamma"){
    # kappa_a<-exp(log.hyper.mu)
    # A<-rgamma(n=n2, shape = s_fixed* kappa_a, rate = s_fixed * kappa_a* exp(-mu))
    k<-exp(log.hyper.mu)
    c<-1
    #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    A<- rggamma(n=n2, theta = sigma,  kappa=k, delta=c)
      } else if(family=="log-Gamma"){
        kappa_a<-exp(log.hyper.mu)
        #rate.a<-1/(1-((1+exp(mu))^(-1/kappa_a)))
        # if(kappa_a>1){
        #   rate.a<- ((kappa_a-1)/(log(1+exp(mu))))-1  #Covariates in mode
        # } else {
        #   rate.a<- 10^10
        # }
        sigma<- qgamma(0.5, shape = kappa_a, rate = 1)/log(1+exp(mu)) ## in median
        rate.a<- 1/sigma
       A<- exp(rgamma(n=n2, shape = kappa_a, rate = rate.a))-1
  } else if(family=="Weibull"){
    # kappa_a<-exp(log.hyper.mu)
    # #lambda<-exp(mu)/gamma(1+(1/kappa_a))
    # lambda<- exp(mu)/(log(2)^(1/kappa_a))  # in median
    # A<-rweibull(n=n2, shape = kappa_a, scale = lambda)
    
    k<-exp(log.hyper.mu)
    c<-k
    #sigma<-exp(mu) * gamma(k/c) / gamma((k+1)/c) ## mean
    sigma<- exp(mu) / ((qgamma(0.5, shape = k/c, rate = 1))^(1/c)) ## median
    A<- rggamma(n=n2, theta = sigma,  kappa=k, delta=c)
    
  } else if (family=="Burr"){
    k<- exp(log.hyper.mu[1])
    c<- exp(log.hyper.mu[2])
    #sigma<- exp(mu)/(((2^(1/k))-1)^(1/c))
    sigma<- exp(mu) / qburr(0.5, shape1 = k, shape2 = c, scale = 1)
    A<- rburr(n=n2, shape1 = k, shape2 = c, scale = sigma)
  } else if(family=="log-normal"){
   # kappa_a<- exp(log.hyper.mu)
    kappa_a<- log.hyper.mu ### original scale
    
    A<-exp(rnorm(n=n2, mean =mu, sd=sqrt(1/kappa_a)))
  }
  return(list(Y=Y, mu=mu, A=A, W1=W1_s, eta=eta, W2=W2_s, A1=A1, A2=A2))
}

####det(Q), quite misleading because of the entries inside the matrix, so the numerical instabilities
## simulating from the multivarite singular normal





