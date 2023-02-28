#############################################
########## Main MCMC function  ##############
#############################################
MCMC.main_function<-function(N.MCMC, Y, ind_NA_Y, A, ind_NA_A, Q, Z1, Z2, by, adapt, burning1, burning2, s.fixed,
                    sigma2_eta, theta_eta, sigma2.mu, theta_mu, sigma2.hyper.mu, theta_hyper.mu, 
                    A1, A2, family, CV1_number, k_fold_CV, print.result=TRUE, traceplot=FALSE, conjugate.hyper.mu=FALSE, conjugate.mu=FALSE, 
                    model, model.base, fixed.kappa_mu, fixed.kappa_eta, CV, true.values=NULL, max.A, simulation)
{  
  n1<-length(Y)
  n2<- length(A)
  n.Q<-dim(Q)[1]
  p<-ncol(Z1)
  q<-ncol(Z2)
  
  ####### Imputation in case of within sample diagnostics
  imputed.Y.WSD<-rep(0, n1)  # Imputed post. mean of Y in case of within sample diganostics=Impute.Y.WSD / N-(burning1+burning2)
  imputed.A.WSD<-rep(0, n2) # Imputed post. mean of Y in case of within sample diganostics=Impute.A.WSD/ N-(burning1+burning2)
  imputed.Y.WSD.samples<-rep(0, n1) # Used for std of Y in case of within sample diganostics
  imputed.A.WSD.samples<-rep(0, n2) # Used for std of A in case of within sample diganostics
  
  ####### Imputation in case of out of sample diagnostics
  imputed.Y.OSD<-rep(0, times=sum(ind_NA_Y))  # Imputed  summation values for Y
  imputed.A.OSD<-rep(0, times=sum(ind_NA_A))  # Imputed  summation values for A
  imputed.Y.OSD.samples<-c() # Imputed  summation^2 values for Y
  imputed.A.OSD.samples<-c() # Imputed  summation^2 values for  A

  ####### posterior summary of latent parameters
  post.sum.mu<-rep(0, times=n2)  ## Post mean for mu=post.sum.mu/N-(burning1+burning2)
  post.sum.eta<-rep(0, times=n1) ## Post mean for eta=post.sum.eta/N-(burning1+burning2)
  post.sum.mu.samples<-c() ## Post variance for mu=post.sum2.mu/N-(burning1+burning2)-(Post mean for mu)^2
  post.sum.eta.samples<-c() ## Post variance for eta=post.sum2.eta/N-(burning1+burning2)-(Post mean for eta)^2
  
  
  ## Storing the tuning parameteres
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt),ncol=3,byrow=TRUE) ## storing the adaptive scale parameters
  sigma.matrix[1,]<-c(sigma2_eta, sigma2.mu, sigma2.hyper.mu)
  
  init.hyper.mu_and_name<-init_fun_log.hyper.mu(family=family)
  cur.samples.log.hyper.mu<-init.hyper.mu_and_name$log.hyper.mu.init #current samples for gamma1
  init.all.other.param_and.name<-init_fun_all_other_param(Z1=Z1, Z2=Z2, A=A, Y=Y, Q=Q)
  init.other.param<-init.all.other.param_and.name$init.all.other.param
  cur.samples.kappa_w1<-init.other.param[1] #current samples for \kappa_w1
  cur.samples.kappa_w2<-init.other.param[2] #current samples for \kappa_w2
  cur.samples.kappa_eta<-init.other.param[3] #current samples for \kappa_w2
  cur.samples.kappa_mu<-init.other.param[4] #current samples for \kappa_w2
  cur.samples.intercept1<-init.other.param[5]
  cur.samples.intercept2<-init.other.param[6]
  cur.samples.beta<-init.other.param[7]
  cur.samples.beta1<-init.other.param[8:(8-1+p)]    #current samples for beta1
  cur.samples.beta2<-init.other.param[(8+p):(8-1+p+q)]    #current samples for beta2
  cur.samples.eta<-init.other.param[(8+p+q):((8-1+p+q+n1))]    #current samples for eta
  cur.samples.w1<-init.other.param[(8+p+q+n1):((8-1+p+q+n1+n.Q))]    #current samples for w1
  cur.samples.mu<-init.other.param[((8+p+q+n1+n.Q)):((8-1+p+q+n1+n.Q+n2))] #current samples for mu
  cur.samples.w2<-init.other.param[((8+p+q+n1+n.Q+n2)):((8-1+p+q+n1+n.Q+n2+n.Q))]  #current samples for w2
  # cur.samples.w1<- cur.samples.w1-mean(cur.samples.w1) ## making sure that sum of W1 is zero
  # cur.samples.w2<- cur.samples.w2-mean(cur.samples.w2) ## making sure that sum of W2 is zero
  ### saving the Chain for all the hyperparameters and some latent parmeters
  samples <- matrix(nrow=floor(N.MCMC/by),ncol=length(cur.samples.log.hyper.mu)+7+p+q+5,byrow=TRUE) # storing the samples
  samples[1,]<- c(cur.samples.log.hyper.mu, init.other.param[c(1:7, 8:(7+p), (8+p):(7+p+q), (8+p+q):(8+p+q+1), (8+p+q)+n1, (8+p+q+n1+n.Q), (8+p+q+n1+n.Q+n2))]) ## saving only few samples
  
  
  
  j<-1
  l<-1
  m<-1
  k<-1
  for (i in 1:(N.MCMC-1)){
    # if((i%%(adapt))-1==0 & i< (burning1+burning2+2)){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
    #   rate.latent<-0
    #   rate.hyper<-0
    # }
    
    if(((i%%(adapt))-1==0) & (i< (burning1+burning2+2))){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
      rate.eta<- 0
      rate.mu<- 0
      rate.hyper.mu<- 0
    }
   
    if((i%%by)-1==0){  ## by is the thinning steps and adapt is the number of iterations after which i update the variance of MALA and random walk algorithms
      par(mfrow=c(6,6),oma=c(0,0,2,0),mar=c(4,5,1,1))
      if(print.result==TRUE){
      if(i< (burning1+burning2+2)){
        print(paste("Iteration: ",i,
                    "; Accep rate eta=",rate.eta/((i%%(adapt))+1),
                    "; sigma eta=",sigma2_eta,
                    "; Accep rate mu=",rate.mu/((i%%(adapt))+1),
                    "; sigma mu=",sigma2.mu,
                    "; Accep rate hyper.mu=",rate.hyper.mu/((i%%(adapt))+1),
                    "; sigma hyper.mu=",sigma2.hyper.mu,sep=""))
      } else{
        print(paste("Iteration: ",i,
                    "; Accep rate eta",rate.eta/(i-(burning1+burning2+2)),
                    "; sigma eta=",sigma2_eta,
                    "; Accep rate mu",rate.mu/(i-(burning1+burning2+2)),
                    "; sigma mu=",sigma2.mu,
                    "; Accep rate hyper.mu",rate.hyper.mu/(i-(burning1+burning2+2)),
                    "; sigma hyper.mu=",sigma2.hyper.mu,sep=""))
      }
      }
      
      if(traceplot==TRUE){
        model.param.name<-c(init.hyper.mu_and_name$log.hyper.mu.init,init.all.other.param_and.name$model.param.name.all.other.param)
        ## Plotting the traceplots
        for (ll  in 1: length(model.param.name)) {
          if(simulation==TRUE){
            plot(by*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            abline(h=true.values[ll], col=2)
          }  else {
            plot(by*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            
          }
        }
      }
      
      if((i%%adapt)-1==0){
        #plot(sigma.matrix[1:k,1],sigma.matrix[1:k,2],xlab="sigma.latent",ylab="sigma.hyper2") #plot for the scale parameters chosen adaptively
        k<-k+1
      }
      
      l<-l+1
    }
    
    # print(rate.hyper.mu)
    # print(i)
    # print(sigma2.hyper.mu)
    
    cur.eta<- cur.samples.eta
    cur.log.hyper.mu<-cur.samples.log.hyper.mu
    cur.mu<-cur.samples.mu
    
    # print(range(cur.mu))
    # print(cur.log.hyper.mu)
    # print(range(cur.eta))
    # 
    # print(sample.mu)
    # print(sample.hyper.mu)
    
    ### Imputation of Y
    # imputed_Y<-impute.NA.Y(ind_NA_Y = ind_NA_Y, eta = cur.eta, CV=CV)
    #   Y[ind_NA_Y]<- imputed_Y   #replacing the NA by the imputed values
    #   ### Imputation of A
    #   imputed_A<-impute.NA.A(ind_NA_A = ind_NA_A, mu=cur.mu, log.hyper.mu = cur.log.hyper.mu, family = family, CV=CV)
    #   A[ind_NA_A]<-imputed_A
    #   
      
      if(CV=="OSD"){
        ### Imputation of Y
        imputed_Y<- impute.NA.Y(ind_NA_Y = ind_NA_Y, eta = cur.eta, CV=CV)
        Y[ind_NA_Y]<- imputed_Y   #replacing the NA by the imputed values

        ### Imputation of A
        imputed_A<-impute.NA.A(ind_NA_A = ind_NA_A, mu=cur.mu, log.hyper.mu = cur.log.hyper.mu, family = family, CV=CV)
        A[ind_NA_A]<-imputed_A

        if(i>burning1+burning2){ ## storing the samples: Posterior predictive mean and standard variability
          imputed.Y.OSD<-imputed.Y.OSD+imputed_Y
          imputed.A.OSD<-imputed.A.OSD+imputed_A
          if((i%% floor((N.MCMC-(burning1+burning2))/250))==0){
            imputed.Y.OSD.samples<- rbind(imputed.Y.OSD.samples, imputed_Y)
          imputed.A.OSD.samples<- rbind(imputed.A.OSD.samples, imputed_A)
          }
        }
      }
         if(CV=="WSD"){
        if(i>burning1+burning2){
          imputed_Y<-impute.NA.Y(ind_NA_Y = ind_NA_Y, eta = cur.eta, CV=CV)
          imputed_A<-impute.NA.A(ind_NA_A = ind_NA_A, mu=cur.mu, log.hyper.mu = cur.log.hyper.mu, family = family, CV=CV)
          imputed.Y.WSD<-imputed.Y.WSD+imputed_Y
          imputed.A.WSD<-imputed.A.WSD+imputed_A
        
        if((i%% floor((N.MCMC-(burning1+burning2))/250))==0){
        imputed.Y.WSD.samples<- rbind(imputed.Y.WSD.samples, imputed_Y)
        imputed.A.WSD.samples<- rbind(imputed.A.WSD.samples, imputed_A)
        }
        }
         }
    
    if(i>burning1+burning2){
      post.sum.mu <- post.sum.mu + cur.mu
      post.sum.eta<- post.sum.eta + cur.eta
    }
      
    
     # for numerical instabilities of A for the cases of extended GPD and GPD cases
  #     if(family=="extended-GPD" | family=="GPD"){
  #     if((min(A)<1)| (max(A)>max.A+10^3)){
  #       ind.A.small<- which(A<1)
  #       ind.A.big<- which(A> max.A+10^3)
  #       A[ind.A.small]<- 1 ## this is only for to avoid problems with extended-GPD, otherwise we dont need it
  #       A[ind.A.big]<- max.A+10^3 ## this is only for to avoid problems with extended-GPD, otherwise we dont need it
  #     }
  # }
    #### Proposing new parameters for kappa_eta (Gibbs steps) #####
    cur.samples.kappa_eta<-kappa_eta_sim(eta = cur.eta, intercept1 = cur.samples.intercept1, 
                                         beta1 = cur.samples.beta1, W1=cur.samples.w1, Z1=Z1, 
                                         A1=A1,  model=model, fixed.kappa_eta=fixed.kappa_eta)
    #### Proposing new parameters for kappa_w1 (Gibbs steps) #####
    cur.samples.kappa_w1<- kappa_w1_sim(W1=cur.samples.w1, node1 = node.set[,1], node2 = node.set[,2])  
    #cur.samples.kappa_w1<- kappa_w1
    ###### Proposing new parameters for kappa_w2 (Gibbs steps) #####
    cur.samples.kappa_w2<- kappa_w2_sim(W2=cur.samples.w2, node1 = node.set[,1], node2 = node.set[,2])  # Proposed hyperparameters using uniform random walks
    
    #cur.samples.kappa_a<- kappa_a  # Proposed hyperparameters using uniform random walks
    #### Proposing new parameters for kappa_a (Gibbs steps) #####
    cur.samples.kappa_mu<- kappa_mu_sim(mu=cur.mu, intercept2 = cur.samples.intercept2, beta2 = cur.samples.beta2,
                                        beta=cur.samples.beta, W1=cur.samples.w1, W2=cur.samples.w2, Z2=Z2, 
                                        A2=A2, model=model, fixed.kappa_mu=fixed.kappa_mu) # Proposed hyperparameters using uniform random walks
    #cur.samples.kappa_a<- kappa_a  # Proposed hyperparameters using uniform random walks
    
    #### Proposing new parameters for hyper.mu MH steps #####
    prop.log.hyper.mu<-prop_log.hyper.mu(A=A, mu=cur.mu, log.hyper.mu = cur.log.hyper.mu, sigma2.hyper.mu = sigma2.hyper.mu,
                                         hyper.mu_fixed=hyper.mu_fixed, conjugate.hyper.mu=conjugate.hyper.mu, family=family)
    
    if(conjugate.hyper.mu==FALSE){
    
    cur.log.post.hyper.mu<- log_post_hyper.mu(A=A, mu=cur.mu, log.hyper.mu = cur.log.hyper.mu, hyper.mu_fixed=hyper.mu_fixed, s.fixed=s.fixed, family = family)
    prop.log.post.hyper.mu<- log_post_hyper.mu(A=A, mu=cur.mu, log.hyper.mu = prop.log.hyper.mu, hyper.mu_fixed=hyper.mu_fixed, s.fixed=s.fixed, family = family)
    
    log.ratio.hyper.mu<-prop.log.post.hyper.mu-cur.log.post.hyper.mu
    
    if (log(runif(1)) < log.ratio.hyper.mu & is.na(log.ratio.hyper.mu)==FALSE){
      cur.samples.log.hyper.mu <- prop.log.hyper.mu
      rate.hyper.mu<-rate.hyper.mu+1
    } else{
      cur.samples.log.hyper.mu<- cur.log.hyper.mu
    }
     # updating the tuning parameters
    sigma2.hyper.mu<-adpative_function(index_MCMC_iter=i, sigma2_adapt=sigma2.hyper.mu, target_accept=0.40, 
                                       rate_adapt=rate.hyper.mu, burning1=burning1, burning2=burning2, adapt_seq=hyper.mu_adapt_seq2, 
                                       adapt=adapt, adpat_param=theta_hyper.mu, lower.acc=0.30, upper.acc=0.50)
    
    } else{
      cur.samples.log.hyper.mu<-prop.log.hyper.mu
      
    }
    ###### Proposing new parameters for beta (Gibbs steps) ######
      if(model.base==TRUE){
        cur.samples.beta=0
      } else {
    cur.samples.beta<-beta_sim(mu=cur.mu, intercept2 = cur.samples.intercept2, kappa_mu = cur.samples.kappa_mu,
                               beta2=cur.samples.beta2, W1=cur.samples.w1, W2=cur.samples.w2, Z2=Z2, A2=A2, t.A2.A2 = t.A2.A2)
      }
      #cur.samples.log.hyper.mu[1]<-log(5)
    #cur.samples.beta<-beta
    ###### Proposing new parameters for intercept1 (Gibbs steps) ######
    cur.samples.intercept1<- intercept1_sim(eta = cur.eta, beta1 = cur.samples.beta1,
                                            W1=cur.samples.w1, kappa_eta = cur.samples.kappa_eta, Z1=Z1, A1=A1)
    # alternative way to estimate intercept2
    #cur.samples.intercept1<-mean(cur.eta-Z1%*%cur.samples.beta1-A1%*%cur.samples.w1)
    
    ###### Proposing new parameters for intercept2 (Gibbs steps) ######
    cur.samples.intercept2<- intercept2_sim(mu=cur.mu, beta2=cur.samples.beta2, beta=cur.samples.beta,
                                            W1=cur.samples.w1, W2=cur.samples.w2, kappa_mu = cur.samples.kappa_mu, Z2=Z2, A2=A2)
    ## alternative way to estimate intercept2
    #cur.samples.intercept2<-mean(cur.samples.mu-Z2%*%cur.samples.beta2-cur.samples.beta*(A2%*%cur.samples.w1)-A2%*%cur.samples.w2)
    
    ###### Proposing new parameters for beta1 (Gibbs steps) ######
    cur.samples.beta1<- beta1_sim(eta = cur.eta, intercept1 = cur.samples.intercept1, W1=cur.samples.w1,
                                  kappa_eta = cur.samples.kappa_eta, Z1=Z1, A1=A1)
    
    ###### updating beta2 using Gibbs #########
    cur.samples.beta2<-beta2_sim(mu=cur.mu, intercept2=cur.samples.intercept2, beta=cur.samples.beta, kappa_mu=cur.samples.kappa_mu, 
                                 W1=cur.samples.w1, W2=cur.samples.w2, Z2=Z2, A2=A2)
    
    ###### updating mu using MALA #########
    prop.mu<-prop_mu(A=A, mu=cur.mu, intercept2 = cur.samples.intercept2, W1=cur.samples.w1, W2=cur.samples.w2, 
                          beta=cur.samples.beta, beta2=cur.samples.beta2,  kappa_mu = cur.samples.kappa_mu,
                          log.hyper.mu = cur.samples.log.hyper.mu, Z2=Z2, A2=A2, sigma2.mu = sigma2.mu, conjugate.mu=conjugate.mu, s.fixed=s.fixed, family=family)
    
    if(conjugate.mu==FALSE){
    cur.log.post.mu<-log_post_mu(A=A, mu=cur.mu, intercept2 = cur.samples.intercept2, W1=cur.samples.w1, W2=cur.samples.w2, 
                                 beta=cur.samples.beta, beta2=cur.samples.beta2,  kappa_mu = cur.samples.kappa_mu,
                                 log.hyper.mu = cur.samples.log.hyper.mu, Z2=Z2, A2=A2, s.fixed=s.fixed, family=family)
    
    prop.log.post.mu<-log_post_mu(A=A, mu=prop.mu, intercept2 = cur.samples.intercept2, W1=cur.samples.w1, W2=cur.samples.w2, 
                                  beta=cur.samples.beta, beta2=cur.samples.beta2,  kappa_mu = cur.samples.kappa_mu,
                                  log.hyper.mu = cur.samples.log.hyper.mu, Z2=Z2, A2=A2, s.fixed=s.fixed, family=family)
    
    cur.mala.log.dens.mu<-MALA_dens_mu(A=A, mu=prop.mu, mu.star = cur.mu, intercept2 = cur.samples.intercept2, W1=cur.samples.w1, W2=cur.samples.w2, 
                                       beta=cur.samples.beta, beta2=cur.samples.beta2,  kappa_mu = cur.samples.kappa_mu,
                                       log.hyper.mu = cur.samples.log.hyper.mu, Z2=Z2, A2=A2, sigma2.mu = sigma2.mu, s.fixed=s.fixed, family=family)
    
    prop.mala.log.dens.mu<-MALA_dens_mu(A=A, mu=cur.mu, mu.star = prop.mu, intercept2 = cur.samples.intercept2, W1=cur.samples.w1, W2=cur.samples.w2, 
                                        beta=cur.samples.beta, beta2=cur.samples.beta2,  kappa_mu = cur.samples.kappa_mu,
                                        log.hyper.mu = cur.samples.log.hyper.mu, Z2=Z2, A2=A2, sigma2.mu=sigma2.mu, s.fixed=s.fixed, family=family)
    
    log.ratio.mu<-prop.log.post.mu + cur.mala.log.dens.mu - cur.log.post.mu - prop.mala.log.dens.mu
    
    if (log(runif(1)) < log.ratio.mu & is.na(log.ratio.mu)==FALSE){
      cur.samples.mu <- prop.mu
      rate.mu<-rate.mu+1
    } else{
      cur.samples.mu<- cur.mu
    }
    # updating the tuning parameters
    sigma2.mu<-adpative_function(index_MCMC_iter=i, sigma2_adapt=sigma2.mu, target_accept=0.57, 
                                       rate_adapt=rate.mu, burning1=burning1, burning2=burning2, adapt_seq=mu_adapt_seq2, 
                                       adapt=adapt, adpat_param=theta_mu, lower.acc=0.50, upper.acc=0.65)
    }else{
       cur.samples.mu<- prop.mu
   }
    
    
    ###### updating eta using MALA #########
    prop.eta<-MALA_prop_eta(Y=Y, eta = cur.eta, intercept1 = cur.samples.intercept1, beta1 = cur.samples.beta1, 
                            kappa_eta = cur.samples.kappa_eta, W1=cur.samples.w1, Z1=Z1, A1=A1, sigma2_eta = sigma2_eta)
    
    cur.log.posterior.eta<-log.post_eta(Y=Y, eta = cur.eta, intercept1 = cur.samples.intercept1, beta1 = cur.samples.beta1, 
                                        kappa_eta = cur.samples.kappa_eta, W1=cur.samples.w1, Z1=Z1, A1=A1)
    
    prop.log.posterior.eta<-log.post_eta(Y=Y, eta = prop.eta, intercept1 = cur.samples.intercept1, beta1 = cur.samples.beta1, 
                                         kappa_eta = cur.samples.kappa_eta, W1=cur.samples.w1, Z1=Z1, A1=A1)
    
    cur.log.proposal.mala.eta<-MALA_dens_eta(Y=Y, eta = prop.eta, eta.star = cur.eta, intercept1 = cur.samples.intercept1, beta1=cur.samples.beta1,
                                             kappa_eta = cur.samples.kappa_eta, W1=cur.samples.w1, Z1=Z1, A1=A1, sigma2_eta = sigma2_eta)
    
    prop.log.proposal.mala.eta<-MALA_dens_eta(Y=Y, eta = cur.eta, eta.star = prop.eta, intercept1 = cur.samples.intercept1, beta1=cur.samples.beta1,
                                              kappa_eta = cur.samples.kappa_eta, W1=cur.samples.w1, Z1=Z1, A1=A1, sigma2_eta = sigma2_eta)
    
    log.ratio.eta<- prop.log.posterior.eta + cur.log.proposal.mala.eta - cur.log.posterior.eta - prop.log.proposal.mala.eta
    
    
    if (log(runif(1)) < log.ratio.eta & is.na(log.ratio.eta)==FALSE){
      cur.samples.eta <- prop.eta
      rate.eta<- rate.eta+1
    } else{
      cur.samples.eta<- cur.eta
    }
    
    sigma2_eta<-adpative_function(index_MCMC_iter=i, sigma2_adapt=sigma2_eta, target_accept=0.57, 
                                 rate_adapt=rate.eta, burning1=burning1, burning2=burning2, adapt_seq=eta_adapt_seq2, 
                                 adapt=adapt, adpat_param=theta_eta, lower.acc=0.50, upper.acc=0.65)
    
    ###### updating W1 using Gibbs #########
    cur.samples.w1<-W1_sim(eta = cur.samples.eta, mu=cur.samples.mu, intercept1 = cur.samples.intercept1, intercept2 = cur.samples.intercept2,
                           beta1 = cur.samples.beta1, beta2 = cur.samples.beta2, W2=cur.samples.w2, kappa_w1 = cur.samples.kappa_w1, 
                           kappa_eta=cur.samples.kappa_eta, kappa_mu=cur.samples.kappa_mu, beta=cur.samples.beta, Q=Q, 
                           Z1=Z1, Z2=Z2, A1=A1, A2=A2, t.A1.A1=t.A1.A1, t.A2.A2=t.A2.A2,  t.one.one, tau.w1_sum.const=0.0001)
    cur.samples.w1<- cur.samples.w1-mean(cur.samples.w1) ## making sure that sum of W1 is zero
    ###### updating W2 using Gibbs #########
    cur.samples.w2<-W2_sim(mu=cur.samples.mu, intercept2 = cur.samples.intercept2, beta2=cur.samples.beta2, W1=cur.samples.w1, 
                           kappa_w2 = cur.samples.kappa_w2, kappa_mu = cur.samples.kappa_mu, beta = cur.samples.beta, 
                           Q=Q, Z2=Z2, A2=A2, t.A2.A2 = t.A2.A2, t.one.one, tau.w2_sum.const=0.0001)
    cur.samples.w2<- cur.samples.w2-mean(cur.samples.w2) ## making sure that sum of W2 is zero
    
    ### Saving the samples after thinning the samples at every by iterations .  
    if((i%%by)-1==0){
      samples[j,]<- c(cur.samples.log.hyper.mu,
                      cur.samples.kappa_w1, 
                      cur.samples.kappa_w2,
                      cur.samples.kappa_eta,
                      cur.samples.kappa_mu,
                      cur.samples.intercept1,
                      cur.samples.intercept2,
                      cur.samples.beta,
                      cur.samples.beta1,
                      cur.samples.beta2,
                      cur.samples.eta[1:2],
                      cur.samples.w1[1], 
                      cur.samples.mu[1],
                      cur.samples.w2[1])
      j=j+1
    }
    ### storing the adaptive tuning parameters
    if((i%%adapt)-1==0){ # to save allexp the scale parameter of the MALA
      sigma.matrix[m,]<- c(sigma2_eta, sigma2.mu, sigma2.hyper.mu)
      m=m+1
    }
  }
  return(list("samples"=samples, ### saving the sample for the hyperparameters to see the traceplots
              "imputed.Y.WSD"=imputed.Y.WSD, "imputed.A.WSD"=imputed.A.WSD, "imputed.Y.WSD.samples"=imputed.Y.WSD.samples, "imputed.A.WSD.samples"=imputed.A.WSD.samples,
              "imputed.Y.OSD"=imputed.Y.OSD, "imputed.A.OSD"=imputed.A.OSD, "imputed.Y.OSD.samples"=imputed.Y.OSD.samples, "imputed.A.OSD.samples"=imputed.A.OSD.samples,
              "post.sum.mu"=post.sum.mu, "post.sum.eta"=post.sum.eta, "post.sum.mu.samples"=post.sum.mu.samples, "post.sum.eta.samples"=post.sum.eta.samples,
              "tuning_param_x_hyper"=sigma.matrix, ## tuning parameters values
              "Acc.rate eta"=rate.eta/(N.MCMC-(burning1+burning2)))) ## acceptance rate for eta
}

