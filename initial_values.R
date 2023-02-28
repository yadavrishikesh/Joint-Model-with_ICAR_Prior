######################################################################################
##########  Initial values and other stuff needed for MCMC steps##############
######################################################################################
init_fun_log.hyper.mu<-function(family){
  if(family=="Gamma-gamma"){
    log.hyper.mu<-c(0.1,0.1)
    log.hyper.mu.name<-c(expression(log(gamma)[1]), expression(log(gamma)[2]))
  } else if (family=="extended-GPD"){
    log.hyper.mu<-c(0.1,0.1)
    log.hyper.mu.name<-c(expression(log(k)[1]), expression(log(xi)[2]))
  } else if (family=="GPD"){
    log.hyper.mu<-c(0.1)
    log.hyper.mu.name<-c(expression(log(xi)))
  } else if(family=="Gamma"){
    log.hyper.mu<-c(0.1)
    log.hyper.mu.name<-c(expression(log(shape)))
  } else if (family=="log-Gamma"){
    log.hyper.mu<-c(2)
    log.hyper.mu.name<-c(expression(log(shape)))
  } else if(family=="Weibull"){
    log.hyper.mu<-c(0.1) 
    log.hyper.mu.name<-c(expression(log(shape)))
  } else if(family=="log-normal"){
    log.hyper.mu<-c(0.1) 
    log.hyper.mu.name<-c(expression(log(kappa)))
  }
  else if(family=="Burr"){
    log.hyper.mu<-c(0.1, 0.1) 
    log.hyper.mu.name<-c(expression(log(k)), expression(log(c)))
  }  else if(family=="generalized-Gamma"){
    log.hyper.mu<-c(0.1, 0.1) 
    log.hyper.mu.name<-c(expression(log(k)), expression(log(c)))
  }
  return(list(log.hyper.mu.init=log.hyper.mu, log.hyper.mu.name=log.hyper.mu.name))
}

#################################################################################################################
####### Initial values of all the other parameters, that does not depends on the choice of mark distribuitions ###
##################################################################################################################
init_fun_all_other_param<-function(Z1, Z2, A, Y, Q){
kappa_w1_init<-1
kappa_w2_init<-1
kappa_eta_init<-1
kappa_mu_init<-1
intercept1.init<-0.1
intercept2.init<-0.1
beta.init<-0.1
beta1.init<-rep(0.1, ncol(Z1))
beta2.init<-rep(0.1, ncol(Z2))
eta_init<-rep(0.01, times=length(Y))
W1_init<-rep(0.01, times=ncol(Q))
mu_init<-rep(0.01, times=length(A))
w2_init<-rep(0.01, times=ncol(Q))
init1<-c(kappa_w1_init, kappa_w2_init, kappa_eta_init, kappa_mu_init,
         intercept1.init, intercept2.init,
         beta.init, beta1.init, beta2.init, eta_init, W1_init, mu_init,  w2_init)
model.param.name<-c(expression(kappa[w1]), expression(kappa[w2]), expression(kappa[eta]), expression(kappa[mu]),
                    "intercept1","intercept2", expression(beta), paste0("beta1_",1:ncol(Z1)), paste0("beta2_",1:ncol(Z2)), 
                    expression(eta[1]), expression(eta[2]), expression(W1[1]), expression(mu[1]), expression(W2[1]))

return(list(init.all.other.param=init1, model.param.name.all.other.param=model.param.name))
}

