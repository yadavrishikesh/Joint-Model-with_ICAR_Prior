rm(list = ls())
## SET THE WORKING DIRECTORY
setwd(this.path::here())
######################################################################################
##########  loading the data and source code ##############
######################################################################################
load("A1_projection.matrix.Rdata") #' @projection matrices that project slop unit to pixel levels for counts (based on real landslides datasets)
load("A2_projection.matrix.Rdata") #' @projection matrices that project slop unit to pixel levels for size (based on real landslides datasets)
load("landslides_data_information.Rdata") #' @data loading the precision matrices, adjacency information and the covariates Z1 and Z2 corresponding to the count and sizes  

################################################################################
######### deciding about the families and types of simulation we needed  ######
################################################################################
families<-c("Gamma-gamma","extended-GPD", "GPD","generalized-Gamma","Gamma","log-Gamma", "Weibull", "Burr", "log-normal") ## name of family
models<-c("Model1","Model2", "Model3", "Model4") ## the four different models to fit
CVs<-c("OSD","WSD")  ## types of cross-validation schemes; OSD: out of sample diagnostics, WSD: within sample diagnostics
model.bases<-c(TRUE, FALSE) ## whether the model is baseline (beta=0) or not
conjugate.hyper.mus<- c(TRUE, FALSE) ## whether we are using the conjugacy of hyperparameter corresponding to size process or not (TRUE only for log-Normal, otherwise FALSE)
conjugate.mus= c(TRUE, FALSE) ## whether we are using the conjugacy of median process parameter (mu) in the size or not (TRUE only for log-Normal, otherwise FALSE)

################# For Gamma-Gamma Simulation set ups
### fixing the model parameter values to run specific model (=1 means the Gamma-Gamma mark for size process)
family.ind=1  ## index that decide which mark distribution would be used as a mark distribution in the landslide size process
model.ind=1 ## index that decide which model is to be fitted (=1 means model 1)
CV.ind=2  ## index that decide which cross-validation scheme to be used, (=2 means within sample)
k_fold_CV=10 ## the number of foldes in the OSD setting
CV1_number=1 ## a indicator number that keep record of which cross-validation number is in the OSD setting (a number between 1 to k_fold_CV)
conjugate.hyper.mu.no=2 ## index that decide whether we are using conjuagte distribution for the hyperparameters on mark distribution(=1 only for log-Normal otherwise this would be 2)
conjugate.mu.no=2 #index that decide whether we are using conjuagte distribution for the median of size processes(=1 only for log-Normal otherwise this would be 2)
model.base.no=2 ## a index that decide whether we are fitting the baseline model (beta=0); (=2 means this is not baseline model)

################# For log-Normal simulation set ups: for this you need to change the just three below parameters
family.ind=9  ## index that decide which mark distribution would be used as a mark distribution in the landslide size process
conjugate.hyper.mu.no=1 ## index that decide whether we are using conjugate distribution for the hyperparameters on mark distribution or not (=1 only for log-Normal otherwise this would be 2)
conjugate.mu.no=1

#### extracting the corresponding model parameters: 
family=families[family.ind]
model<-models[model.ind]
CV<-CVs[CV.ind]
model.base<- model.bases[model.base.no]
conjugate.hyper.mu<-conjugate.hyper.mus[conjugate.hyper.mu.no]
conjugate.mu<-conjugate.mus[conjugate.mu.no]

#########################################################################
##########  Simulating the Data #########################################
#########################################################################
source("simulating_data.R") ## loading the function that simulate from the joint mark distributions for given family of distributions and model parameters: We simulate from model M1 but with different mark distributions listed in Table 1 in manuscript
### Fixed parameters
adjac_inform<-fun_proj_matrix(Q)
kappa_w1<- 1 # precision parameter in CAR prior W1
kappa_w2<- 1  #precision parameter in CAR prior W2
kappa_eta<- 0.5 #precision parameter in independent random effects \varepsilon_{\eta}
kappa_mu<- 0.5 #precision parameter in independent random effects \varepsilon_{\mu}
intercept1<- -0.5 ## intercept for counts
intercept2<- 0.5 ## intercept for sizes
beta<- 1 ## shared parameter between the count and size processes
other.hyper<- c(kappa_w1, kappa_w2, kappa_eta, kappa_mu, intercept1, intercept2,  beta) ### collecting all the parameteres into a vectors
beta1<-c(rep(0.2,3), rep(-0.2,3), rep(0.15,5))  ## fixed covariates coefficient for counts
p<-length(beta1)  
beta2<-c(rep(0.15,3), rep(-0.1,3), rep(0.2,5)) ## fixed covariates coefficient for sizes
q<-length(beta2)

#### Hyperprameters associated with all the 9 Mark distributions: for more details about hyperparameters of the mark distributions see Table 1 in main manuscript 
fun.hyper.mu<- function(family){
  if(family=="Gamma-gamma"){
    log.hyper.mu<-log(c(10, 10))  ## the true values in log-scale for the Gamma-Gamma mark distribution
  } else if (family=="extended-GPD"){
    log.hyper.mu<-c(log(5), log(0.5)) ## the true values in log-scale for the extended-GPD mark distribution
  } else if (family=="GPD"){
    log.hyper.mu<-log(0.5) ## the true values in log-scale for the GPD mark distribution
  } else if(family=="generalized-Gamma"){
    log.hyper.mu<-log(c(10, 10)) ## the true values in log-scale for the generalized-Gamma mark distribution
  } else if(family=="Gamma"){
    log.hyper.mu<-log(c(10)) ## the true values in log-scale for the Gamma mark distribution
  } else if(family=="log-Gamma"){
    log.hyper.mu<-log(c(10)) ## the true values in log-scale for the log-Gamma mark distribution
  } else if(family=="Weibull"){
    log.hyper.mu<-log(c(10)) ## the true values in log-scale for the Weibull mark distribution
  } else if(family=="Burr"){
    log.hyper.mu<-log(c(10, 10)) ## the true values in log-scale for the Burr mark distribution
  } else if(family=="log-normal"){
    log.hyper.mu<- 5 ## the true values in for the log-normal mark distribution
  }
  return(log.hyper.mu)
}

source("Sim_singular_MVN_precision.R")  ## R functions that simulate from the singular multivariate distribution
source("Other_functions.R") ## other R functions that includes, for example the density of the extended-GPD and generalized-Gamma distributions
#' @sim_mark_function: function that simulate the landslides counts (Y) and sizes (A) given the model parameters; see simulating_data.R for the details about this function
set.seed(1)
Sim_mark_data<-sim_mark_function(Q=Q,
                                 other.hyper=other.hyper,
                                 beta1=beta1, beta2=beta2,
                                 Z1=as.matrix(Z1), Z2=as.matrix(Z2), A1=as.matrix(A1), A2=as.matrix(A2), family=family)

Y<-Sim_mark_data$Y  ## the simulated counts 
A<-Sim_mark_data$A  ## the simulated sizes 

### saving the true values to compare it with the posterior predictive distributions
true.values<-c(fun.hyper.mu(family = family), other.hyper, beta1, beta2, Sim_mark_data$eta[1:2], Sim_mark_data$W1[1], Sim_mark_data$mu[1], Sim_mark_data$W2[2])


######## Uploading all the other important R functions 
source("update.param.R")  #' @updata.parameters to update the model parameters:
source("likel_and_grad.R") #' @likelhood_and_gradients likelihood and gradients for the simulation using MALA and Metropolis-Hasting steps
source("initial_values.R") #' @initial_values R function that contain the simulation of the initial values for the model hyperparameters
source("MCMC_function_CV.R") #' @main_MCMC_function this is the most important MCMC function where we call several other function to finally draw Markov samples 
source("update_tuning_param_function.R") #' @update_tuning update tuning parameters of the adaptive MALA and adaptive Random walk algorithm (RWM)
source("Imputation_at_missing_locations.R") #' @update_tuning ### imputation at unobserved locations: needed to impute at locations where we don't have observations: useful in  OSD cross-validation setting

#########Main MCMC parameters ##############
N.MCMC<- 1*10^5 ## Total Number of MCMC iterations
adapt<-500  ## after this number of MCMC iterations, we update the tuning parameters of MALA and RWM
eta_adapt_seq2<-seq(from=adapt, to=N.MCMC, by=3*adapt) ##  sequences where we adapt the eta parameters (intensity of Poisson distribution) 
mu_adapt_seq2<-seq(from=2*adapt, to=N.MCMC, by=3*adapt) ##  sequences where we adapt the mu parameters (median of the Marked process) 
hyper.mu_adapt_seq2<-seq(from=3*adapt, to=N.MCMC, by=3*adapt) ##  sequences where we adapt the hyperparameters (hyper.mu) of the mark distribution (\theta_A)  
sigma2_eta<-0.0001 ## tuning parameters for eta 
sigma2.mu<-0.001  ##  tuning parameters for mu
sigma2.hyper.mu<-0.001 #  tuning parameters for hyper.mu
theta_eta<-0.4  ## scale in the adaptive function for eta
theta_mu<-0.4 ## scale in the adaptive function for mu 
theta_hyper.mu<-0.4 ## scale in the adaptive function for hyper.mu
by<-100  ## the thinning steps 

fixed.kappa_mu<-1000 ## fixed precision hyperparameters in the independent random effects \varepsilon_{\mu} (useful for the models M2, M3, and M4)
fixed.kappa_eta<-1000 ## fixed precision hyperparameters in the random effects \varepsilon_{\eta} (useful for the models M2, M3, and M4)

## burning samples
burning1<-floor(N.MCMC/4)  ## burning1 is the number of iteration we use for first burn in period 
burning2<-floor(N.MCMC/2) ## burning2 is the number of iteration we use for second burn in period 
burnin<-burning1+burning2 ## total number of burn-in samples


#### Fixed hyperparameters in all of the hyperpriors
m.bar<- 5.39 ## avaerage number of neighbour, to adjust the priors in the likelihood: see the manuscript Section 2.3 for more details about the prior choice
hyper_fixed<-c(0.75, 3, # kappa_eta ~ Gamma(shape=hyper_fixed[1], rate=hyper_fixed[2])
               0.75, 3/(m.bar*0.7^2), # kappa_w1 ~ Gamma(shape=hyper_fixed[3], rate=hyper_fixed[4])
               0.75, 3/(m.bar*0.7^2), # kappa_w2 ~ Gamma(shape=hyper_fixed[5], rate=hyper_fixed[6])
               0.75, 3, # kappa_mu ~ Gamma(shape=hyper_fixed[7], rate=hyper_fixed[8])
               0.01,     # beta ~ norm(mean=0, precision=hyper_fixed[13]))
               0.01,     # beta1 ~ mvtnorm_p(mean=0, precision=hyper_fixed[14] I_p))
               0.01,     # beta2 ~ mvtnorm_q(mean=0, precision=hyper_fixed[15] I_q))
               0.01,     # intercept1~ rnorm(1,mean=0, precision=hyper_fixed[16]) 
               0.01)     # intercept2~ rnorm(1,mean=0, precision=hyper_fixed[17])
hyper.mu_fixed<- c(0.25,0.25)  ## hyperparameters in the hyper.mu parameters 


######### Divide the datasets for the OSD setting 
set.seed(1)
size.subset<- floor(length(unique(idx_car_cox))/k_fold_CV) ## number of slope units that should be on average in each of the fold of OSD setting
SU.rand<- sample(unique(idx_car_cox)) ## randomly arranging all the slope units
subset.SU<- split(SU.rand, ceiling(seq_along(SU.rand)/size.subset)) ### slope units in each foldes 
### check if all the slop unite are selected and its disjoint
#sort(unique(stack(subset.SU)$values), decreasing = FALSE)
##stack(subset.SU)$ind
if(CV=="OSD"){ ## creating the data when the setting is OSD
  subset.SU.stak<-stack(subset.SU)
  sample.SU<- stack(subset.SU)$values[stack(subset.SU)$ind==CV1_number]
  Y_data_frame_cox<-data.frame(cbind(idx_car_cox=idx_car_cox, original_Y=Y, Y_with_NA=Y))
  Y_data_frame_cox[idx_car_cox%in%sample.SU,3]<-NA
  
  A_data_frame_size<-data.frame(cbind(idx_car_size=idx_car_size, original_A=A, A_with_NA=A))
  A_data_frame_size[idx_car_size%in%sample.SU,3]<-NA
  colnames(A_data_frame_size)<-c("idx_car_size", "original_A", "A_with_NA")
} else if (CV=="WSD"){ ## creating the data when the setting is WSD
  Y_data_frame_cox<-data.frame(cbind(idx_car_cox=idx_car_cox, original_Y=Y, Y_with_NA=Y))
  A_data_frame_size<-data.frame(cbind(idx_car_size=idx_car_size, original_A=A, A_with_NA=A))
  colnames(A_data_frame_size)<-c("idx_car_size", "original_A", "A_with_NA")
}

##################################################################################
############################ Running the MCMC  algorithm #############################
##################################################################################

start_time<-Sys.time()
set.seed(1)
#' @MCMC.main_function: 
#' @Iinput this function take input as counts Y, sizes A, design matrices Z1, Z2, projection matrices A1, A2 and several other information of the MCMC setting as well as the model related information 
#' @Output The posterior samples for each of the model parameters with a thinning size equal to by
MCMC1<-MCMC.main_function(N.MCMC = N.MCMC, Y = Y_data_frame_cox$Y_with_NA, ind_NA_Y=is.na(Y_data_frame_cox$Y_with_NA), 
                          A = A_data_frame_size$A_with_NA, ind_NA_A=is.na(A_data_frame_size$A_with_NA), 
                          Q = Q, Z1 = as.matrix(Z1), Z2 = as.matrix(Z2), by = by, adapt = adapt, burning1 = burning1, 
                          burning2 = burning2, s.fixed=1, sigma2_eta = sigma2_eta, theta_eta = theta_eta, 
                          sigma2.mu = sigma2.mu, theta_mu = theta_mu, sigma2.hyper.mu = sigma2.hyper.mu, 
                          theta_hyper.mu = theta_hyper.mu, A1 = as.matrix(A1), A2 = as.matrix(A2), family = family, CV1_number = CV1_number, 
                          k_fold_CV=k_fold_CV, print.result = TRUE,  traceplot = TRUE, conjugate.hyper.mu = conjugate.hyper.mu, 
                          conjugate.mu = conjugate.mu, model=model, model.base=model.base, fixed.kappa_mu=fixed.kappa_mu, 
                          fixed.kappa_eta=fixed.kappa_eta,CV=CV, true.values=true.values, max.A=max(A_data_frame_size$original_A, na.rm = TRUE),  simulation=TRUE)
end_time<-Sys.time()
run_time<-end_time-start_time
#### saving the outputs of MCMC  
# name.file<-paste0("RData_simulation_results/Mydata1_",family, "_",CV1_number,"_", CV,model,"_BaseModel_",model.base,".Rdata")
# save.image(name.file)



