#####################################################################################################
########## Adaptive function for updating the variance componenets of MH and MALA algorithm #########
#####################################################################################################
#' @index_MCMC_iter: index of MCMC iterations
#' @sigma2_adapt: the tunuining parameter in the MH and MALA algorithm
#' @target_accept: target acceptance probabilty, which is 0.224 for random wala MH and 0.57 for MALA\
#' @rate_adapt: acceptance rate in MCMC iterations calculated at every fixed number of iterations
#' @burning1: the first burning period in MCMC algorithm
#' @burning2: the second burning periord in MCMC algorithm
#' @adapt_seq: the sequence, at which we update the @sigma2_adapt
#' @adpat: the number of MCMC iterations after which we upadte the @sigmna2_adapt
#' @adapt_param: the scale parameter in adaptive algorithm
#' @lower.acc: The lower acceptance probability for RWM and MALA which is 0.15, and 0.50, respectively
#' @upper.acc: The lower acceptance probability for RWM and MALA which is 0.30, and 0.65, respectively
adpative_function<-function(index_MCMC_iter, sigma2_adapt, target_accept, rate_adapt, burning1, burning2, adapt_seq, adapt, adpat_param, lower.acc, upper.acc){
  if(index_MCMC_iter %in% adapt_seq){
  if (index_MCMC_iter<(burning1+1)){
      sigma2_adapt<- exp(((rate_adapt/adapt)-target_accept)/adpat_param)*sigma2_adapt
    }
    if (((burning1+1) < index_MCMC_iter)  & (index_MCMC_iter< (burning1+burning2+2))){
    if (((rate_adapt/(adapt))>upper.acc) | ((rate_adapt/(adapt))<lower.acc)){ #burning+2 to calculate the accpertance rate after the burning samples
        sigma2_adapt<- exp(((rate_adapt/adapt)-target_accept)/adpat_param) * sigma2_adapt
      }
    }
    } else {
      sigma2_adapt<-sigma2_adapt
    }
  return(sigma2_adapt)
}

