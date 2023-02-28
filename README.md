## For more details about the mdoeling framework and RCode contact rishikesh.yadav@kaust.edu.sa

## Simulation study and inference for the: ``Joint modeling of landslide counts and sizes using spatial marked point processes with sub-asymptotic mark distributions (2023)''

This set of R codes simulate and do the inference of the Joint marked point process model defiend by Eqn. (2) and Eqn. (3) in the manuscript.

The description of all the attached files are provided below:
### running_MCMC_algorithm_simulation.R
This R file calls all the other function to run the MCMC and generate the final outputs  

### MCMC_function_CV.R
This is the main MCMC function that combines several other R functions to simulate the MCMC samples of each of the model parameters

### likel_and_grad.R
Likelihood and gradient for all the model and mark distributions 

### Imputation_at_missing_locations.R
Set fo function used to impute at the missing locations if needed

### initial_values.R
Set the initial values for all the model parameters

### simulating_data.R
R functions that simulate the data from the model detailed in Eqn. (2) and Eqn. (3) for certain combinations of model parameters 

### update_tuning_param_function.R
R function that tune the tuning parameters of the random walk Metropolis-Hastings (RWMH) and Metropolis-Hastings adjusted Langevin(MALA)

### update.param.R
R functions that are used as a proposal distribution for all the model parameters using either by RWMH or MALA 

### Sim_singular_MVN_precision.R
Functions to simulate from the multivaraite Gaussian density with precision parameterization 

### Other_functions.R 
Other helpful function 

### landslides_data_information.Rdata
Rdata consiing of useful information such as covariates and location infomations of the data (Note that the real data are not provided due to the confidentiality)

### A1_projection.matrix.Rdata
Projection matrix that project the slope unite (coarser) resolution to the pixel unit (finer) resolution  for counst data 

### A2_projection.matrix.Rdata
Projection matrix that project the slope unite (coarser) resolution to the pixel unit (finer) resolution  for counst data 


