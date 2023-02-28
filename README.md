## For more details about the mdoeling framework and RCode contact rishikesh.yadav@kaust.edu.sa

## Simulation study and inference for the: ``Joint modeling of landslide counts and sizes using spatial marked point processes with sub-asymptotic mark distributions (2023)''

This set of R codes simulate and do the inference of the Joint marked point process model defiend by Eqn. (2) and Eqn. (3) in the manuscript.

The description of all the attached files are provided below:
# Data_and_initial_values.R
R functions that simulate the data from the above product mixture model for certain combinations of model parameters and also set the initial values for all the model parameters

Data_nsites_20_ntime_100.Rdata
The R data that is resulting from some parameter combinations mentioned in Data_and_initial_values.R for 20 spatial locations and 100 independent temporal replicates

MCMC_main_function.R
This is the main MCMC function that combines several other R functions to simulate the MCMC samples of each of the model parameters

other-function.R
Other function that are used in MCMC_main_function.R

update_tuning_param_function.R
R function that tune the tuning parameters of the RWM and SGLD algorithms

update.param.R
R functions that are used as a proposal distribution for all the model parameters
