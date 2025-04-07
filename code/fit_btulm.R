library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)
library(foreach)
library(doParallel)
library(Matrix)
library(microbenchmark)
library(mvtnorm)
library(LaplacesDemon)
library(BayesLogit)
source("code/models/BTULM.R")

##########################
### Simulation Settings ##
##########################
set.seed(555)

gibbs.runs <- 2000
gibbs.burn <- 1000
discrete<-T

save_dir <- "data/"

load(paste0(save_dir, "empirical_samples.RData"))
Nsim <- length(samples)

numCores <- 2

cl<-makeCluster(numCores, type="FORK", outfile="")
registerDoParallel(cl)

pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")

btulm_fit <-foreach(i=1:Nsim, .combine="rbind", .verbose=TRUE,
	              .packages=c('dplyr', 'microbenchmark')) %dopar% {
############################
### Generate Sample Data ###
############################
  set.seed(i + 15) 

  sigma2_beta <- 10000
  ### Process sample
  HPS_sample_long <- samples[[i]]

  areas <- HPS_sample_long$AREA
  n_areas <- nlevels(areas)
  Psi <- Matrix(model.matrix(~areas-1))
  colnames(Psi) <- c()

  covar_names <- HPS_pop_long %>% select(starts_with("COVAR")) %>%
                             names() %>%
                             paste(collapse="+")
  covar_formula <- as.formula(paste("~", covar_names))
  
  X <- Matrix(model.matrix(covar_formula, data=HPS_sample_long))
  colnames(X) <- c()
  weights <- HPS_sample_long$SCALE_WEIGHT
  weeks <- HPS_sample_long$WEEK
 
  model_data <- list(X=X, Y=HPS_sample_long$RESPONSE, Psi=Psi, 
                     sigma2_beta=sigma2_beta, weights=weights, 
                     weeks=weeks, areas=areas)
  
  ### Fit Gibbs
  gibbs_dep_data <- model_data
  gibbs_dep_data$n_gibbs <- gibbs.runs
  gibbs_dep_data$n_burn <- gibbs.burn

  gibbs_dep_runtime <- microbenchmark(
                       gibbs_dep_out <- do.call(fit_dep_gibbs, gibbs_dep_data),
                    times = 1
                   )

  print(paste("Dependent Gibbs sampler", i, "fit in", 
              gibbs_dep_runtime$time/(10^9), "seconds"))
  
  return(gibbs_dep_out)
}
stopCluster(cl)

save(btulm_fit, file=paste0(save_dir,"simulation_btulm_results_binary.RData"))

