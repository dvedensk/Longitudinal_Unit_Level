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
source("code/models/GTULM.R")

##########################
### Simulation Settings ##
##########################
set.seed(555)

gibbs.runs <- 2000
gibbs.burn <- 500
 
save_dir <- "data/"

load(paste0(save_dir, "empirical_samples_gauss.RData"))
Nsim <- length(samples)

numCores <- 2

cl<-makeCluster(numCores, type="FORK", outfile="")
registerDoParallel(cl)

gtulm_fit <-foreach(i=1:Nsim, .combine="rbind", .verbose=TRUE,
	              .packages=c('dplyr', 'microbenchmark')) %dopar% {
############################
### Generate Sample Data ###
############################
  set.seed(i + 15) 

  sigma2_beta <- 10000
  ### Process sample
  gibbs_dep_data <- params[[i]]

  ### Fit Gibbs
  gibbs_dep_data$n_gibbs <- gibbs.runs
  gibbs_dep_data$n_burn <- gibbs.burn

  gibbs_dep_runtime <- microbenchmark(
                       gibbs_dep_out <- do.call(fit_dep_gibbs, gibbs_dep_data),
                    times = 1
                   )

  print(paste("Dependent Gibbs sampler", i, "fit in", 
              gibbs_dep_runtime$time/(10^9), "seconds"))
  
  return(gibbs_dep_out, gibbs_dep_runtime)
}
stopCluster(cl)

save(gtulm_fit, file=paste0(save_dir,"simulation_gtulm_results.RData"))

