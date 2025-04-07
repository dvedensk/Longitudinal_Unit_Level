library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)
library(foreach)
library(doParallel)
library(microbenchmark)
library(mvtnorm)
library(Matrix)
library(LaplacesDemon)
library(BayesLogit)
source("code/models/BBULM.R")

##########################
### Simulation Settings ##
##########################
set.seed(555)

Nsim <- 2
gibbs.runs <- 100
gibbs.burn <- 50
discrete<-T

save_dir <- "data/"

load(paste0(save_dir, "empirical_samples.RData"))
Nsim <- length(samples)

numCores <- 2

#Parallelization parameters 
cl<-makeCluster(numCores, type="FORK", outfile="")
registerDoParallel(cl)

pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")
week_covars <-  gsub("\\+COVAR.NEW", "", pop_covars)

pop_grouped <- HPS_pop_long %>%  
                  group_by(WEEK, AREA, COVAR.1, COVAR.2,
                           COVAR.3, COVAR.4, COVAR.NEW) %>% 
                  summarize(popsize=n()) 

bbulm_fit <-foreach(i=1:Nsim, .combine="rbind", .verbose=TRUE,
	              .packages=c('dplyr', 'microbenchmark')) %dopar% {
############################
### Generate Sample Data ###
############################
  set.seed(i + 15) 

  sigma2_beta <- 10000
  HPS_sample_long <- samples[[i]]

  areas <- HPS_sample_long$AREA
  n_areas <- nlevels(areas)
  n_weeks <- max(HPS_sample_long$WEEK)

  Psi <- Matrix(model.matrix(~areas-1))
  colnames(Psi) <- c()

  covar_names <- HPS_pop_long %>% select(starts_with("COVAR")) %>%
                             names() %>%
                             paste(collapse="+")
  covar_formula <- as.formula(paste("~", covar_names))
  
  gibbs_indep_results <- list()
  for(t in 1:n_weeks){
    print(paste("Fitting week", t))
    HPS_week <- filter(HPS_sample_long, WEEK==t) %>% select(-COVAR.NEW)
    X_week <- model.matrix(as.formula(paste("~", week_covars)), data=HPS_week)
    N_week <- nrow(HPS_week)
    scale_weights <- HPS_week$PWEIGHT * N_week/sum(HPS_week$PWEIGHT)

    gibbs_indep_results[[t]]  <- PLMBmcmc(Y=HPS_week$RESPONSE, X=X_week, 
                                     Psi=model.matrix(~AREA-1, data=HPS_week), 
                                     wgt=scale_weights, sig2b=sigma2_beta,
                                     iter=gibbs.runs, burn=gibbs.burn)
  }
  return(gibbs_indep_results)
}
stopCluster(cl)

save(bbulm_fit,
     file=paste0(save_dir,"simulation_bbulm_results_binary.RData"))

