library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)
library(foreach)
library(doParallel)
library(microbenchmark)
library(mvtnorm)
library(ggplot2)
library(sampling)
library(Matrix)
source("code/helper_functions.R")
source("code/models/GBULM.R")

##########################
### Simulation Settings ##
##########################
set.seed(555)

Nsim <- 100
gibbs.runs <- 2000
gibbs.burn <- 500

numCores <- 2
save_dir <- "data/"

load(paste0(save_dir, "empirical_samples_gauss.RData"))

#Parallelization parameters 
cl<-makeCluster(numCores, type="FORK", outfile="")
registerDoParallel(cl)

pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")
n_weeks <- max(HPS_pop_long$WEEK)
n_areas <- nlevels(HPS_pop_long$AREA)

sim_results <-foreach(i=1:Nsim, .combine="rbind", .verbose=TRUE,
	              .packages=c('dplyr', 'microbenchmark')) %dopar% {
############################
### Generate Sample Data ###
############################
  set.seed(i + 15) 

  ### Process sample
  print(paste("Taking sample",i))
  HPS_sample_long <- samples[[i]] 

  start_time <- Sys.time()
  ##Fit non-time-dependent
  gibbs_indep_ests <- gibbs_indep_lows <- gibbs_indep_highs <- c()
  for(t in 1:n_weeks){
    HPS_week <- filter(HPS_sample_long, WEEK==t)
    X_week <- model.matrix(as.formula(paste("~", pop_covars)), data=HPS_week)
    N_week <- nrow(HPS_week)
    scale_weights <- HPS_week$PWEIGHT * N_week/sum(HPS_week$PWEIGHT)

    mod <- model1(y=HPS_week$RESPONSE, X=X_week, 
                  Psi=model.matrix(~AREA-1, data=HPS_week), 
                  wgt=scale_weights, iter=gibbs.runs, burn=gibbs.burn)
    popWeek <- filter(HPS_pop_long, WEEK==t)
    popXweek <- model.matrix(as.formula(paste("~", pop_covars)), data=popWeek)
    popPsiweek <- model.matrix(~AREA-1, data=popWeek)
    preds <- cbind(popXweek, popPsiweek) %*% t(mod$Beta)

    sigmas <- matrix(rep(sqrt(mod$sig2), nrow(popXweek)), 
                     nrow=nrow(popXweek), byrow=T)

    preds <- rnorm(n=length(preds),
                   mean=c(preds),
                   sd=c(sigmas))
    preds <- matrix(preds,
                    nrow=nrow(popXweek))
    
    #append to running total
    gibbs_indep_agg <- data.frame(AREA=popWeek$AREA, preds) %>%
                       group_by(AREA) %>%
                       summarize_all(mean) %>%
                       arrange(AREA)

    rm(popXweek, popPsiweek, popWeek, preds, sigmas)
    gc(verbose=F)

    gibbs_indep_ests <- c(gibbs_indep_ests, rowMeans(gibbs_indep_agg[,-1]))
    gibbs_indep_lows <- c(gibbs_indep_lows, apply(gibbs_indep_agg[,-1], 1, quantile, probs=0.025))
    gibbs_indep_highs <- c(gibbs_indep_highs, apply(gibbs_indep_agg[,-1], 1, quantile, probs=0.975))

  }

  gibbs_indep_results <- data.frame(AREA=gibbs_indep_agg$AREA, 
                                  WEEK=sort(rep(1:n_weeks, n_areas)),
                                  MEAN=gibbs_indep_ests,
                                  CI_LOWER=gibbs_indep_lows,
                                  CI_UPPER=gibbs_indep_highs)

  runtime <- Sys.time() - start_time                        
  return(gibbs_indep_results, runtime)
}
stopCluster(cl)

save(sim_results, HPS_pop_long,
     file=paste0(save_dir,"simulation_gbulm_results.RData"))

