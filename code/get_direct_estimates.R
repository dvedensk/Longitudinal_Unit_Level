library(dplyr)
library(mase)
source("code/helper_functions.R")
case <- "gaussian"

save_dir <- "data/"

if(case=="binary"){
  load(paste0(save_dir, "empirical_samples.RData"))  
}else{
  load(paste0(save_dir, "empirical_samples_gauss.RData"))  
}
Nsim <- length(samples)

population_counts_by_time <- HPS_pop_long %>% 
                                group_by(WEEK) %>% 
                                summarize(N=n()) %>% 
                                select(N) %>% 
                                unlist()

n_areas <- nlevels(HPS_pop_long$AREA)
n_weeks <- max(HPS_pop_long$WEEK)

direct_estimates <- list()
for(i in 1:Nsim){
  print(paste0("Sample #", i))
  HPS_sample_long <- samples[[i]]

  direct_estimates[[i]] <- get_direct_estimates(HPS_sample_long,
                                                population_counts_by_time,
                                                n_areas=n_areas, n_weeks=n_weeks)
}

if(case=="binary"){
  save(direct_estimates, file=paste0(save_dir, "direct_estimates.RData"))
}else{
  save(direct_estimates, file=paste0(save_dir, "direct_estimates_gauss.RData"))
}
