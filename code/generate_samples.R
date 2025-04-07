library(dplyr)
library(purrr)
library(tidyr)
library(Matrix)
library(sampling)
source("code/helper_functions.R")

##########################
### Simulation Settings ##
##########################
set.seed(555)
case <- "gaussian"

Nsim <- 100

save_dir <- "data/"

if(case=="binary"){
  load(paste0(save_dir, "HPS_empirical_pop_binary_df.RData"))
}else{
  load(paste0(save_dir, "HPS_empirical_pop_gauss_df.RData"))
}
HPS_pop_long <- HPS_df_long
HPS_pop_wide <- HPS_df_wide

#add quadratic term for age
HPS_pop_long <- HPS_pop_long %>% mutate(COVAR.4=COVAR.2^2)
HPS_pop_wide <- HPS_pop_wide %>% mutate(COVAR.4=COVAR.2^2)

state_exclude <- c("Alaska", "Hawaii")

HPS_pop_long <- filter(HPS_pop_long, !AREA %in% state_exclude) %>% 
                     mutate(AREA=droplevels(AREA))
HPS_pop_wide <- filter(HPS_pop_wide, !AREA %in% state_exclude) %>% 
                     mutate(AREA=droplevels(AREA))

n_weeks <- max(HPS_pop_long$WEEK)
n_areas <- nlevels(HPS_pop_long$AREA)

pop_size <- nrow(HPS_pop_wide)
sample_size <- floor(.02 * pop_size)

if(case == "binary"){
  new_covar <- addNA(factor(HPS_pop_long$PREV_RESPONSE))
  levels(new_covar) <- c("PREV_EQ_NO", "PREV_RESP_YES", "PREV_RESP_NA")
  new_covar <- relevel(new_covar, ref="PREV_RESP_NA")
  HPS_pop_long$COVAR.NEW <- new_covar
}

samples <- list()
params <- list()
for(i in 1:Nsim){
  set.seed(i + 15) 

  print(paste("Taking sample...", i))
  sample_out <- get_sample(HPS_pop_long=HPS_pop_long, HPS_pop_wide=HPS_pop_wide, 
                         sample_size=sample_size, n_weeks=n_weeks, n_areas=n_areas)
  HPS_sample_long <- sample_out$HPS_sample_long
  HPS_sample_long <- HPS_sample_long %>% group_by(WEEK) %>%
                                         mutate(SCALE_WEIGHT=
                                                PWEIGHT*n()/sum(PWEIGHT)) %>%
                                         ungroup() %>%
                                         arrange(WEEK)
  
  samples[[i]] <- HPS_sample_long
  params[[i]] <- sample_out$params
}

if(case=="binary"){
  save(samples, HPS_pop_long, file=paste0(save_dir, "empirical_samples.RData"))
}else{
  save(samples, params, HPS_pop_long, file=paste0(save_dir, "empirical_samples_gauss.RData"))
}

