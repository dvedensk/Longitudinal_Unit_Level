library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)
library(foreach)
library(doParallel)
library(microbenchmark)
library(mvtnorm)
library(ggplot2)
library(Matrix)
library(LaplacesDemon)
library(BayesLogit)
library(sampling)
source("code/helper_functions.R")
source("code/models/BTULM.R")

##########################
### Simulation Settings ##
##########################
set.seed(555)

gibbs.runs <- 8000
gibbs.burn <- 1000
discrete<-T

save_dir <- "data/"

#######################
### Set True Values ###
#######################

load(paste0(save_dir, "HPS_empirical_pop_binary_df.RData"))
HPS_pop_long <- HPS_df_long
HPS_pop_wide <- HPS_df_wide

HPS_pop_long <- HPS_pop_long %>% mutate(COVAR.4=COVAR.2^2)
HPS_pop_wide <- HPS_pop_wide %>% mutate(COVAR.4=COVAR.2^2)

state_exclude <- c("Alaska", "Hawaii")

HPS_pop_long <- filter(HPS_pop_long, !AREA %in% state_exclude) %>% 
                     mutate(AREA=droplevels(AREA))
HPS_pop_wide <- filter(HPS_pop_wide, !AREA %in% state_exclude) %>% 
                     mutate(AREA=droplevels(AREA))

HPS_pop_long <- HPS_pop_long %>% group_by(WEEK) %>%
                                        mutate(SCALE_WEIGHT=
                                               ORIG_WEIGHT*n()/sum(ORIG_WEIGHT)) %>%
                                        ungroup() %>%
                                        arrange(WEEK)

new_covar <- addNA(factor(HPS_pop_long$PREV_RESPONSE))
levels(new_covar) <- c("PREV_EQ_NO", "PREV_RESP_YES", "PREV_RESP_NA")
new_covar <- relevel(new_covar, ref="PREV_RESP_NA")
HPS_pop_long$COVAR.NEW <- new_covar

popPsi <- Matrix(model.matrix(~factor(HPS_pop_long$AREA) -1))
colnames(popPsi) <- c()
pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")
popX <- model.matrix(as.formula(paste("~", pop_covars)), data=HPS_pop_long)
popX <- Matrix(data=as.matrix(popX))
colnames(popX) <- c()

sigma2_beta <- 10^4

areas <- HPS_pop_long$AREA
weights <- HPS_pop_long$SCALE_WEIGHT
weeks <- HPS_pop_long$WEEK

model_data <- list(X=popX, Y=HPS_pop_long$RESPONSE, Psi=popPsi, 
                   sigma2_beta=sigma2_beta, weights=weights, 
                   weeks=weeks, areas=areas)

### Fit Gibbs
gibbs_dep_data <- model_data
gibbs_dep_data$n_gibbs <- gibbs.runs
gibbs_dep_data$n_burn <- gibbs.burn

gibbs_dep_out <- do.call(fit_dep_gibbs, gibbs_dep_data)
print(paste0("Ran in ", gibbs_dep_out$runtime, " seconds"))
print(paste0(gibbs.runs/gibbs_dep_out$runtime, " iterations per second"))
print(paste0(gibbs_dep_out$runtime/gibbs.runs, " seconds per iteration"))

save(HPS_pop_long, HPS_pop_wide, gibbs_dep_out,
     file=paste0(save_dir,"binary_data_analysis_results.RData"))
