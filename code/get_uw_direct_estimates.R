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

UW_direct_estimates <- list()
for(i in 1:Nsim){
  print(paste0("Sample #", i))
  HPS_sample_long <- samples[[i]]

  if(case=="binary"){
    UW_direct_estimates[[i]] <- HPS_sample_long %>%
        group_by(WEEK, AREA) %>%
        summarize(prop=mean(RESPONSE),
                  n=n(),
                  SE=sqrt(prop*(1-prop)/n)) %>%
        mutate(CI_LOWER = prop - 1.96*SE,
               CI_UPPER = prop + 1.96*SE)
  }else{
    UW_direct_estimates[[i]] <- HPS_sample_long %>%
        group_by(WEEK, AREA) %>%
        summarize(prop=mean(RESPONSE),
                  n=n(),
                  SE=sd(RESPONSE)/sqrt(n)) %>%
        mutate(CI_LOWER = prop - 1.96*SE,
               CI_UPPER = prop + 1.96*SE)
  }
}

if(case=="binary"){
  save(UW_direct_estimates, file=paste0(save_dir, "UW_direct_estimates_binary.RData"))
}else{
  save(UW_direct_estimates, file=paste0(save_dir, "UW_direct_estimates_gauss.RData"))
}
