library(dplyr)
library(foreach)
library(doParallel)

save_dir <- "data/"
load(paste0(save_dir,"simulation_gtulm_results.RData"))
load(paste0(save_dir,"empirical_samples_gauss.RData"))

pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")
n_weeks <- max(HPS_pop_long$WEEK)
n_areas <- nlevels(HPS_pop_long$AREA)
popX <- model.matrix(as.formula(paste("~", pop_covars)), data=HPS_pop_long)
popPsi <- model.matrix(~AREA-1, data=HPS_pop_long)
popweeks <- HPS_pop_long$WEEK
popN <- nrow(HPS_pop_long)
Nsim <- nrow(gtulm_fit)

gibbs_dep_results <- c()

numCores <- 4
cl<-makeCluster(numCores, type="FORK", outfile="")
registerDoParallel(cl)

gtulm_post_preds <-foreach(i=1:Nsim, .combine="rbind", .verbose=TRUE,
                      .packages=c('dplyr')) %dopar% {
  print(paste0("Sim #", i))
  gibbs_dep_out <- gtulm_fit[i,]
  
  betas <- gibbs_dep_out$betas
  sigmas <- sqrt(gibbs_dep_out$sigma2)
  eta <- gibbs_dep_out$eta
  Xbeta <- popX%*%t(betas)
  eta_fs <- apply(eta, 1,
                  function(x) (popPsi%*%x)[cbind(seq_along(popweeks), popweeks)])
  pred.means <- as.matrix(Xbeta + eta_fs)
  rm(Xbeta, eta_fs)
  gc(verbose=F)
  print("step 1 done")
  
  Sigma <- matrix(rep(sigmas, popN), nrow=popN, byrow=T)
  posterior.pred <- rnorm(n=length(pred.means),
                                mean=c(pred.means),
                                sd=c(Sigma))
  posterior.pred <- matrix(posterior.pred, nrow=popN)
  print("step 2 done")

  #tabulate by week area
  gibbs_dep_agg <- data.frame(AREA=HPS_pop_long$AREA,
                              WEEK=HPS_pop_long$WEEK,
                              posterior.pred) %>%
                          group_by(AREA, WEEK) %>%
                          summarize_all(mean) %>%
                          arrange(WEEK, AREA)

  gibbs_dep_results <- data.frame(AREA=gibbs_dep_agg$AREA,
                                  WEEK=gibbs_dep_agg$WEEK,
                                  MEAN=rowMeans(gibbs_dep_agg[,-c(1,2)]),
                                  CI_LOWER=apply(gibbs_dep_agg[,-c(1,2)], 1, quantile,.025),
                                  CI_UPPER=apply(gibbs_dep_agg[,-c(1,2)], 1, quantile,.975))
  
  rm(Sigma, posterior.pred, pred.means, gibbs_dep_agg)
  gc(verbose=F)

  return(gibbs_dep_results)
}
stopCluster(cl)

save(gtulm_post_preds,
     file=paste0(save_dir,"simulation_gtulm_postpreds.RData"))

