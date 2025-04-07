library(Matrix)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
source("code/helper_functions.R")

save_dir <- "data/"

load(paste0(save_dir, "simulation_bbulm_results_binary.RData"))
load(paste0(save_dir, "simulation_btulm_results_binary.RData"))
load(paste0(save_dir, "empirical_samples.RData"))

Nsim <- length(samples)

#Parallelization parameters 
numCores <- 5
cl<-makeCluster(numCores, type="FORK", outfile="") #
registerDoParallel(cl)

n_weeks <- max(HPS_pop_long$WEEK)
n_areas <- nlevels(HPS_pop_long$AREA)

pop_covars <- select(HPS_pop_long, starts_with("COVAR")) %>% names %>% paste(collapse="+")
week_covars <-  gsub("\\+COVAR.NEW", "", pop_covars)

pop_grouped <- HPS_pop_long %>%
                group_by(WEEK, AREA, COVAR.1, COVAR.2,
                         COVAR.3, COVAR.4, COVAR.NEW) %>%
                         summarize(popsize=n()) %>%
               ungroup()

popPsi <- Matrix(model.matrix(~factor(pop_grouped$AREA) -1))
popPsiTime <- Matrix(model.matrix(~factor(pop_grouped$AREA):factor(pop_grouped$WEEK) -1))
popweeks <- as.integer(pop_grouped$WEEK)
pop_covars <- select(pop_grouped, starts_with("COVAR")) %>% names %>% paste(collapse="+")
week_covars <-  gsub("\\+COVAR.NEW", "", pop_covars)
popX <- model.matrix(as.formula(paste("~", pop_covars)), data=pop_grouped)
popX <- Matrix(data=as.matrix(popX))
popN <- nrow(popX)

week_areas <- pop_grouped %>%
                ungroup %>%
                select(WEEK,AREA)%>%
                unique 
    
n_cells <- nrow(week_areas)

truth <- HPS_pop_long %>% 
           group_by(AREA, WEEK) %>% 
           summarize(P=mean(RESPONSE)) %>%
           arrange(WEEK, AREA)

postpreds <-foreach(sim_num=1:Nsim, .combine="rbind", .verbose=TRUE,
                      .packages=c('dplyr', 'microbenchmark','Matrix')) %dopar% {

  print(paste0("Processing simulation #", sim_num))
  #for BTULM
  print("Start BTULM")
  gibbs_dep_out <- btulm_fit[sim_num,]
  gibbs_indep_out <- bbulm_fit[sim_num,]

  n_iter <- length(gibbs_dep_out$phi)

  #posterior preds
  betas <- gibbs_dep_out$betas
  eta <- gibbs_dep_out$eta
  Xbeta <- popX%*%t(betas)
  eta <- matrix(eta,nrow=(n_weeks*n_areas), byrow=T) 
  Psi_eta <- popPsiTime %*% eta
  pred.probs <- plogis(as.matrix(Xbeta + Psi_eta))

  ##Also save posterior preds for a few areas
  ## as well as Xbeta + Psieta for a few responses
  linpred_traceplots <- as.matrix(Xbeta + Psi_eta)[c(12,24,100),]

  rm(Xbeta, Psi_eta, betas)
  gc(verbose=F)

  posterior.pred <- rbinom(n=length(pred.probs),
                           size=pop_grouped$popsize,
                           prob=c(pred.probs))
  rm(pred.probs)
  gc(verbose=F)
  posterior.pred <- matrix(posterior.pred, nrow=popN)
  
  #tabulate by week area
  gibbs_dep_agg <- data.frame(AREA=pop_grouped$AREA, 
                              WEEK=pop_grouped$WEEK, 
                              posterior.pred) %>%
                          group_by(AREA, WEEK) %>%
                          summarize_all(sum) %>%
                          arrange(WEEK, AREA)

  pop_totals <- pop_grouped %>% group_by(WEEK, AREA) %>%  summarize(N=sum(popsize))
  gibbs_dep_agg[,-c(1,2)] <- gibbs_dep_agg[,-c(1,2)]/pop_totals$N

  #these correspond to MA Week 1, Minn. Week 3, Connecticut week 10
  area_pred_traceplots <- gibbs_dep_agg[c(20,120,447),-c(1,2)]  %>% as.matrix(nrow=3)
  gibbs_dep_traceplots <- rbind(linpred_traceplots, area_pred_traceplots)
  rownames(gibbs_dep_traceplots) <- c("individual_12", "individual_24", "individual_100",
                                       "MA_Week1", "MN_Week3", "CT_Week10")
  rm(posterior.pred)
  gc(verbose=F)

  gibbs_dep_results <- data.frame(AREA=gibbs_dep_agg$AREA, 
                                  WEEK=gibbs_dep_agg$WEEK, 
                                  MEAN=rowMeans(gibbs_dep_agg[,-c(1,2)]), 
                                  CI_LOWER=apply(gibbs_dep_agg[,-c(1,2)], 1, quantile,.025),
                                  CI_UPPER=apply(gibbs_dep_agg[,-c(1,2)], 1, quantile,.975))
  print("Start independent")
  #get BBULM post-pred
  gibbs_indep_agg <- c()
  for(tt in 1:n_weeks){
     print(paste0("Week ", tt))
     mod <- gibbs_indep_out[[tt]]

     pop_grouped_week <- filter(pop_grouped, WEEK==tt)
     n_pop_week <- nrow(pop_grouped_week)
     predX_week <- model.matrix(as.formula(paste("~", week_covars)), data=pop_grouped_week)
     predPsi_week <-  model.matrix(~AREA-1, data=pop_grouped_week)
     pred.probs <- plogis(cbind(predX_week, predPsi_week) %*% t(cbind(mod$Beta, mod$Eta)))

     preds.tmp <- rbinom(n=length(pred.probs),
                         size=pop_grouped_week$popsize,
                         prob=pred.probs)

     preds <- matrix(preds.tmp, nrow=n_pop_week)
   
     tmp <- data.frame(AREA = pop_grouped_week$AREA,
                       WEEK = tt,
                       preds) %>%
             group_by(AREA,WEEK) %>%
             summarize_all(sum) %>% 
             arrange(WEEK,AREA)

     pop_totals_week <- filter(pop_totals, WEEK==tt)$N

     tmp[,-c(1,2)] <- tmp[,-c(1,2)]/pop_totals_week
     gibbs_indep_agg <- rbind(gibbs_indep_agg, tmp) 
   }
   gibbs_indep_agg <- ungroup(gibbs_indep_agg)
   gibbs_indep_traceplots <- cbind(
                          gibbs_indep_agg %>% filter(WEEK==6, AREA=="Georgia") %>% select(X3),
                          gibbs_indep_agg %>% filter(WEEK==4, AREA=="Louisiana") %>% select(X2),
                          gibbs_indep_agg %>% filter(WEEK==11, AREA=="Virginia") %>% select(X2)
                       )

  gibbs_indep_results <- data.frame(AREA=gibbs_indep_agg$AREA, 
                                    WEEK=gibbs_indep_agg$WEEK, 
                                    MEAN=rowMeans(gibbs_indep_agg[,-c(1,2)]), 
                                    CI_LOWER=apply(gibbs_indep_agg[,-c(1,2)], 1, quantile,.025),
                                    CI_UPPER=apply(gibbs_indep_agg[,-c(1,2)], 1, quantile,.975))

  
  return(list(gibbs_dep_post_pred=gibbs_dep_results, 
              gibbs_dep_traceplots=gibbs_dep_traceplots,
              gibbs_indep_post_pred=gibbs_indep_results,
              gibbs_indep_traceplots=gibbs_indep_traceplots
              )
         )
}
stopCluster(cl)

save(postpreds, truth,
     file=paste0(save_dir,"binary_simulation_results_postpreds.RData"))
