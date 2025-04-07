get_sample <- function(HPS_pop_long=HPS_pop_long, 
                       HPS_pop_wide=HPS_pop_wide,
                       sample_size=2500, n_weeks, n_areas){
  #add a size variable for weights 
  individual_means <- rowMeans(HPS_pop_wide[ , c("RESPONSE_1", "RESPONSE_2", "RESPONSE_3")], na.rm=T)
  HPS_pop_wide$INDIV_MEAN <- individual_means

  HPS_pop_wide$SIZE_VAR <- as.numeric(exp(.1*scale(HPS_pop_wide$ORIG_WEIGHT) +
                                      .2*scale(HPS_pop_wide$INDIV_MEAN)))

  ###These are *true* probabilities of selection
  inclusion_probs <- inclusionprobabilities(HPS_pop_wide$SIZE_VAR, sample_size)
  inclusion_probs <- inclusion_probs/sum(inclusion_probs) * sample_size
  
  HPS_pop_wide$PWEIGHT <- 1/inclusion_probs #These are true weights
  HPS_pop_wide$PROB_SELECTION <- inclusion_probs #These are true prob. selection
  indices <- UPpoisson(inclusion_probs)
  SCRAM_to_sample <- HPS_pop_wide$SCRAM[indices==1]
  
  HPS_sample_wide <- filter(HPS_pop_wide, SCRAM %in% SCRAM_to_sample)
  #Take the sample in long format
  #Need to join first because HPS_pop_wide now has the sampling weights we'll 
  #need later
  HPS_sample_long <- HPS_pop_long %>% 
                        filter(SCRAM %in% SCRAM_to_sample)

  reduced_wide <- select(HPS_sample_wide, SCRAM, PWEIGHT, PROB_SELECTION)
  HPS_sample_long <- merge(HPS_sample_long, reduced_wide, 
                           all=TRUE, by="SCRAM") %>% tibble

  y1_y2.df <- get_y1_y2(HPS_sample_long)
  ###Fit model
  y_1.df <- y1_y2.df$y_1.df
  y_2.df <- y1_y2.df$y_2.df
  N_1 <- y1_y2.df$N_1
  N_2 <- y1_y2.df$N_2
  X_1 <- y1_y2.df$X_1
  X_2 <- y1_y2.df$X_2
  
  areas_1 <- y_1.df$AREA
  areas_2 <- y_2.df$AREA
  weeks_1 <- as.numeric(y_1.df$WEEK)
  weeks_2 <- as.numeric(y_2.df$WEEK)
  weights_1 <- as.numeric(y_1.df$SCALE_WEIGHT)
  weights_2 <- as.numeric(y_2.df$SCALE_WEIGHT)

  #return just the matrices for responses
  y_1 <- y_1.df %>% select(RESPONSE) %>% as.matrix()
  y_2 <- y_2.df %>% select(PREV_RESPONSE, RESPONSE) %>% as.matrix()

  params <- list(N_1=N_1, N_2=N_2,
                 y_1=y_1, y_2=y_2,
                 X_1=X_1, X_2=X_2,
                 n_areas=n_areas, n_weeks=n_weeks,
                 areas_1=areas_1, areas_2=areas_2,
                 weeks_1=weeks_1, weeks_2=weeks_2,
                 weights_1=weights_1, weights_2=weights_2) 
  return(list(params=params, HPS_sample_wide=HPS_sample_wide,
              HPS_sample_long=HPS_sample_long))
}

get_y1_y2 <- function(HPS_sample_long){
  #need to scale weights so that at each week
  #they sum to that week's sample size
  HPS_sample_long <- HPS_sample_long %>% group_by(WEEK) %>% 
                                    mutate(SCALE_WEIGHT=
                                              PWEIGHT*n()/sum(PWEIGHT)) %>% 
                                   ungroup() %>%
                                   arrange(WEEK)

  ##### First, get first-time respondents
  y_1.df <- filter(HPS_sample_long, !IS_FOLLOWUP)
  
  #get the format we need. i.e. list of N_{t,a} to track indices
  N_1 <- y_1.df %>% arrange(WEEK) %>% 
                 group_by(WEEK) %>% 
                 summarize(N_1=n()) %>%
                 select(N_1) %>% ungroup()
  N_1 <- N_1$N_1
  
  covar_names <-  y_1.df %>% select(starts_with("COVAR")) %>% 
                             names() %>% 
                             paste(collapse="+") 
  covar_formula <- as.formula(paste("~", covar_names))
            
  X_1 <- model.matrix(covar_formula, data=y_1.df)

  ###Get repeat respondents in a N_2 x 2 matrix
  ##for each week from 2 to T, check if there was a response in the previous week, 
  ##which we store as entry [i,1]
  #first handle respondents who answered twice
  y_2.df <- filter(HPS_sample_long, IS_FOLLOWUP)
  
  N_2 <- y_2.df %>% arrange(WEEK) %>% 
                 group_by(WEEK) %>% 
                 summarize(N_2=n()) %>%
                 select(N_2)

  N_2 <- c(0, N_2$N_2) # pad with 0 so we can treat it as a list of dim n_weeks
                       # obviously there are no repeat respondents in week 1.
  
  X_2 <- model.matrix(covar_formula, data=y_2.df)

  return(list(y_1.df=y_1.df, N_1=N_1, X_1=X_1,
              y_2.df=y_2.df, N_2=N_2, X_2=X_2))
}

get_direct_estimates <- function(HPS_sample_long, population_counts_by_time,
                                 n_areas, n_weeks){

  HPS_sample_long <-  HPS_sample_long %>% mutate(WEEK=factor(WEEK, levels=1:12))
  y_1 <- filter(HPS_sample_long, !IS_FOLLOWUP)
  y_2 <- filter(HPS_sample_long, IS_FOLLOWUP)
  #Direct estimate does not use scaled weights PWEIGHT instead
  # (sum w_i*y_i)/N_{population,t} where N_population is # unique respondents at time t
  response_1 <- y_1 %>% select(AREA, WEEK, RESPONSE, PWEIGHT) %>% 
                               mutate(INCLUSION_PROB=1/PWEIGHT)
  response_2 <- y_2  %>% select(AREA, WEEK, RESPONSE, PWEIGHT) %>% 
                               mutate(INCLUSION_PROB=1/PWEIGHT)
  response_by_week_area <- bind_rows(response_1, response_2)

  response_by_week_area <- data.frame(WEEK=factor(1:12), 
                                      N=population_counts_by_time) %>% 
                               right_join(response_by_week_area, by="WEEK") 

  response_by_week_area<-   response_by_week_area %>%           
                                   group_by(WEEK) %>%
                                   mutate(SCALE_WEIGHT=
                                          PWEIGHT*N/sum(PWEIGHT)) %>%
                                   mutate(INCLUSION_PROB=1/SCALE_WEIGHT) %>%
                                   ungroup() 

  response_by_week_area <- response_by_week_area %>%
                             group_by(AREA, WEEK) 
  
  
  ##Add in NA for WEEK/AREA not represented in sample?
  tmp <- expand.grid(AREA=levels(HPS_sample_long$AREA),
                     WEEK=levels(HPS_sample_long$WEEK)) #cover all possibilities regardless of sample
  dir_est_table <- HPS_sample_long %>% group_by(AREA, WEEK) %>% summarize(N=n())
  dir_est_table <- left_join(tmp, dir_est_table, by=c("AREA", "WEEK"))
  
  K <- nrow(dir_est_table)
  for(k in 1:K){
    if(k %% 50 == 0){print(paste0("Direct estimate #",k,"/",K))}
    curr_area <- dir_est_table[[k,"AREA"]]
    curr_week <-  dir_est_table[[k,"WEEK"]]
    curr_resps <- filter(response_by_week_area, AREA==curr_area & WEEK==curr_week)
    #try LinHB, if that doesn't work try LinHH
    HT <- mase::horvitzThompson(y=curr_resps$RESPONSE, 
                                pi=curr_resps$INCLUSION_PROB, 
                                var_est=T, var_method = "LinHTSRS") 
    dir_est_table[k,"DIR_EST_MEAN"] = HT$pop_mean
    dir_est_table[k,"DIR_EST_VAR"] = HT$pop_mean_var
    dir_est_table[k,"CI_LOWER"] = HT$pop_mean - 1.96*sqrt(HT$pop_mean_var)
    dir_est_table[k,"CI_UPPER"] = HT$pop_mean + 1.96*sqrt(HT$pop_mean_var)
  }
  
  dir_est_table$WEEK <- as.double(dir_est_table$WEEK)
  return(dir_est_table)  
}
