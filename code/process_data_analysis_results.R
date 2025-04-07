library(Matrix)
library(dplyr)
library(ggplot2)
library(stringr)
library(survey)
library(readr)

save_dir <- "output/figures/"
data_dir <- "data/"

make_areal_plot <- function(plot.df, model_name, sex_sub, race_sub, age_sub, age_str){
 plot.df <- filter(plot.df, model==model_name, SEX==sex_sub, RACE==race_sub, AGE_CAT==age_sub)
 states <- tolower(levels(plot.df$AREA))
 state_poly <- map_data("state") %>% filter(region %in%  states)
 plot.df$region<- tolower(plot.df$AREA)
 plot.df <- merge(plot.df, state_poly, by="region")

 lower_lim <- .95*min(plot.df$prop)
 upper_lim <- 1.05*max(plot.df$prop)

 plot.df %>%
   ggplot(mapping = aes(x=long, y=lat,
                        fill=prop, group=group)) +
     geom_polygon(color="black", size=0.15) +
     scale_fill_viridis_c(name="Estimated proportion",
                          limits=c(lower_lim, upper_lim), oob=scales::squish) +
     coord_map() +
     facet_wrap(~WEEK) +
            ggthemes::theme_map() +
            theme(legend.position="right",
                  plot.margin=grid::unit(c(0,0,0,0), "mm"))

 filepath <- paste0(save_dir,
                    model_name,"_",
                    "estimated_proportion_",
                    tolower(sex_sub), "_", 
                    tolower(race_sub), "_", 
                    gsub(" ", "_", age_str),
                    ".pdf")

 ggsave(filepath, dpi=300)
     system2(command = "pdfcrop",
              args    = c(filepath,
                          filepath))

}

load(paste0(data_dir, "binary_data_analysis_results.RData"))
census.raw <- readr::read_csv(paste0(data_dir,"census_tables.csv"))

n_weeks <- 12
#age_breaks <- c(18,25,40,55,65, 100)

census.table <- census.raw %>%
                   select(STATE, NAME, SEX, ORIGIN, RACE, AGE,
		         POPESTIMATE2021) %>% #does choice of estimate matter much?
                  filter(NAME %in% levels(HPS_pop_long$AREA)) %>%
	          filter(SEX!=0, AGE >= 18, ORIGIN!=0) %>%
                          mutate(SEX=recode(as.factor(SEX), "1" = "MALE", "2" = "FEMALE")) %>%
                  mutate(RACE=recode(as.factor(RACE), "1" = "White", "2" = "Black", "4" = "Asian",
                                                      "3" = "Other", "5" = "Other", "6" = "Other")) %>%
		  group_by(NAME, SEX, AGE, RACE) %>%
		  summarize(COUNT=sum(POPESTIMATE2021)) %>%
                  mutate(AGE.2=AGE^2) %>% 
                  mutate(AREA=as.factor(NAME)) %>%
                  ungroup() %>% 
                  select(AREA, SEX, AGE, RACE, AGE.2, COUNT)

cell_pops <- census.table$COUNT #keep a count of how many people are in each cell to start
n_cells <- nrow(census.table)

prev_resp_fact <- c("PREV_RESP_NA", "PREV_EQ_NO", "PREV_RESP_YES")

census.table <- merge(census.table, 
                      factor(prev_resp_fact, levels=prev_resp_fact)) %>%
                rename(COVAR.NEW=y) %>%
                mutate(COUNT=ifelse(COVAR.NEW != "PREV_RESP_NA", 0, COUNT))

betas <-gibbs_dep_out$betas[,]
etas <- gibbs_dep_out$eta[,,]
rm(gibbs_dep_out)

popPsi <- Matrix(model.matrix(~factor(census.table$AREA) -1))
names(popPsi) <- c()
#
#Order should be Female, Age, Black, Asian, Other, Age^2, COVAR.NEW
popX <- Matrix(model.matrix(~SEX+AGE+RACE+AGE.2+COVAR.NEW, 
                            data=census.table))

#need to swap RACEASIAN and RACEOTHER to match betas from model
popX <- popX[,c(1:4,6,5,7:ncol(popX))]
n_gibbs <- nrow(betas)
age_breaks <- c(17, seq(25,65,5),100)

#This loop generates the aggregated posterior predictions
for(tt in 1:n_weeks){
  print(paste0("tt = ", tt))
  etas_t <- etas[,,tt]
  beta_eta <- t(cbind(betas, etas_t)) #each col an MCMC draw

  if(tt == 1){
    #only PREV_RESP_NA is relevant for week 1, so take this subset
    idx <- 1:n_cells #week 1 has different indices
    #every MCMC iteration has the same population just in this step
    counts_t <- replicate(n_gibbs, cell_pops) 
  }else{
    id_no <- which(census.table$COVAR.NEW=="PREV_EQ_NO")
    id_yes <- which(census.table$COVAR.NEW=="PREV_RESP_YES")
    idx <- c(id_no, id_yes)
  }

  popX_t <- popX[idx,]
  popPsi_t <- popPsi[idx,]
  #update eta
  linpred <- cbind(popX_t, popPsi_t) %*% beta_eta #get Xbeta + Psi*eta
  pred.probs <- plogis(as.matrix(linpred))
  
  posterior.pred <- rbinom(n=length(pred.probs), 
                           size=c(counts_t),
                           prob=c(pred.probs))
  posterior.pred <- matrix(posterior.pred, ncol=n_gibbs)

  aggregate_yes <- data.frame(AREA=census.table[idx,]$AREA, 
		              SEX=census.table[idx,]$SEX, 
		              AGE=census.table[idx,]$AGE, 
		              RACE=census.table[idx,]$RACE, 
                              MCMC=posterior.pred) %>% 
                   mutate(AGE_CAT=cut(AGE,   
                                      breaks=age_breaks))  %>%
                   select(AREA, SEX, AGE, AGE_CAT, everything()) %>%
                   group_by(AREA, SEX, RACE, AGE_CAT) %>% 
                   summarize_all(sum) %>%
                   select(-AGE)

  ##aggregate by gender and age-group
  aggregate_counts <- data.frame(AREA=census.table[idx,]$AREA, 
                                 SEX=census.table[idx,]$SEX, 
                                 AGE=census.table[idx,]$AGE, 
                                 RACE=census.table[idx,]$RACE, 
                                 MCMC=counts_t) %>% 
                   mutate(AGE_CAT=cut(AGE,   
                                      breaks=age_breaks))  %>%
                   select(AREA, SEX, AGE, AGE_CAT, everything()) %>%
                   group_by(AREA, SEX, RACE, AGE_CAT) %>% 
                   summarize_all(sum) %>%
                   select(-AGE)

  props <- aggregate_yes[,-c(1,2,3,4)]/aggregate_counts[,-c(1,2,3,4)] 
  total <- aggregate_yes[, -c(1:4)]
  agg_results <- data.frame(AREA=aggregate_yes$AREA,
                            SEX=aggregate_yes$SEX, 
                            RACE=aggregate_yes$RACE, 
                            AGE_CAT=aggregate_yes$AGE_CAT, 
                            WEEK=tt,
                     #       MCMC=props) %>% group_by(AREA)
                            MCMC=total) %>% group_by(AREA)
  #divide total yes's by (yes+no) for state totals

  ##UPDATE COUNTS
  if(tt == 1){
  #posterior pred is same as matrix of yes counts
  #counts_t - posterior.pred is same as a matrix of no counts
  #paste these together to get a popsize matrix for the next
  #time step, which will have different pop sizes for each MCMC step
    yes_count <- posterior.pred
  }else{
    #append as we go, except in week 1 when we don't have prev results
    agg_results <- rbind(prev_week_agg, agg_results) 
    #at time tt >=2 yes counts are distributed between the two halves of the matrix
    #need to add the first half of the matrix to the second half to get a matrix of yes counts
    #then subtract this from cell_pops??
    yes_count <- posterior.pred[1:n_cells,] + posterior.pred[(n_cells+1):(2*n_cells),]
  }
  no_count <- cell_pops - yes_count
  counts_t <- rbind(no_count,
                    yes_count)
  prev_week_agg <- agg_results
}

save(agg_results, 
     popX, 
     popPsi,
     etas,
     betas,
     file=paste0(data_dir,"binary_data_analysis_postpreds.RData"))

####Generate direct estimates with survey library and replicate weights
load(paste0(data_dir,"HPS_empirical_pop_binary_df.RData"))
dir_est.df <- HPS_df_long %>% 
                   filter(AREA!="Alaska", AREA!="Hawaii") %>%
                   mutate(AREA=droplevels(AREA)) %>%
                   mutate(AGE_CAT=cut(COVAR.2,
                                      breaks=age_breaks)) %>%
                   rename(SEX=COVAR.1, AGE=COVAR.2, RACE=COVAR.3) %>%
                   select(SCRAM, WEEK, AREA, ORIG_WEIGHT, 
                          RESPONSE, SEX, AGE, AGE_CAT, RACE) %>% 
                   unique %>%
                   left_join( select(census.table, -AGE.2),
                               by = c("AREA"="AREA",
                                      "SEX"="SEX",
                                      "RACE"="RACE",
                                      "AGE"="AGE")) %>%
                   group_by(WEEK) %>% 
                   mutate(INCLUSION_PROB=1/ORIG_WEIGHT) %>%
                   ungroup() %>% 
                   select(-COUNT,-COVAR.NEW, -AGE) %>%
                   unique()
                          
dir_est_results <- c()
for(week_num in 1:12){
  print(paste0("Dir Est week ", week_num))
  week_str <- str_pad(week_num,2,pad=0)
  repwgt <- read_csv(paste0("data/HPS_PUFs/HPS_Week", week_str,
                            "_PUF_CSV/pulse2020_repwgt_puf_", week_str, ".csv"))

  svy.tst <- left_join(filter(dir_est.df, WEEK==week_num), repwgt, by="SCRAM")
  hps.svy = svrepdesign(data=svy.tst,
                         weights=~ORIG_WEIGHT,
                         repweights = "PWEIGHT[1-9]+",
                         type="Fay",
                         rho=0.5)

  tmp <- svyby(~RESPONSE, ~AREA+SEX+RACE+AGE_CAT,
               design=hps.svy, FUN=svymean, na.rm=TRUE)
  rownames(tmp) <- c()
  tmp$WEEK <- week_num
  dir_est_results <- rbind(dir_est_results, tmp)  
}
dir_est_results <- tibble(dir_est_results) %>%
                   select(AREA, WEEK, SEX, RACE, AGE_CAT,
                          prop=RESPONSE, std_err=se) %>%
                   mutate(std_err=ifelse(std_err==0, NA, std_err))

agg_mcmc <- select(ungroup(agg_results), starts_with("MCMC")) %>% select(-MCMC.AGE_CAT)

#plot results
agg.df<- data.frame(AREA=agg_results$AREA, 
                         WEEK=agg_results$WEEK, 
                         SEX=agg_results$SEX,
                         RACE=agg_results$RACE,
                         AGE_CAT=agg_results$AGE_CAT, 
                         prop=rowMeans(agg_mcmc),
                         std_err=apply(agg_mcmc, 1, sd))

agg.df$RACE <- as.factor(toupper(agg.df$RACE))

backup <- dir_est_results
full.df <- rbind(cbind(agg.df, model="TULM"), 
                (cbind(dir_est_results, model="DirEst"))) 

 dir_est.df %>%
    group_by(WEEK, AREA, SEX, RACE, AGE_CAT) %>%
    tally %>%
    right_join(full.df, by=c("AREA", "SEX" , "RACE", "AGE_CAT", "WEEK")) %>%
    select(-std_err) %>% 
    tidyr::pivot_wider(names_from=model, values_from=prop) %>%
    filter(!is.na(TULM/DirEst))%>%
    ggplot(aes(x=n, y=DirEst/TULM)) +
    geom_point() +
    geom_line(aes(y=1), col="red") + 
    xlab("Sample size") +
    ylab("Ratio") +
    theme_classic()
ggsave("sample_size_comparison.pdf")

##
sex_sub <- "MALE"
sex_str <- tolower(sex_sub)
race_sub <- "WHITE"
race_str <- "White"
age_sub <- "(17,25]"
age_str <- gsub("\\(|\\[|\\]","",age_sub) %>%
             gsub(","," to ", .)

plot.subset.df <- filter(agg.df, SEX==sex_sub,
                                 RACE==race_sub,  
                                 AGE_CAT==age_sub) %>%
#                  select(-prop) %>%
                  rename(TULM_std_err=std_err) %>%
                  rename(TULM_prop=prop) %>%
                  arrange(AREA, WEEK)

DirEst_std_err <- filter(dir_est_results, SEX==sex_sub,
                                        RACE==race_sub,
                                        AGE_CAT==age_sub)  %>%
                  arrange(AREA, WEEK)

plot.subset.df$DirEst_std_err <- DirEst_std_err$std_err
plot.subset.df$DirEst_prop <- DirEst_std_err$prop

plot.subset.df <- mutate(plot.subset.df, ratio=TULM_std_err/DirEst_std_err)

#plots a faceted time series of estimated proportions for TULM and DE
make_areal_plot(full.df, "TULM", sex_sub, race_sub, age_sub, age_str)
make_areal_plot(full.df, "DirEst", sex_sub, race_sub, age_sub, age_str)

##Plots ratio of TULM to DE standard errors
##Figure 7 in the manuscript
states <- tolower(levels(plot.subset.df$AREA))
state_poly <- map_data("state") %>% filter(region %in%  states)
plot.subset.df$region<- tolower(plot.subset.df$AREA)
plot.subset.df <- merge(plot.subset.df, state_poly, by="region")

lower_lim <- .95*min(plot.subset.df$ratio)
upper_lim <- 1.05*max(plot.subset.df$ratio)
        
plot.subset.df %>% filter(ratio<.6 | is.na(ratio)) %>%
  ggplot(mapping = aes(x=long, y=lat,
                       fill=ratio, group=group)) +
    geom_polygon(color="black", size=0.15) +
    scale_fill_viridis_c(name="(TULM SE)/(DirEst SE)",
                         limits=c(lower_lim, upper_lim), oob=scales::squish) +
    labs(fill="Estimate") +
    coord_map() +
    facet_wrap(~WEEK)+
           ggthemes::theme_map() +
           theme(legend.position="right",
                 plot.margin=grid::unit(c(0,0,0,0), "mm"))

filepath <- paste0(save_dir,
                     tolower(race_sub), "_" , tolower(sex_sub), "_", gsub(" ", "_", age_str),
                "_SE_ratio.pdf")
  ggsave(filepath)

  system2(command = "pdfcrop",
          args    = c(filepath,
                      filepath))

'%!in%' <- function(x,y)!('%in%'(x,y))
library(covdata)
library(ggrepel)
library(tidycensus)

total_population_20 <- get_decennial(
  geography = "state", 
  variables = "P1_001N",
  year = 2020
)
total_pop <- total_population_20 %>% mutate(pop_100k = value/100000)
total_pop <- select(total_pop, state=NAME, pop_100k)

us_cov_deaths <- nytcovstate %>%  filter(date >= "2020-04-08" & date <= "2020-07-21")

#undo cumulative
us_cov_deaths <- us_cov_deaths %>% 
  arrange(state, date) %>% 
  group_by(state) %>% 
  mutate(daily_deaths=c(deaths[1],diff(deaths))) %>%
  mutate(daily_cases=c(cases[1],diff(cases))) %>%
  filter(date != "2020-04-08") #date not part of HPS window, just using it for diff

us_cov_deaths <- filter(us_cov_deaths, !is.na(state))
us_cov_deaths <- us_cov_deaths %>% left_join(total_pop, 
                                             by=join_by(state==state))

#end dates for collection
HPS_phases <- data.frame(week=1:12,
                         date=c(as.Date("2020-04-23"),
                                as.Date("2020-05-07"),
                                as.Date("2020-05-14"),
                                as.Date("2020-05-21"),
                                as.Date("2020-05-28"),
                                as.Date("2020-06-04"),
                                as.Date("2020-06-11"),
                                as.Date("2020-06-18"),
                                as.Date("2020-06-25"),
                                as.Date("2020-07-02"),
                                as.Date("2020-07-09"),
                                as.Date("2020-07-16")))

states_to_keep <- c("Texas", "Florida", "Georgia", "California","Indiana", "California", "Nevada")

HPS_phases <- plot.subset.df %>% 
                select(state=AREA, week=WEEK, TULM_prop) %>% 
                left_join(HPS_phases, by=join_by(week==week)) %>%
                filter(state %in% states_to_keep)

lower_lim <- .95*min(HPS_phases$TULM_prop)
upper_lim <- 1.05*max(HPS_phases$TULM_prop)

us_cov_deaths %>%
                filter(state %in% states_to_keep) %>%
   #filter(state %nin% c("New York", "New Jersey","Florida", "Connecticut",
   #                     "Kansas", "Alaska", "Hawaii")) %>% 
   group_by(state) %>%
   ggplot() +
     geom_line(aes(x=date,y=daily_cases/pop_100k)) +
     geom_vline(data=HPS_phases, aes(xintercept=date, alpha=.65, color=TULM_prop))+# , linetype="dashed") +
     scale_color_viridis_c(limits=c(lower_lim, upper_lim), oob=scales::squish) +
     geom_text(data = HPS_phases, aes(label = week, x = date, y = -Inf, alpha=.5), 
               angle = 0, inherit.aes = F, hjust = -.5, vjust = -.5, size=.25)+ 
     facet_wrap(~state) +
     guides(alpha='none') +
     labs(y="Cases per 100k", 
          x="Date", 
          color="Proportion expected job loss")  +
    theme_classic()

filepath <- paste0(save_dir,"cases_per_100k.pdf")
ggsave(filepath,dpi=300)

system2(command = "pdfcrop",
          args  = c(filepath,
                    filepath))




###Scatterplot for std err ratio
plot.subset.df %>%
    mutate(DirEst_std_err=ifelse(DirEst_std_err == 0, NA, DirEst_std_err)) %>%
    select(AREA, WEEK, TULM_std_err, DirEst_std_err, ratio) %>%
    group_by(WEEK) %>%
    mutate(col_group=ifelse(ratio>=1, "NA", "Est."),
           ratio=ifelse(is.na(ratio), 0, ratio)) %>%
    ggplot() +
      geom_point(aes(y=DirEst_std_err, x=TULM_std_err, color=col_group)) +
      geom_abline(slope=1, intercept=0) +
      ylim(c(0,.45)) +# xlim(c(0,1)) +
      facet_wrap(~WEEK) +
      scale_color_manual(values=c("black", "red"), guide="none") +
      theme_classic() +
      xlab("TULM SE") +
      ylab("Direct estimate SE")
ggsave(paste0(save_dir, "SE_ratio_scatterplot_classic.pdf"))




load(paste0(data_dir,"HPS_empirical_pop_binary_df.RData"))
HPS_df_long <- HPS_df_long %>% 
                  filter(AREA!="Alaska", AREA!="Hawaii") %>%
                  mutate(AREA=droplevels(AREA)) %>%
                  filter(COVAR.2 <= 85) %>% 
                  mutate(AGE_CAT=cut(COVAR.2,
                                      breaks=age_breaks)) %>%
                  rename(SEX=COVAR.1, AGE=COVAR.2, RACE=COVAR.3) 

tmp <- HPS_df_long %>% 
         select(SCRAM, WEEK, AREA, ORIG_WEIGHT, 
                RESPONSE, SEX, AGE, AGE_CAT, RACE) %>% 
         unique %>%
         left_join( select(census.table, -AGE.2),
                      by = c("AREA"="AREA",
                             "SEX"="SEX",
                             "RACE"="RACE",
                             "AGE"="AGE")) %>%
         group_by(WEEK) %>% 
         mutate(INCLUSION_PROB=1/ORIG_WEIGHT) %>%
         ungroup()

tmp <- unique(select(tmp,-COUNT,-COVAR.NEW))

full_resps <- tmp %>% 
                group_by(AREA, WEEK, SEX, RACE, AGE_CAT) %>%
                  summarize(DIR_EST_MEAN=mase::horvitzThompson(y=RESPONSE, pi=1/ORIG_WEIGHT)$pop_mean,
                            DIR_EST_VAR=mase::horvitzThompson(y=RESPONSE, pi=1/ORIG_WEIGHT, 
                                                              var_est=T, 
                                                              var_method = "LinHTSRS")$pop_mean_var,
                  UW_RESP=mean(RESPONSE))

dir_est_table <- select(full_resps, -UW_RESP)
dir_est_table <- dir_est_table %>%
    rename(prop=DIR_EST_MEAN) %>%
    mutate(std_err=sqrt(DIR_EST_VAR)) %>%
    select(-DIR_EST_VAR)
