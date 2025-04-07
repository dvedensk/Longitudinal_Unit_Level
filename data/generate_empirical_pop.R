library(readr)
library(purrr)
library(tidyr)
library(dplyr)

discrete <- F

#Data can be downloaded with a command like e.g.
# for i in {1..12}
# do if [ $i -lt 10 ]
# then wget https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk$i/HPS_Week0$i'_PUF_CSV.zip'
# else wget https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk$i/HPS_Week$i'_PUF_CSV.zip'
# fi
# done
# inconsistent about leading zeros, so handle 10-12 separately.
# unzip HPS*zip
# then run the code below

#following https://stackoverflow.com/a/65001063
filenames <- list.files(pattern="pulse2020_puf_[0-9]+.csv$",recursive=T)
filenames <- filenames[1:12]
combined_data <- purrr::map_df(filenames, ~read_csv(.x) %>% mutate(filename = .x))
                                
n_weeks <- length(filenames) #12

if(discrete){
  HPS_df <- select(combined_data, SCRAM, WEEK, AREA=EST_ST, TBIRTH_YEAR, 
                   ORIG_WEIGHT=PWEIGHT, EGENDER, RESPONSE=EXPCTLOSS, RACE=RRACE)
}else{ 
  HPS_df <- select(combined_data, SCRAM, WEEK, AREA=EST_ST, TBIRTH_YEAR, 
                   ORIG_WEIGHT=PWEIGHT, EGENDER, RESPONSE=TSPNDFOOD, RACE=RRACE)
}


#pare things down:
HPS_df <- filter(HPS_df, RESPONSE>0) %>% 
              mutate(COVAR.1=factor(EGENDER, levels=c(1,2),
                                             labels=c("MALE","FEMALE"))) %>%
              mutate(COVAR.2=as.integer(2020-TBIRTH_YEAR)) %>% 
              mutate(COVAR.3=factor(RACE, levels=c(1,2,3,4),
                                          labels=c("WHITE", "BLACK", "ASIAN", "OTHER"))) %>% 
              select(-TBIRTH_YEAR, -EGENDER, -RACE) 

if(discrete){
  #1) Yes, 2) No, so map 2 to 0
  HPS_df <- HPS_df %>% mutate(RESPONSE=as.integer(ifelse(RESPONSE==2, 0, 1)))
}else{
  boxcox.response <- forecast::BoxCox(HPS_df$RESPONSE, lambda="auto")
  HPS_df$RESPONSE <- boxcox.response
}

#use list from provided excel sheet
HPS_df$AREA <- recode_factor(HPS_df$AREA,
                    '01'='Alabama', '02'='Alaska', '04'='Arizona',
                    '05'='Arkansas','06'='California','08'='Colorado',
                    '09'='Connecticut', '10'='Delaware', '11'='District of Columbia',
                    '12'='Florida','13'='Georgia','15'='Hawaii',
                    '16'='Idaho','17'='Illinois','18'='Indiana',
                    '19'='Iowa','20'='Kansas','21'='Kentucky',
                    '22'='Louisiana','23'='Maine','24'='Maryland',
                    '25'='Massachusetts','26'='Michigan','27'='Minnesota',
                    '28'='Mississippi','29'='Missouri','30'='Montana',
                    '31'='Nebraska','32'='Nevada','33'='New Hampshire',
                    '34'='New Jersey','35'='New Mexico','36'='New York',
                    '37'='North Carolina','38'='North Dakota','39'='Ohio', 
                    '40'='Oklahoma','41'='Oregon','42'='Pennsylvania',
                    '44'='Rhode Island','45'='South Carolina','46'='South Dakota',
                    '47'='Tennessee','48'='Texas','49'='Utah',
                    '50'='Vermont','51'='Virginia','53'='Washington',
                    '54'='West Virginia','55'='Wisconsin','56'='Wyoming')

HPS_df_long <- HPS_df %>% arrange(WEEK) %>% 
                          group_by(SCRAM) %>% 
                          mutate(RESPONSE_NUMBER=1:n()) %>% 
                          #Take only the first WEIGHT, which we'll 
                          #use to make a SIZE_VAR later:
                          mutate(ORIG_WEIGHT=ORIG_WEIGHT[1]) %>%
                          arrange(SCRAM) %>% 
                          ungroup()

#Want to add a logical flag for whether a response is a follow-up
#and a PREV_RESPONSE column for response at time t-1 if follow-up (o.w. NA)
HPS_df_long <- HPS_df_long %>% 
                 mutate(IS_FOLLOWUP=ifelse(RESPONSE_NUMBER>1, TRUE, FALSE)) %>%
                 mutate(PREV_RESPONSE=ifelse(IS_FOLLOWUP, lag(RESPONSE), NA))

#Also want a wide data frame to sample from. 
HPS_df_wide <- HPS_df_long %>%  pivot_wider(values_from=c(RESPONSE, WEEK), 
                                            names_from=c(RESPONSE_NUMBER),
                                            id_cols=SCRAM)

#Join back to get covariate values
HPS_df_wide <- HPS_df_wide %>% left_join(HPS_df_long, 
                                         by=c("SCRAM"="SCRAM", 
                                              "RESPONSE_1"="RESPONSE", 
                                              "WEEK_1"="WEEK")) %>%
                               select(-RESPONSE_NUMBER)

if(!discrete){
  save(HPS_df_long, HPS_df_wide, file="data/HPS_empirical_pop_gauss_df.RData")
}else{
  save(HPS_df_long, HPS_df_wide, file="data/HPS_empirical_pop_binary_df.RData")
}

