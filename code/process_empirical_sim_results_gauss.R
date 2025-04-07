library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(xtable)

save_dir <- "data/"

load(paste0(save_dir, "simulation_gtulm_postpreds.RData"))
load(paste0(save_dir, "simulation_gbulm_results.RData"))
load(paste0(save_dir, "direct_estimates_gauss.RData"))
load(paste0(save_dir, "UW_direct_estimates_gauss.RData"))
load(paste0(save_dir, "empirical_samples_gauss.RData"))

case <- "gaussian"

table_dir <- paste0("output/tables/",case,"/")
plot_dir <- paste0("output/figures/",case,"/")

ind_model <- "G-BULM"
dep_model <- "G-TULM"

Nsim <- length(direct_estimates)
#test_num <- 20

truth <- HPS_pop_long %>%
           group_by(AREA, WEEK) %>%
           summarize(P=mean(RESPONSE)) %>% 
           arrange(WEEK, AREA) %>% 
           ungroup()

#make tables
gibbs_dep_full <- gtulm_post_preds

fullEsts<-matrix(gibbs_dep_full$MEAN, ncol=Nsim)
fullL<-matrix(gibbs_dep_full$CI_LOWER, ncol=Nsim)
fullU<-matrix(gibbs_dep_full$CI_UPPER, ncol=Nsim)


gibbs_indep_full <- gbulm_fit

nonARests<-matrix(gibbs_indep_full$MEAN, ncol=Nsim)
nonARU<-matrix(gibbs_indep_full$CI_UPPER, ncol=Nsim)
nonARL<-matrix(gibbs_indep_full$CI_LOWER, ncol=Nsim)

direct_est_full <- map_df(direct_estimates,
                          bind_rows, 
                          .id="sim_num")%>%
                  tibble() %>% 
    mutate(sim_num=as.integer(gsub("result\\.","",sim_num))) %>%
    arrange(sim_num, WEEK, AREA) %>%
    rename(MEAN=DIR_EST_MEAN)

dirEsts<-matrix(direct_est_full$MEAN, ncol=Nsim)
dirU<-matrix(direct_est_full$CI_UPPER, ncol=Nsim)
dirL<-matrix(direct_est_full$CI_LOWER, ncol=Nsim)

UW_direct_est_full <- map_df(UW_direct_estimates,
                          bind_rows, 
                          .id="sim_num")%>%
                  tibble() %>% 
    mutate(sim_num=as.integer(gsub("result\\.","",sim_num)))

missing_rows <- direct_est_full %>%
    anti_join(UW_direct_est_full, by=c("WEEK", "AREA","sim_num")) %>%
    select(sim_num, WEEK, AREA, MEAN, N, DIR_EST_VAR, CI_LOWER, CI_UPPER)
UW_direct_est_full <- rbind(UW_direct_est_full,
                            setNames(missing_rows, names(UW_direct_est_full)))

UW_direct_est_full <- UW_direct_est_full %>%
    arrange(sim_num, WEEK, AREA) %>%
    rename(MEAN=prop)

uwDirEsts<-matrix(UW_direct_est_full$MEAN, ncol=Nsim)
uwDirU<-matrix(UW_direct_est_full$CI_UPPER, ncol=Nsim)
uwDirL<-matrix(UW_direct_est_full$CI_LOWER, ncol=Nsim)


##OVERALL TABLE
MSE <- c(
  (mean((truth$P - uwDirEsts)^2, na.rm=T)),
  (mean((truth$P - dirEsts)^2, na.rm=T)),
  (mean((truth$P - nonARests)^2)),
  (mean((truth$P - fullEsts)^2))
)
# Coverage

COV <- c(
  mean((truth$P < uwDirU) & (truth$P > uwDirL), na.rm=T),
  mean((truth$P < dirU) & (truth$P > dirL), na.rm=T),
  mean((truth$P < nonARU) & (truth$P > nonARL)),
  mean((truth$P < fullU) & (truth$P > fullL))
)

BIAS <- c(
  mean(abs(rowMeans(uwDirEsts) - truth$P), na.rm=T),
  mean(abs(rowMeans(dirEsts) - truth$P), na.rm=T),
  mean(abs(rowMeans(nonARests) - truth$P)),
  mean(abs(rowMeans(fullEsts) - truth$P))
)

int_score <- function(U,L,Q){
  return(
    (U - L) + 2/0.05*(Q < L)*(L - Q) + 2/0.05*(Q > U)*(Q - U)
  )
}

INT_SCORE <- c(
  mean(int_score(uwDirU, uwDirL, truth$P), na.rm=T),
  mean(int_score(dirU, dirL, truth$P), na.rm=T),
  mean(int_score(nonARU, nonARL, truth$P)),
  mean(int_score(fullU, fullL, truth$P))
)

tabDF <- data.frame(Method = c("UW Direct", "Direct", "Independent", "Dependent"),
                    MSE=MSE,
                    `Abs Bias`=BIAS,
                    Coverage=COV,
                    `Interval Score`=INT_SCORE)

print(xtable(tabDF, digits=-3,
             caption=paste("Overall results for direct- and model-based estimates averaged over",
                            Nsim,"simulations in the", case, "case")), 
      include.rownames=F, math.style.exponents=TRUE, 
      sanitize.colnames.function=function(x)gsub("\\."," ",x),
      file=paste0(table_dir,"overall_table.tex"))

#MSE by week
dir_est_mse <- tibble(truth,  sqr_err=(truth$P-dirEsts)^2) %>% 
   group_by(WEEK) %>% 
   select(-AREA) %>%
   summarize_all(mean, na.rm=T) %>% 
   select(WEEK, `Direct estimate`=sqr_err)

nonAR_mse <- tibble(truth,  sqr_err=(truth$P-nonARests)^2) %>% 
   group_by(WEEK) %>% 
   select(-AREA) %>%
   summarize_all(mean) %>% 
   select(WEEK, `Independent`=sqr_err)

AR_mse <- tibble(truth,  sqr_err=(truth$P-fullEsts)^2) %>% 
   group_by(WEEK) %>% 
   select(-AREA) %>%
   summarize_all(mean) %>% 
   select(WEEK, `Dependent`=sqr_err)

all_mses <- merge(merge(dir_est_mse, nonAR_mse, by="WEEK"), AR_mse, by="WEEK")
overall_mse <- colMeans(all_mses)
#could also do
#Reduce(function(x,y) merge(x,y, by="WEEK", all=T), list(dir_est_mse, nonAR_mse, AR_mse))
colnames(all_mses) <- c("WEEK", "Direct Est.", ind_model, dep_model)

all_mses %>% 
   rename(Week=WEEK) %>%
   pivot_longer(cols=-Week, names_to="Estimator", values_to="MSE") %>% 
   #make sure dir est comes first
   mutate(Estimator = forcats::fct_relevel(Estimator, "Direct Est.")) %>% 
   ggplot(aes(x=Week, y=MSE, color=Estimator, shape=Estimator)) + 
   geom_line() + 
   geom_point(size=2.5) + 
   scale_x_discrete(limits=1:12) + 
   scale_y_continuous(labels= function(x) format(x, scientific = TRUE)) +
   scale_color_brewer(palette = "Dark2") +
   theme_classic() +
   theme(aspect.ratio=.5) 
#   ggtitle(paste0("Time series of estimator MSEs in ",
#                   case, " case"))
ggsave(paste0(plot_dir,"MSE_time_series_plot.pdf"), dpi=300)

#fix margins for including in poster
   system2(command = "pdfcrop", 
           args    = c(paste0(plot_dir,"MSE_time_series_plot.pdf"),
                        paste0(plot_dir,"MSE_time_series_plot.pdf")))

##Table
all_mses <- rbind(all_mses, overall_mse)
all_mses[nrow(all_mses),"WEEK"] = "Overall"
print(xtable(all_mses, digits=5,
             caption=paste("Mean square error for all models averaged over",Nsim,"simulations in", 
                           case, " case")),
      include.rownames=F, file=paste0(table_dir,"MSE_table.tex"))

#Bias by week
dir_est_mab <- tibble(truth,  bias=abs(truth$P-rowMeans(dirEsts))) %>% 
   group_by(WEEK) %>% 
   summarize_all(mean, na.rm=T) %>% 
   select(WEEK, `Direct estimate`=bias)

nonAR_mab <- tibble(truth,  bias=abs(truth$P-rowMeans(nonARests))) %>% 
   group_by(WEEK) %>% 
   summarize_all(mean) %>% 
   select(WEEK, `Independent`=bias)

AR_mab <- tibble(truth,  bias=abs(truth$P-rowMeans(fullEsts))) %>% 
   group_by(WEEK) %>% 
   summarize_all(mean) %>% 
   select(WEEK, `Dependent`=bias)

all_mabs <- merge(merge(dir_est_mab, nonAR_mab, by="WEEK"), AR_mab, by="WEEK")
overall_mab <- colMeans(all_mabs)
all_mabs <- rbind(all_mabs, overall_mab)
all_mabs[nrow(all_mabs),"WEEK"] = "Overall"
print(xtable(all_mabs, digits=5, 
            caption=paste("Mean absolute bias for all models averaged over",Nsim,"simulations in", case, "case" )), 
      include.rownames=F, file=paste0(table_dir,"MAB_table.tex"))


#CI coverage
Q <- truth$P
gibbs_dep_CIs <- cbind(gibbs_dep_full, Q)
gibbs_dep_CIs$COVERED <- NA

gibbs_indep_CIs <- cbind(gibbs_indep_full, Q)
gibbs_indep_CIs$COVERED <- NA

dir_est_CIs <- cbind(direct_est_full, Q)
dir_est_CIs$COVERED <- NA

M <- nrow(gibbs_dep_full)
for(i in 1:M){
  print(paste0("m=",i))
  gibbs_dep_CIs$COVERED[i] <- between(gibbs_dep_CIs$Q[i], 
                                      gibbs_dep_CIs$CI_LOWER[i],
                                      gibbs_dep_CIs$CI_UPPER[i])

  gibbs_indep_CIs$COVERED[i] <- between(gibbs_indep_CIs$Q[i], 
                                        gibbs_indep_CIs$CI_LOWER[i],
                                        gibbs_indep_CIs$CI_UPPER[i])

  dir_est_CIs$COVERED[i] <- between(dir_est_CIs$Q[i], 
                                    dir_est_CIs$CI_LOWER[i],
                                    dir_est_CIs$CI_UPPER[i])
}

AR_CIs <- gibbs_dep_CIs %>% 
  group_by(WEEK) %>% 
  select(`Dependent`=COVERED) %>% 
  summarize_all(mean)

nonAR_CIs <- gibbs_indep_CIs %>% 
  group_by(WEEK) %>% 
  select(`Independent`=COVERED) %>% 
  summarize_all(mean)

DE_CIs <- dir_est_CIs %>% 
  group_by(WEEK) %>% 
  select(`Direct estimate`=COVERED) %>% 
  summarize_all(mean)

all_CIs <- merge(merge(DE_CIs, nonAR_CIs, by="WEEK"), AR_CIs, by="WEEK")
overall_CI <- colMeans(all_CIs, na.rm=T)
all_CIs <- rbind(all_CIs, overall_CI)
all_CIs[nrow(all_CIs),"WEEK"] = "Overall"
print(xtable(all_CIs, digits=5, 
            caption=paste("Mean CI coverage for all models averaged over",Nsim,"simulations in", case, "case")), 
      include.rownames=F, file=paste0(table_dir,"CI_cov.tex"))

overall_table <- cbind(MSE=overall_mse, 
                       MAB=overall_mab, 
                       `CI coverage`=overall_CI)[-1,]

dir_est_score <- tibble(truth,  score=int_score(dirU,dirL,truth$P)) %>%
   group_by(WEEK) %>%
   summarize_all(mean, na.rm=T) %>%
   select(WEEK, `Direct estimate`=score)

nonAR_score <- tibble(truth, score=int_score(nonARU,nonARL,truth$P)) %>%
   group_by(WEEK) %>%
   summarize_all(mean) %>%
   select(WEEK, `Independent`=score)

AR_score <- tibble(truth, score=int_score(fullU,fullL,truth$P)) %>%
   group_by(WEEK) %>%
   summarize_all(mean) %>%
   select(WEEK, `Dependent`=score)

all_scores <- merge(merge(dir_est_score, nonAR_score, by="WEEK"), AR_score, by="WEEK")
overall_score <- colMeans(all_scores)
all_scores <- rbind(all_scores, overall_score)
all_scores[nrow(all_scores),"WEEK"] = "Overall"
print(xtable(all_scores, digits=5,
       caption=paste("Interval scores for models and direct estimates averaged over",
                     Nsim,"simulations in", case, "case")),
      include.rownames=F, file=paste0(table_dir,"score_table.tex"))


##plots
MSE.plot.df <- bind_rows(
                 tibble(truth,  sqr_err=(truth$P-dirEsts)^2, MODEL="Dir. Est."),
                 tibble(truth,  sqr_err=(truth$P-nonARests)^2, MODEL=ind_model),
                 tibble(truth,  sqr_err=(truth$P-fullEsts)^2, MODEL=dep_model)
               )  %>%
               group_by(WEEK,AREA,MODEL) %>% 
               summarize(MSE=mean(sqr_err, na.rm=T))

MSE.plot.df$MODEL <- factor(MSE.plot.df$MODEL, levels=c("Dir. Est.", ind_model, dep_model))

states <- tolower(levels(HPS_pop_long$AREA))
state_poly <- map_data("state") %>% filter(region %in%  states)


##Relative MSE
relative_MSE.df <- MSE.plot.df %>% 
                    group_by(AREA,MODEL) %>% 
                    summarize(mean(MSE)) %>% 
                    pivot_wider(names_from=MODEL, values_from="mean(MSE)") %>%  
                    ##annoying syntax to evaluate string argument in dplyr https://stackoverflow.com/a/56546176
                    summarize(ratio=!! rlang::sym(dep_model)/ !! rlang::sym(ind_model)) %>%
                    select(AREA, ratio)

relative_MSE.df$region<- tolower(relative_MSE.df$AREA)
relative_MSE.df <- merge(relative_MSE.df, state_poly, by="region")

lower_lim <- .95*min(relative_MSE.df$ratio)
upper_lim <- 1.05*max(relative_MSE.df$ratio)

#https://stackoverflow.com/a/28963405
centroids <- setNames(do.call("rbind.data.frame",
                              by(relative_MSE.df,
                                 relative_MSE.df$group,
                                 function(x) {Polygon(x[c('long', 'lat')])@labpt})), c('long', 'lat'))
centroids$label <- relative_MSE.df$AREA[match(rownames(centroids),
                                              relative_MSE.df$group)]
centroids$label <- state.abb[match(centroids$label, state.name)]
centroids <- centroids%>% arrange(label)
centroids <- centroids[-c(17, 19, 22, 28, 30, 37, 39, 40, 51, 52, 55:58),] #hack to remove duplicates
centroids[8, 1] <- centroids[8, 1] + .6

relative_MSE.df %>%
  ggplot(mapping = aes(x=long, y=lat,
                       fill=ratio, group=group)) +
    geom_polygon(color="black", size=0.15) +
  with(centroids, annotate(geom="text", x = long, y=lat, label = label, size = 5)) +
#    viridis::scale_fill_viridis(name=paste0("MSE ratio\n(", dep_model,"/", ind_model, ")"), 
#                               discrete=F)+
    scale_fill_viridis_c(name=paste0("MSE ratio\n(", dep_model,"/", ind_model, ")"),
                                     limits=c(lower_lim, upper_lim), oob=scales::squish) +
    #ggtitle(paste0("Ratio of ", ind_model, " MSE to ", dep_model, " MSE")) +
    labs(fill="Estimate")+
    coord_map() +
    ggthemes::theme_map() +
#   theme(legend.background=element_blank())
    theme(legend.position="left",
             plot.margin=grid::unit(c(1,0,0,0), "mm"))
ggsave(paste0(plot_dir,"MSE_reduction_labeled.pdf"))

system2(command = "pdfcrop", 
        args    = c(paste0(plot_dir,"MSE_reduction_labeled.pdf"),
                    paste0(plot_dir,"MSE_reduction_labeled.pdf")))                    


