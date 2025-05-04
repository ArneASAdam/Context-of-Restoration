# Script is used to analyse importance of selection of variables to restoration benefits for scenarios with/without disturbances using random forest.
# Additionally, partial dependency plots are contructed using normalised and raw data (corresponding to figure 4 in manuscript)
# internal note, script is derived from Rwinners_VR_1dist_1int_191124.R

library(tidyverse) ## for data wrangling
library(cowplot) ## for multi-panel plots using ggplot
library(lme4) ## for (generalised) mixed-effects models
library(lmerTest) ## for getting p-values from lmer models
library(emmeans) ## for posthoc multiple comparisons
library(readxl) ## to import xl
library(MASS)
library(forecast)
library(r2glmm)
library(broom)
library(ggplot2)
library('car')
library('corrplot')
library('caret')
library('pdp')
library('gridExtra')
library('randomForest')

# set directory
setwd('C:/Users/uqaadam5/OneDrive - The University of Queensland/Desktop/UQ postdoc/Reef_emulator/Reef_emulator_Jan_2024/rme_ml_2024_01_08/1int_1cyc_170824/')
setwd('C:/Users/uqaadam5/OneDrive - The University of Queensland/Desktop/UQ postdoc/Reef_emulator/Reef_emulator_Jan_2024/rme_ml_2024_01_08/1int_1dhw_170824/')

# depending on what scenarios you want to analyses select the proper input compiled data file
data_int_scen<-read.csv("Extract_final_OC_VR_1cyc_1int_3LS_nodist_norm_3SD_3PAR_170824.txt") #scenarios with no disturbances
data_int_scen<-read.csv("Extract_final_OC_VR_1cyc_1int_170824_3larv_3SD_3PAR_cor_norm.txt") #scenarios with only cyclones
data_int_scen<-read.csv("Extract_final_OC_VR_1dhw_1int_3LS_5_10_20_200824_cor_norm.txt") #scenarios with only bleaching
head(data_int_scen)

data_int_scen_norm<-data_int_scen[,c(107:212)] #only normalised data
data_int_scen_raw= data_int_scen[,c(1:106)] # only raw data
head(data_int_scen_norm)
head(data_int_scen_raw)

ssp19_raw= data_int_scen[,c(1:106)] # only raw data

#include all informative variables from 1 disturbance/1 intervention perspective
#no disturbances
Nodist_1cyc_zscore_only_nocorr_170824 <- data_int_scen_norm[ , c("zcoral_cover_int_ref_censusyr_mean","ztot_larv_in_param_index","zcoral_cover_startyr","zsize_dist_index","zself_ret","zspecies_comp_index","zsurface_WQ","zprop_rest","zstockdens","zyears_postdeploy")] #option 1

# cyclones
cyc_zscore_only_nocorr_170824 <- data_int_scen_norm[ , c("zcoral_cover_int_ref_censusyr_mean","ztot_larv_in_param_index","zcoral_cover_startyr","zsize_dist_index","zself_ret","zspecies_comp_index","zsurface_WQ","zprop_rest","zstockdens","ztot_cycl_till_startyr","zfirst_cyc_aft_first_int","zyears_postdeploy")] #option 1

# bleaching event
dhw_zscore_only_nocorr_5_10_20_200824 <- data_int_scen_norm[ , c("zcoral_cover_int_ref_censusyr_mean","ztot_larv_in_param_index","zcoral_cover_startyr","zsize_dist_index","zself_ret","zspecies_comp_index","zsurface_WQ","zprop_rest","zstockdens","zheattol","ztot_bleaching_till_startyr","zfirst_bleaching_aft_first_int","zyears_postdeploy")] #option 1

# test with raw data focussing on data with no disturbances
Nodist_1cyc_raw_only_nocorr_170824 <- data_int_scen[ , c("coral_cover_int_ref_censusyr_mean","tot_larv_in_param_index","coral_cover_startyr","size_dist_index","self_ret","species_comp_index","surface_WQ","prop_rest","stockdens","years_postdeploy")] #option 1

# Run model to assess variable importance on normalised data for each of the three analyses (nodist/cyc/bleach)

# Nodist
model <- randomForest(
  formula = zcoral_cover_int_ref_censusyr_mean ~ .,
  data = Nodist_1cyc_zscore_only_nocorr_170824
)

saveRDS(model,"Nodist_cyc_nocorr_rf_3selfret_3larval_3stockdens_3PAR_norm_suppl_noHT_170824.R")

#with cyclone, was analysed on the HPC as it is computationally heavy

model <- randomForest(
  formula = zcoral_cover_int_ref_censusyr_mean ~ .,
  data = cyc_zscore_only_nocorr_170824
)

saveRDS(model,"cyc_nocorr_rf_3SD_3LV_3SD_3PAR_5_10_20_zscore_180824.R")

model<-readRDS("cyc_nocorr_rf_3SD_3LV_3SD_3PAR_5_10_20_zscore_180824.R") 


# with bleaching

model <- randomForest(
  formula = zcoral_cover_int_ref_censusyr_mean ~ .,
  data = dhw_zscore_only_nocorr_5_10_20_200824
)

saveRDS(model,"dhw_nocorr_rf_3SD_3LV_3SD_3PAR_5_10_20_zscore_220824.R")

par(mfrow=c(1,1))

# Visualise variable importance for all 3 scenarios
t<-varImpPlot(model, 
              sort=T, 
              bg = "skyblue",
              main="Intervention scenarios (no disturbance) Variable Importance")
##################################################################################

# Plot variable importance using partial dependency plots to assess how a variable changes restoration benefits while keep the other metrics and constant levels

library("pdp")


par(mfrow=c(4,2))

plotlist = list()


#08/08/24 & 17/08/24

raw_var_names=sub('.', '', rownames(model$importance)) #remove z from the values

for (i in 1:8) { #no disturbance, for every unique variable calculate partial dependency relative to the response variable witch is coral cover benefit
  
  plot_var_name <- str_c(c("ggplot", raw_var_names[i]), collapse = "_")  
  pred.grid_final=expand_grid(unique(Nodist_1cyc_zscore_only_nocorr_170824[noquote(rownames(model$importance)[i])]),unique(Nodist_1cyc_zscore_only_nocorr_170824$zyears_postdeploy ))
  
  colnames(pred.grid_final)=c(rownames(model$importance)[i],"zyears_postdeploy")
  forest_stats<-partial(model,pred.data=Nodist_1cyc_zscore_only_nocorr_170824,pred.var=c(rownames(model$importance)[i],"zyears_postdeploy"),pred.grid = pred.grid_final)
  
  forest_stats[!duplicated(forest_stats), ]

  # calculate mean and standard dev from raw data
  mean_var=mean(data_int_scen_raw[[(raw_var_names[i])]])
  std_var= sd(data_int_scen_raw[[(raw_var_names[i])]])
  
  mean_resp=mean(data_int_scen_raw$coral_cover_int_ref_censusyr_mean)
  std_resp=sd(data_int_scen_raw$coral_cover_int_ref_censusyr_mean)
  
  mean_cens=mean(data_int_scen_raw$years_postdeploy)
  std_cens=sd(data_int_scen_raw$years_postdeploy)

  v<-unique(forest_stats[noquote(rownames(model$importance)[i])])
  
  forest_stats$zyears_postdeploy=rep(0:14,nrow(v))

  forest_stats_non_norm=forest_stats
  forest_stats_non_norm[[(raw_var_names[i])]]=(forest_stats[[noquote(rownames(model$importance)[i])]]*std_var)+mean_var
  
  forest_stats_non_norm$coral_cover_int_ref_censusyr_mean=(forest_stats$yhat* std_resp)+mean_resp
  max_coralcover=max(forest_stats_non_norm$coral_cover_int_ref_censusyr_mean)
  
  if(i==1){
    forest_stats_non_norm$tot_larv_in_param_index=as.factor(forest_stats_non_norm$tot_larv_in_param_index)
    a1=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = tot_larv_in_param_index), linewidth=0.75)+geom_point(aes(color = tot_larv_in_param_index),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Yearly external larval supply (larvae/m2)') +

      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    
    assign(plot_var_name, a1)}
  
  else if(i==2){ 
    forest_stats_non_norm$coral_cover_startyr=as.factor(forest_stats_non_norm$coral_cover_startyr)
    a2=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = coral_cover_startyr), linewidth=0.75)+geom_point(aes(color = coral_cover_startyr),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Coral cover year before deployment (%)') +
      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    assign(plot_var_name, a2)}
  
  else if(i==3){ 
    forest_stats_non_norm$size_dist_index=as.factor(forest_stats_non_norm$size_dist_index)
    a3=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = size_dist_index), linewidth=0.75)+geom_point(aes(color = size_dist_index),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Size distribution index') +
      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    assign(plot_var_name, a3)}
  
  else if(i==4){ 
    forest_stats_non_norm$self_ret=forest_stats_non_norm$self_ret*100
    forest_stats_non_norm$self_ret=as.factor(forest_stats_non_norm$self_ret)
    
    a4=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = self_ret), linewidth=0.75)+geom_point(aes(color = self_ret),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Self retention (%)') +
      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    assign(plot_var_name, a4)}
  
  else if(i==5){ 
    forest_stats_non_norm$species_comp_index=as.factor(forest_stats_non_norm$species_comp_index)
    a5=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = species_comp_index), linewidth=0.75)+geom_point(aes(color = species_comp_index),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Species composition index') +
      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    assign(plot_var_name, a5)}
  
  else if(i==6){ 
    forest_stats_non_norm$surface_WQ=as.factor(forest_stats_non_norm$surface_WQ)
    a6=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = surface_WQ), linewidth=0.75)+geom_point(aes(color = surface_WQ),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Sedimentation (mg/l)') +
      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    assign(plot_var_name, a6)}
  
  else if(i==7){ 
    forest_stats_non_norm$prop_rest=as.factor(forest_stats_non_norm$prop_rest)
    a7=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = prop_rest), linewidth=0.75)+geom_point(aes(color = prop_rest),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Proportion of reef restored (%)') +
      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    assign(plot_var_name, a7)}
  
  else if(i==8){ 
    forest_stats_non_norm$stockdens=as.factor(forest_stats_non_norm$stockdens)
    a8=ggplot(forest_stats_non_norm,aes(zyears_postdeploy,coral_cover_int_ref_censusyr_mean))+geom_line(aes(color = stockdens), linewidth=0.75)+geom_point(aes(color = stockdens),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      
      labs(x = "Years post deployment", y = "Mean coral cover benefit (%)")+ labs(color='Stocking density (ind./m2)') +
      theme(axis.text.x = element_text( colour="black",size = 11))+
      theme(axis.text.y = element_text(colour="black", size = 11))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 12))+
      theme(axis.title.y = element_text(face="bold", size = 12))+
      theme (legend.title=element_text (face="bold", size=11))+
      theme (legend.title=element_blank ())+
      guides(colour = guide_legend(title.position = "top", label.hjust=0.5))+
      theme(legend.position="none")+ 
      expand_limits(y = c(0, 5))+ scale_y_continuous(breaks=c(1,2,3, 4,5, 6,7, 8))+
      
      
      scale_x_continuous(breaks=c(0,2, 4, 6, 8,10,12,14))
    
    assign(plot_var_name, a8)}

}

margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))


#reordered 19/11/24 5,10,20%
p <- grid.arrange(grobs=list(ggplot_tot_larv_in_param_index,
                             ggplot_species_comp_index,
                             ggplot_stockdens,
                             ggplot_prop_rest,
                             ggplot_self_ret,
                             ggplot_coral_cover_startyr,
                             ggplot_surface_WQ,
                             ggplot_size_dist_index),nrow=4,ncol=2)


#########################################################
######### Test relationships on raw data ################


#08/08/24 and 17/08/24- raw data

#########################################################
######### Test relationships on raw data (Supps figure S3) ################
Nodist_1cyc_raw_only_nocorr_140424<-Nodist_1cyc_raw_only_nocorr_080824
Nodist_1cyc_raw_only_nocorr_140424<-Nodist_1cyc_raw_only_nocorr_170824

# create subset of larval supply
names<-colnames(Nodist_1cyc_raw_only_nocorr_140424)

raw_regress=table()

for (var in 2:9) {
  plot_var_name <- str_c(c("ggplot", names[var]), collapse = "_")  
  
  raw_regress=data.frame(expand_grid(unique(Nodist_1cyc_raw_only_nocorr_170824[noquote(names[var])]),unique(Nodist_1cyc_raw_only_nocorr_170824$years_postdeploy )[c(1:15)]))
  colnames(raw_regress)[2]="year_postdeploy"
  for (comb in 1:nrow(raw_regress)) {
    #   for (comb in 1:5) {
    Nodist_1cyc_raw_only_nocorr_170824_LS_subset=Nodist_1cyc_raw_only_nocorr_170824[Nodist_1cyc_raw_only_nocorr_170824[noquote(names[var])]==raw_regress[noquote(names[var])][comb,],]
    nrow(Nodist_1cyc_raw_only_nocorr_170824_LS_subset)
    
    Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub=Nodist_1cyc_raw_only_nocorr_170824_LS_subset[Nodist_1cyc_raw_only_nocorr_170824_LS_subset$years_postdeploy==raw_regress$year_postdeploy[comb],]
    unique(Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub$years_postdeploy)
    
    # calculate mean of all variables
    #[[noquote(rownames(model$importance)[i])]]
    Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub_mean<-sapply(Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub,FUN=mean)
    Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub_sd<-sapply(Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub,FUN=sd)
    
    raw_regress$Mean_benefit[comb]<-Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub_mean[1]
    raw_regress$std_benefit[comb]<-Nodist_1cyc_raw_only_nocorr_170824_LS_sub_census_sub_sd[1]
  }
  i=var-1
  if(i==1){
    raw_regress$tot_larv_in_param_index=as.factor(raw_regress$tot_larv_in_param_index)
    a1=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = tot_larv_in_param_index), linewidth=0.75)+geom_point(aes(color = tot_larv_in_param_index),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Yearly external larval supply (larvae/m2)') +

      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 12))+ 
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    
    
    assign(plot_var_name, a1)}
  
  else if(i==2){ 
    raw_regress$coral_cover_startyr=as.factor(raw_regress$coral_cover_startyr)
    a2=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = coral_cover_startyr), linewidth=0.75)+geom_point(aes(color = coral_cover_startyr),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Coral cover year before deployment (%)') +
      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    assign(plot_var_name, a2)}
  
  else if(i==3){ 
    raw_regress$size_dist_index=as.factor(raw_regress$size_dist_index)
    a3=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = size_dist_index), linewidth=0.75)+geom_point(aes(color = size_dist_index),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Size distribution index') +
      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    assign(plot_var_name, a3)}
  
  else if(i==4){ 
    raw_regress$self_ret=raw_regress$self_ret*100
    raw_regress$self_ret=as.factor(raw_regress$self_ret)
    
    a4=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = self_ret), linewidth=0.75)+geom_point(aes(color = self_ret),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Self retention (%)') +
      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    assign(plot_var_name, a4)}
  
  else if(i==5){ 
    raw_regress$species_comp_index=as.factor(raw_regress$species_comp_index)
    a5=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = species_comp_index), linewidth=0.75)+geom_point(aes(color = species_comp_index),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Species composition index') +
      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    assign(plot_var_name, a5)}
  
  else if(i==6){ 
    raw_regress$surface_WQ=as.factor(raw_regress$surface_WQ)
    a6=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = surface_WQ), linewidth=0.75)+geom_point(aes(color = surface_WQ),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Sedimentation (mg/l)') +
      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    assign(plot_var_name, a6)}
  
  else if(i==7){ 
    raw_regress$prop_rest=as.factor(raw_regress$prop_rest)
    a7=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = prop_rest), linewidth=0.75)+geom_point(aes(color = prop_rest),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Proportion of reef restored (%)') +
      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    assign(plot_var_name, a7)}
  
  else if(i==8){ 
    raw_regress$stockdens=as.factor(raw_regress$stockdens)
    a8=ggplot(raw_regress,aes(year_postdeploy ,Mean_benefit))+geom_line(aes(color = stockdens), linewidth=0.75)+geom_point(aes(color = stockdens),shape=19, size=2)+scale_color_viridis_d()+ theme_bw()+
      geom_errorbar(aes(ymin=Mean_benefit-std_benefit, ymax=Mean_benefit+std_benefit), width=.2,
                    position=position_dodge(0.05))+  
      
      labs(x = "Years post deployment", y = "Coral cover benefit (%)")+ labs(color='Stocking density (ind./m2)') +
      theme(axis.text.x = element_text( colour="black",size = 9))+
      theme(axis.text.y = element_text(colour="black", size = 9))+
      theme (legend.text=element_text (size=11))+
      theme(axis.title.x = element_text(face="bold", size = 10))+
      theme(axis.title.y = element_text(face="bold", size = 10))+
      theme (legend.title=element_text (size=11))+
      theme(legend.position="none")+ expand_limits(y = c(0, 16))+ scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
      scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14))
    
    assign(plot_var_name, a8)}
  
}

margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))


#Reordered 17/08/24

p <- grid.arrange(grobs=list(ggplot_tot_larv_in_param_index,
                             ggplot_species_comp_index,
                             ggplot_prop_rest,
                             ggplot_stockdens,
                             ggplot_self_ret,
                             ggplot_coral_cover_startyr,
                             ggplot_surface_WQ,
                             ggplot_size_dist_index),nrow=4,ncol=2)
#################################################################
  