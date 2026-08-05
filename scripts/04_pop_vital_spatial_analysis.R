################################################################################
# 
#  Population vital rates spatial analysis      
#
#  by Anonymized                           
#
################################################################################

# This script maps mussel population vital rate predictions from the structural
# equation modelling onto spatial rasters. These rasters are used in 
# correlation analyses with predicted mussel density rasters from Annis et al.
# 2022.

################################################################################
##  House Keeping  ##
################################################################################

### Clear workspace  ###
rm(list = ls())

### Load in packages  ###
library(raster)
library(dplyr)
library(readxl)
library(corrplot)
library(spdep)


################################################################################
##  Load data  ##
################################################################################

### Predicted mussel density rasters ###########################################
# Load in predicted mussel density rasters from Annis et al., 2022 for use as
# response in correlation analysis
cp_sdm <- 
  raster::raster("Mussel SDM rasters/mussel_densities_predicted_cannons.tif")  
lp_sdm <- 
  raster::raster("Mussel SDM rasters/mussel_densities_predicted_sapelo.tif")  


### SEM outputs ################################################################
# Load in outputs from SEM analysis to create prediction rasters of population
# vital rates. These rates will be used in correlation analysis with predicted
# mussel densities
growth_out <- readRDS("model_out/growth_SEM_final_output.rds")
pred_out <- readRDS("model_out/predation_SEM_final_output.rds")
recruit_out <- readRDS("model_out/recruitment_SEM_final_output.rds")


###  Explanatory rasters #######################################################
# Load in marsh-level explanatory rasters and convert into meters (1/3.281).
# These rasters will be used with SEM model exports to create prediction rasters
# of population vital rates. Prediction rasters will need to be create from
# both study sites.

## Cannons point ##
cp_dem <- 
  raster::raster("Mussel SDM rasters/elevation_distance_cannons.tif")/3.281  
cp_head <- 
  raster::raster("Mussel SDM rasters/head_distance_cannons.tif ")/3.281  
cp_main <- 
  raster::raster("Mussel SDM rasters/subtidal_distance_cannons.tif")/3.281  


## Sapelo island ##
lp_dem <- 
  raster::raster("Mussel SDM rasters/elevation_distance_sapelo.tif")/3.281  
lp_head <- 
  raster::raster("Mussel SDM rasters/head_distance_sapelo.tif")/3.281  
lp_main <- 
  raster::raster("Mussel SDM rasters/subtidal_distance_sapelo.tif")/3.281  


### Non-spatial explanatory variables ##########################################
# Explanatory variables such as temperature and mussel size did not vary 
# spatially in SEMs. These either varied across time or randomly across 
# experimental units. As such, mean values of these variables will be used in 
# vital rate rasters.

## Raw data  ##
growth_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 1)
pred_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 2)
recruit_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 3)


## growth  ##
growth_temp <- growth_data %>% 
  mutate(site = substr(plot,1,1)) %>% 
  group_by(site) %>% 
  summarise(temp = mean(mean_water_temp,na.rm = T),
            size = mean(mussel_size,na.rm = T))

# Cannon's Point
cp_gtemp <- growth_temp$temp[growth_temp$site == "c"]  # water temperature
cp_gsize <- growth_temp$size[growth_temp$site == "c"]  # mussel size

# Sapelo island
lp_gtemp <- growth_temp$temp[growth_temp$site == "l"]  # water temperature
lp_gsize <- growth_temp$size[growth_temp$site == "l"]  # museel size


## Predation  ##
pred_temp <- pred_data %>% 
  mutate(site = substr(plot,1,1)) %>% 
  group_by(site) %>% 
  summarise(temp = mean(temp,na.rm = T))

# Cannon's Point
cp_temp <- pred_temp$temp[pred_temp$site == "c"]  # water temperature

# Sapelo island
lp_temp <- pred_temp$temp[pred_temp$site == "l"]  # water temperature


## Recruit  ##
recruit_temp <- recruit_data %>% 
  mutate(site = substr(plot,1,1)) %>% 
  group_by(site) %>% 
  summarise(temp = mean(mean_water_temp,na.rm = T))

# Cannon's Point
cp_rtemp <- recruit_temp$temp[recruit_temp$site == "c"]  # water temperature

# Sapelo island
lp_rtemp <- recruit_temp$temp[recruit_temp$site == "l"]  # water temperature


################################################################################
##  Vital rate prediction rasters  ##
################################################################################

# Prediction rasters of each vital rate are created using the SEM outputs and 
# spatial predictors. Seasonal and mussel specific predictors are held at mean
# values. Rasters were used to create Figure 5 of the main text in Qgis.

### Growth Rate ################################################################

## Cannon's Point ##

# Salinity prediction
cp_gsal <- growth_out$mean$b5.0 + 
  cp_gtemp*growth_out$mean$b5.2+ cp_gtemp*cp_gtemp*growth_out$mean$b5.2q +
  cp_main*growth_out$mean$b5.14

# Submergence prediction
cp_gsub <- growth_out$mean$b3.0 + 
  growth_out$mean$b3.12*cp_dem + 
  growth_out$mean$b3.13*cp_head

# Flow prediction
cp_gflow <- growth_out$mean$b4.0 + 
  growth_out$mean$b4.3*cp_gsub + 
  growth_out$mean$b4.13*cp_head + 
  growth_out$mean$b4.14*cp_main + 
  growth_out$mean$b4.12*cp_dem   

# Growth prediction
cp_growth <- growth_out$mean$b0.0 + 
  growth_out$mean$b0.2*cp_gtemp +
  growth_out$mean$b0.3*cp_gsub +
  growth_out$mean$b0.4*cp_gflow + 
  growth_out$mean$b0.5*cp_gsal + 
  growth_out$mean$b0.9*cp_gsize

# Plot results
plot(cp_growth)

# # Export Growth Raster
# raster::writeRaster(
#   cp_growth,
#   "Prediction rasters/growth_prediction_cannons.tif", 
#   overwrite=TRUE
#   )


## Sapelo island  ##

# Salinity prediction
lp_gsal <- growth_out$mean$b5.0 + 
  lp_gtemp*growth_out$mean$b5.2+ lp_gtemp*lp_gtemp*growth_out$mean$b5.2q +
  lp_main*growth_out$mean$b5.14

# Submergence prediction
lp_gsub <- growth_out$mean$b3.0 + 
  growth_out$mean$b3.12*lp_dem + 
  growth_out$mean$b3.13*lp_head

# Flow prediction
lp_gflow <- growth_out$mean$b4.0 + 
  growth_out$mean$b4.3*lp_gsub + 
  growth_out$mean$b4.13*lp_head + 
  growth_out$mean$b4.14*lp_main + 
  growth_out$mean$b4.12*lp_dem   

# Growth prediction
lp_growth <- growth_out$mean$b0.0 + 
  growth_out$mean$b0.2*lp_gtemp +
  growth_out$mean$b0.3*lp_gsub +
  growth_out$mean$b0.4*lp_gflow + 
  growth_out$mean$b0.5*lp_gsal + 
  growth_out$mean$b0.9*lp_gsize

# Plot results
plot(lp_growth)

# # Export growth raster
# raster::writeRaster(
#   lp_growth,
#   "Prediction rasters/growth_prediction_sapelo.tif",
#   overwrite=TRUE
#   )

### Predation ##################################################################

## Cannon's Point ##

# Salinity prediction
cp_psal <- pred_out$mean$b6.0 + 
  cp_temp*pred_out$mean$b6.2+ cp_temp*cp_temp*pred_out$mean$b6.2q +
  cp_main*pred_out$mean$b6.5

# SA height prediction
cp_psah <- pred_out$mean$b10.0 +   
  pred_out$mean$b10.5*cp_main

# Submergance prediction
cp_psub <- pred_out$mean$b11.0 + 
  pred_out$mean$b11.2*cp_temp + pred_out$mean$b11.2q*cp_temp*cp_temp +
  pred_out$mean$b11.3*cp_dem + 
  pred_out$mean$b11.4*cp_head  

# Predation prediction
cp_pred <- (1/(1 + exp(-(pred_out$mean$b0.0 +
                           pred_out$mean$b0.2*cp_temp + 
                           pred_out$mean$b0.6*cp_psal +
                           pred_out$mean$b0.10*cp_psah +
                           pred_out$mean$b0.11*cp_psub)  )))

# Plot results
plot(cp_pred)

# # Export
# raster::writeRaster(
#   cp_pred,
#   "Prediction rasters/predation_prediction_cannons.tif"
#   )


## Sapelo island ##

# Salinity prediction
lp_psal <- pred_out$mean$b6.0 + 
  lp_temp*pred_out$mean$b6.2+ lp_temp*lp_temp*pred_out$mean$b6.2q +
  lp_main*pred_out$mean$b6.5

# SA height prediction
lp_psah <- pred_out$mean$b10.0 +   
  pred_out$mean$b10.5*lp_main

# Submergance prediction
lp_psub <- pred_out$mean$b11.0 + 
  pred_out$mean$b11.2*lp_temp + pred_out$mean$b11.2q*lp_temp*lp_temp +
  pred_out$mean$b11.3*lp_dem + 
  pred_out$mean$b11.4*lp_head  

# Predation prediction
lp_pred <- (1/(1 + exp(-(pred_out$mean$b0.0 +
                           pred_out$mean$b0.2*lp_temp + 
                           pred_out$mean$b0.6*lp_psal +
                           pred_out$mean$b0.10*lp_psah +
                           pred_out$mean$b0.11*lp_psub)  )))

# plot results
plot(lp_pred)

# # Export
# raster::writeRaster(
#   lp_pred,
#   "Prediction rasters/predation_prediction_sapelo.tif"
#   )


### Recruitment ################################################################

## Cannon's Point ##

# Salinity prediction
cp_rsal <- recruit_out$mean$b5.0 +
  recruit_out$mean$b5.2*cp_rtemp + recruit_out$mean$b5.2q*cp_rtemp*cp_rtemp +
  recruit_out$mean$b5.14*cp_main

# SA density prediction
cp_rsad <- exp(recruit_out$mean$b10.0 +
  recruit_out$mean$b10.13*cp_head +
  recruit_out$mean$b10.14*cp_main)

# SA height prediction
cp_rsah <- recruit_out$mean$b11.0 +   
  recruit_out$mean$b11.12*cp_dem+
  recruit_out$mean$b11.14*cp_main+  
  recruit_out$mean$b11.2*cp_rtemp

# Crab burrow prediction
cp_rcrab <- exp(recruit_out$mean$b15.0 +   
  recruit_out$mean$b15.12*cp_dem+
  recruit_out$mean$b15.13*cp_head+
  recruit_out$mean$b15.14*cp_main +
  recruit_out$mean$b15.10*cp_rsad +  
  recruit_out$mean$b15.11*cp_rsah)

# Recrtuitment Prediction 
cp_recruit <- exp(recruit_out$mean$b0.0 + 
  recruit_out$mean$b0.2*cp_rtemp + 
  recruit_out$mean$b0.5*cp_rsal + 
  recruit_out$mean$b0.10*cp_rsad +
  recruit_out$mean$b0.11*cp_rsah +
  recruit_out$mean$b0.15*cp_rcrab +
  recruit_out$mean$b0.14*cp_main)

# Plot results
plot(cp_recruit)

# # Export
# raster::writeRaster(
#   cp_recruit,
#   "Prediction rasters/recruitment_prediction_cannons.tif"
#   )


## Sapelo island ##

# Salinity prediction
lp_rsal <- recruit_out$mean$b5.0 +
  recruit_out$mean$b5.2*lp_rtemp + recruit_out$mean$b5.2q*lp_rtemp*lp_rtemp +
  recruit_out$mean$b5.14*lp_main

# SA density prediction
lp_rsad <- exp(recruit_out$mean$b10.0 +
                 recruit_out$mean$b10.13*lp_head +
                 recruit_out$mean$b10.14*lp_main)

# SA height prediction
lp_rsah <- recruit_out$mean$b11.0 +   
  recruit_out$mean$b11.12*lp_dem+
  recruit_out$mean$b11.14*lp_main+  
  recruit_out$mean$b11.2*lp_rtemp

# Crab burrow prediciton
lp_rcrab <- exp(recruit_out$mean$b15.0 +   
                  recruit_out$mean$b15.12*lp_dem+
                  recruit_out$mean$b15.13*lp_head+
                  recruit_out$mean$b15.14*lp_main +
                  recruit_out$mean$b15.10*lp_rsad +  
                  recruit_out$mean$b15.11*lp_rsah)

# Recruitment prediction
lp_recruit <- exp(recruit_out$mean$b0.0 + 
                    recruit_out$mean$b0.2*lp_rtemp + 
                    recruit_out$mean$b0.5*lp_rsal + 
                    recruit_out$mean$b0.10*lp_rsad +
                    recruit_out$mean$b0.11*lp_rsah +
                    recruit_out$mean$b0.15*lp_rcrab +
                    recruit_out$mean$b0.14*lp_main)

# Plot results
plot(lp_recruit)

# # Export
# raster::writeRaster(
#   lp_recruit,
#   "Prediction rasters/recruitment_prediction_cannons.tif"
#   )


################################################################################
##  Vital rate prediction correlation analysis  ##
################################################################################

### Prepare raster data for analysis ###########################################

## Resample rasters  ##
# Sapelo island's predicted mussel density raster does not line up with
# explanatory rasters. Resample this raster to enable correlation analysis
lp_sdm_resampled <- resample(lp_sdm, lp_growth, method = "bilinear")


## Combine both sites into one raster  ##
sdm<- mosaic(cp_sdm,lp_sdm_resampled, fun = "mean")  # prediction mussel density
g<- mosaic(cp_growth,lp_growth, fun = "mean")  # predicted growth
p<- mosaic(cp_pred,lp_pred, fun = "mean")  # predicted predation
r<- mosaic(cp_recruit,lp_recruit, fun = "mean")  # predicted recruitment


## Create dataframe with values from all rasters  ##
cor_df <- cbind(values(sdm),values(g),values(p),values(r))
colnames(cor_df) <- c("Density","Growth","Predation","Recrutiment")


### Analysis ###################################################################

## Correlation matrix  ##
cor_mat_sp <-cor(cor_df,use = "na.or.complete",method = "spearman")


## Individual pairwise correlations and significance testing  ##
cor.test(values(sdm), values(g),use = "na.or.complete",method = "spearman")
cor.test(values(sdm), values(r),use = "na.or.complete",method = "spearman")
cor.test(values(sdm), values(p),use = "na.or.complete",method = "spearman")

### Plotting  ##################################################################
# This code creates the correlation plot used in Figure 6 of the main text
png(filename = "plots/cor/spear_cor_plot.png", width = 800, height = 800)
corrplot(cor_mat_sp, diag = F, type = "lower")
dev.off()


################################################################################
##  Moran's I   ##
################################################################################

# Test for spatial autocorrelation in residuals of each vital rate SEM

## load in Site coordinates  ##
coors <- read.csv("plot_parameters/site_coors.csv")

## load in residuals  ##
resid_files <- list.files("model_out","resid")
resid_list <- lapply(
  resid_files, 
  function(i) readRDS(file.path("model_out",i))
  )
names(resid_list) <- resid_files

## Run Moran I on each model  ##
moran_list <-lapply(resid_list[-3], function(r){
  
  # add resids to points
  sites <- coors %>% 
    mutate(plot = casefold(plot)) %>% 
    left_join(r,by="plot")  %>% 
    summarise(
      resid = mean(resid,na.rm=T),
      .by = c(plot,lat,lon)
    )

  # Find distance based neareast neighbor
  coords <- as.matrix(sites[, c("lon","lat")]) 
  lapply(1:10, function(k){
    knn <- knearneigh(coords, k = k)
    nb  <- knn2nb(knn)
    lw  <- nb2listw(nb, style = "W")
    
    # Moran I
    set.seed(123)
    mi <- moran.mc(
      sites$resid, 
      lw, 
      nsim = 999, 
      zero.policy = TRUE,
      alternative ="two.sided"
      )
    
    data.frame(
      K = k,
      I = mi$statistic,
      rank = mi$parameter,
      p_value = mi$p.value
    )
  }
  ) %>% 
    bind_rows()
  
})

# Export
lapply(1:length(moran_list), function(i){
  name <- paste0(
    "plots/moran_tables/",
    substr(names(moran_list)[i],1,1),
    "_moran_table.csv"
    )
  write.csv(moran_list[[i]],name)
})





