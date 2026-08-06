################################################################################
# 
#  Total effect plots for marsh-level predictors      
#
#  by Anonymized                                 
#
################################################################################

# This script maps predicts the total effect marsh-level variables and temperature
# on each vital rate accounting for all direct and indirect effects, while 
# holding all other variables at mean value. Marsh-level predictors include 
# elevation,distance to subtidal creek (main), and distance to creek head (head).
# These predictors are used to creek total effect plots for each vital rate
# and predictor (4x3=12 plots). These plots are used to create figure 4 in the
# main text


################################################################################
##  House Keeping  ##
################################################################################


### Clear workspace  ###
rm(list = ls())

### Load in packages  ###
library(dplyr)
library(ggplot2)
library(readxl)
library(scales)


################################################################################
##  Load in data  ##
################################################################################


### Growth  ###
growth_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 1)
growth_data$flow_rate <- growth_data$flow_rate
growth_out <- readRDS("model_out/growth_SEM_final_output.rds")
growth_out <- as.data.frame(growth_out$sims.list)

### Predation ###
pred_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 2)
pred_out <- readRDS("model_out/predation_SEM_final_output.rds")
pred_out <- as.data.frame(pred_out$sims.list)

### Recruitment  ###
recruit_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 3)
recruit_out <- readRDS("model_out/recruitment_SEM_final_output.rds")
recruit_out <- as.data.frame(recruit_out$sims.list)

################################################################################
### Growth Rate prediction ####
################################################################################


### Temperature  ###############################################################

### Create vector of values across the range of predictor  ###
temp <- seq(min(growth_data$mean_water_temp,na.rm = T), 
            max(growth_data$mean_water_temp,na.rm = T),
            length.out = 1000)

### Empty Vector to contain predictions  ##
growth.predict.temp.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(growth_out),1)
  growth.predict.temp <- growth_out$b0.0[x] + 
    growth_out$b0.2[x]*temp +
    growth_out$b0.3[x]*mean(growth_data$sub_prop) + 
    growth_out$b0.4[x]*mean(growth_data$flow_rate) +
    growth_out$b0.5[x]*(growth_out$b5.0[x] + 
                          growth_out$b5.2[x]*temp + 
                          growth_out$b5.2q[x]*temp^2+
                          growth_out$b5.14[x]*mean(growth_data$main)) + 
    growth_out$b0.9[x]*mean(growth_data$mussel_size) 
  
  # Bind predictions from each sample into one dataframe
  growth.predict.temp.all <- rbind(growth.predict.temp.all,growth.predict.temp)
  
}

### Mean and credible intervals  ###
growth.temp.md <- apply(growth.predict.temp.all, 2, quantile, probs=0.5)
growth.temp.up <- apply(growth.predict.temp.all, 2, quantile, probs=0.975)
growth.temp.lo <- apply(growth.predict.temp.all, 2, quantile, probs=0.025)

### Elevation ##################################################################

### Create vector of values across the range of predictor  ###
ele <- seq(min(growth_data$ele,na.rm = T), 
           max(growth_data$ele,na.rm = T),
           length.out = 1000)

### Empty Vector to contain predictions  ##
growth.predict.ele.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(growth_out),1)
  growth.predict.ele <- growth_out$b0.0[x] + 
    growth_out$b0.2[x]*mean(growth_data$mean_water_temp,na.rm = T) +
    growth_out$b0.3[x]*(growth_out$b3.0[x] + 
                            growth_out$b3.12[x]*ele + 
                            growth_out$b3.13[x]*mean(growth_data$head)) + 
    growth_out$b0.4[x]*(growth_out$b4.0[x] + 
                            growth_out$b4.3[x]*(growth_out$b3.0[x] + 
                                                    growth_out$b3.12[x]*ele + 
                                                    growth_out$b3.13[x]*mean(growth_data$head)) + 
                            growth_out$b4.12[x]*ele + 
                            growth_out$b4.13[x]*mean(growth_data$head)+
                            growth_out$b4.14[x]*mean(growth_data$main)) +
    growth_out$b0.5[x]*mean(growth_data$mean_sal) + 
    growth_out$b0.9[x]*mean(growth_data$mussel_size)
  
  # Bind predictions from each sample into one dataframe
  growth.predict.ele.all <- rbind(growth.predict.ele.all,growth.predict.ele)
  
}

### Mean and credible intervals  ###
growth.ele.md <- apply(growth.predict.ele.all, 2, quantile, probs=0.5)
growth.ele.up <- apply(growth.predict.ele.all, 2, quantile, probs=0.975)
growth.ele.lo <- apply(growth.predict.ele.all, 2, quantile, probs=0.025)

### Head  ######################################################################

### Create vector of values across the range of predictor  ###
head <- seq(min(growth_data$head,na.rm = T), 
            max(growth_data$head,na.rm = T),
            length.out = 1000)

### Empty Vector to contain predictions  ##
growth.predict.head.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(growth_out),1)
  growth.predict.head<- growth_out$b0.0[x] + 
    growth_out$b0.2[x]*mean(growth_data$mean_water_temp,na.rm = T) +
    growth_out$b0.3[x]*(growth_out$b3.0[x] + 
                          growth_out$b3.12[x]*mean(growth_data$ele) + 
                          growth_out$b3.13[x]*head) + 
    growth_out$b0.4[x]*(growth_out$b4.0[x] + 
                          growth_out$b4.3[x]*(growth_out$b3.0[x] + 
                                                growth_out$b3.12[x]*mean(growth_data$ele) + 
                                                growth_out$b3.13[x]*head) + 
                          growth_out$b4.12[x]*mean(growth_data$ele) + 
                          growth_out$b4.13[x]*head+
                          growth_out$b4.14[x]*mean(growth_data$main)) +
    growth_out$b0.5[x]*mean(growth_data$mean_sal) + 
    growth_out$b0.9[x]*mean(growth_data$mussel_size)
  
  # Bind predictions from each sample into one dataframe
  growth.predict.head.all <- rbind(growth.predict.head.all,growth.predict.head)
  
}

### Mean and credible intervals  ###
growth.head.md <- apply(growth.predict.head.all, 2, quantile, probs=0.5)
growth.head.up <- apply(growth.predict.head.all, 2, quantile, probs=0.975)
growth.head.lo <- apply(growth.predict.head.all, 2, quantile, probs=0.025)

### Main  ######################################################################

### Create vector of values across the range of predictor  ###
main <- seq(min(growth_data$main,na.rm = T), 
            max(growth_data$main,na.rm = T),
            length.out = 1000)

### Empty Vector to contain predictions  ##
growth.predict.main.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(growth_out),1)
  growth.predict.main <- growth_out$b0.0[x] + 
    growth_out$b0.2[x]*mean(growth_data$mean_water_temp,na.rm = T) +
    growth_out$b0.3[x]*mean(growth_data$sub_prop) + 
    growth_out$b0.4[x]*(growth_out$b4.0[x] + 
                            growth_out$b4.3[x]*mean(growth_data$sub_prop) + 
                            growth_out$b4.12[x]*mean(growth_data$ele) + 
                            growth_out$b4.13[x]*mean(growth_data$head)+
                            growth_out$b4.14[x]*main) +
    growth_out$b0.5[x]*(growth_out$b5.0[x] + 
                            growth_out$b5.2[x]*mean(growth_data$mean_water_temp,na.rm = T) + 
                            growth_out$b5.2q[x]*mean(growth_data$mean_water_temp,na.rm = T)^2+
                            growth_out$b5.14[x]*main) + 
    growth_out$b0.9[x]*mean(growth_data$mussel_size)
  
  # Bind predictions from each sample into one dataframe
  growth.predict.main.all <- rbind(growth.predict.main.all,growth.predict.main)
  
}

### Mean and credible intervals  ###
growth.main.md <- apply(growth.predict.main.all, 2, quantile, probs=0.5)
growth.main.up <- apply(growth.predict.main.all, 2, quantile, probs=0.975)
growth.main.lo <- apply(growth.predict.main.all, 2, quantile, probs=0.025)


################################################################################
### Predation prediction ###
################################################################################


### Temperature  ###############################################################

### Create vector of values across the range of predictor  ###
temp.p <- seq(min(pred_data$temp,na.rm = T), 
              max(pred_data$temp,na.rm = T),
              length.out = 1000)

### Empty Vector to contain predictions  ##
pred.predict.temp.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(pred_out),1)
  pred.predict.temp <- pred_out$b0.0[x]+
    pred_out$b0.2[x]*temp.p + 
    pred_out$b0.6[x]*(pred_out$b6.0[x] +
                          pred_out$b6.2[x]*temp.p + 
                          pred_out$b6.2q[x]*temp.p^2+
                          pred_out$b6.5[x]*mean(pred_data$main))  +
    pred_out$b0.10[x]*mean(pred_data$sa_height,na.rm = T) +
    pred_out$b0.11[x]*(pred_out$b11.0[x] +
                           pred_out$b11.2[x]*temp.p +
                           pred_out$b11.2q[x]*temp.p^2 +
                           pred_out$b11.3[x]*mean(pred_data$ele) +
                           pred_out$b11.4[x]*mean(pred_data$head))	  
  pred.predict.temp <- (1/(1+1/(exp(pred.predict.temp))))
  
  # Bind predictions from each sample into one dataframe
  pred.predict.temp.all <- rbind(pred.predict.temp.all,pred.predict.temp)
  
}

### Mean and credible intervals  ###
pred.temp.md <- apply(pred.predict.temp.all, 2, quantile, probs=0.5)
pred.temp.up <- apply(pred.predict.temp.all, 2, quantile, probs=0.975)
pred.temp.lo <- apply(pred.predict.temp.all, 2, quantile, probs=0.025)

### Elevation  #################################################################

### Create vector of values across the range of predictor  ###
ele.p <- seq(min(pred_data$ele,na.rm = T), 
             max(pred_data$ele,na.rm = T),
             length.out = 1000)

### Empty Vector to contain predictions  ##
pred.predict.ele.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(pred_out),1)
  pred.predict.ele <- pred_out$b0.0[x] +
    pred_out$b0.2[x]*mean(pred_data$temp,na.rm = T) + 
    pred_out$b0.6[x]*mean(pred_data$sal,na.rm = T)  +
    pred_out$b0.10[x]*mean(pred_data$sa_height,na.rm = T) +
    pred_out$b0.11[x]*(pred_out$b11.0[x] +
                           pred_out$b11.2[x]*mean(pred_data$temp,na.rm = T) +
                           pred_out$b11.2q[x]*mean(pred_data$temp,na.rm = T)^2 +
                           pred_out$b11.3[x]*ele.p +
                           pred_out$b11.4[x]*mean(pred_data$head))  
  pred.predict.ele <- (1/(1+1/(exp(pred.predict.ele))))
  
  # Bind predictions from each sample into one dataframe
  pred.predict.ele.all <- rbind(pred.predict.ele.all,pred.predict.ele)
  
}

### Mean and credible intervals  ###
pred.ele.md <- apply(pred.predict.ele.all, 2, quantile, probs=0.5)
pred.ele.up <- apply(pred.predict.ele.all, 2, quantile, probs=0.975)
pred.ele.lo <- apply(pred.predict.ele.all, 2, quantile, probs=0.025)

### Head #######################################################################

### Create vector of values across the range of predictor  ###
head.p <- seq(min(pred_data$head,na.rm = T), 
              max(pred_data$head,na.rm = T),
              length.out = 1000)

### Empty Vector to contain predictions  ##
pred.predict.head.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(pred_out),1)
  pred.predict.head <- pred_out$b0.0[x] +
    pred_out$b0.2[x]*mean(pred_data$temp,na.rm = T) + 
    pred_out$b0.6[x]*mean(pred_data$sal,na.rm = T)  +
    pred_out$b0.10[x]*mean(pred_data$sa_height,na.rm = T) +
    pred_out$b0.11[x]*(pred_out$b11.0[x] +
                         pred_out$b11.2[x]*mean(pred_data$temp,na.rm = T) +
                         pred_out$b11.2q[x]*mean(pred_data$temp,na.rm = T)^2 +
                         pred_out$b11.3[x]*mean(pred_data$ele) +
                         pred_out$b11.4[x]*head.p)  
  pred.predict.head <- (1/(1+1/(exp(pred.predict.head))))
  
  # Bind predictions from each sample into one dataframe
  pred.predict.head.all <- rbind(pred.predict.head.all,pred.predict.head)
  
}

### Mean and credible intervals  ###
pred.head.md <- apply(pred.predict.head.all, 2, quantile, probs=0.5)
pred.head.up <- apply(pred.predict.head.all, 2, quantile, probs=0.975)
pred.head.lo <- apply(pred.predict.head.all, 2, quantile, probs=0.025)


### Main  ######################################################################

### Create vector of values across the range of predictor  ###
main.p <- seq(min(pred_data$main,na.rm = T), 
              max(pred_data$main,na.rm = T),
              length.out = 1000)

### Empty Vector to contain predictions  ##
pred.predict.main.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(pred_out),1)
  pred.predict.main <- pred_out$b0.0[x] +
    pred_out$b0.2[x]*mean(pred_data$temp,na.rm = T) + 
    #sal
    pred_out$b0.6[x]*(pred_out$b6.0[x] +
                          pred_out$b6.2[x]*mean(pred_data$temp,na.rm = T) + 
                          pred_out$b6.2q[x]*mean(pred_data$temp,na.rm = T)^2+
                          pred_out$b6.5[x]*main.p)  +
    pred_out$b0.10[x]*(pred_out$b10.0[x] +
                           pred_out$b10.5[x]*main.p) +
    pred_out$b0.11[x]*mean(pred_data$sub_time,na.rm = T)	  
  pred.predict.main <- (1/(1+1/(exp(pred.predict.main))))
  
  # Bind predictions from each sample into one dataframe
  pred.predict.main.all <- rbind(pred.predict.main.all,pred.predict.main)
  
}

### Mean and credible intervals  ###
pred.main.md <- apply(pred.predict.main.all, 2, quantile, probs=0.5)
pred.main.up <- apply(pred.predict.main.all, 2, quantile, probs=0.975)
pred.main.lo <- apply(pred.predict.main.all, 2, quantile, probs=0.025)

# combine as df
pred.main <- as.data.frame(cbind(main.p,pred.main.md,pred.main.up,pred.main.lo))


################################################################################
### Recruitment prediction ###
################################################################################


### Temperature  ###############################################################

### Create vector of values across the range of predictor  ###
temp.r <- seq(min(recruit_data$mean_water_temp,na.rm = T), 
              max(recruit_data$mean_water_temp,na.rm = T),
              length.out = 1000)

### Empty Vector to contain predictions  ##
recruit.predict.temp.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(recruit_out),1)
  recruit.predict.temp <- exp(recruit_out$b0.0[x] + 
                                recruit_out$b0.2[x]*temp.r +
                                recruit_out$b0.5[x]*(recruit_out$b5.0[x] + 
                                                       recruit_out$b5.2[x]*temp.r + 
                                                       recruit_out$b5.2q[x]*temp.r^2+
                                                       recruit_out$b5.14[x]*mean(recruit_data$main,na.rm = T)) + 
                                recruit_out$b0.10[x]*mean(recruit_data$sa_density,na.rm = T) +
                                recruit_out$b0.11[x]*(recruit_out$b11.0[x] + 
                                                        recruit_out$b11.12[x]*mean(recruit_data$ele,na.rm = T) + 
                                                        recruit_out$b11.14[x]*mean(recruit_data$main,na.rm = T) + 
                                                        recruit_out$b11.2[x]*temp.r) + 
                                recruit_out$b0.14[x]*mean(recruit_data$main,na.rm = T) +
                                recruit_out$b0.15[x]*mean(recruit_data$burrow,na.rm = T))
  
  # Bind predictions from each sample into one dataframe
  recruit.predict.temp.all <- rbind(recruit.predict.temp.all,recruit.predict.temp)
  
}

### Mean and credible intervals  ###
recruit.temp.md <- apply(recruit.predict.temp.all, 2, quantile, probs=0.5)
recruit.temp.up <- apply(recruit.predict.temp.all, 2, quantile, probs=0.975)
recruit.temp.lo <- apply(recruit.predict.temp.all, 2, quantile, probs=0.025)

### Elevation  #################################################################

### Create vector of values across the range of predictor  ###
ele.r <- seq(min(recruit_data$ele,na.rm = T), 
             max(recruit_data$ele,na.rm = T),
             length.out = 1000)

### Empty Vector to contain predictions  ##
recruit.predict.ele.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(recruit_out),1)
  recruit.predict.ele  <-  exp(recruit_out$b0.0[x] + 
                                 recruit_out$b0.2[x]*mean(recruit_data$mean_water_temp,na.rm = T) +
                                 recruit_out$b0.5[x]*mean(recruit_data$mean_sal,na.rm = T) + 
                                 recruit_out$b0.10[x]*mean(recruit_data$sa_density,na.rm = T) +
                                 recruit_out$b0.11[x]*(recruit_out$b11.0[x] + 
                                                         recruit_out$b11.12[x]*ele.r + 
                                                         recruit_out$b11.14[x]*mean(recruit_data$main,na.rm = T) + 
                                                         recruit_out$b11.2[x]*mean(recruit_data$mean_water_temp,na.rm = T)) + 
                                 recruit_out$b0.14[x]*mean(recruit_data$main,na.rm = T) +
                                 recruit_out$b0.15[x]*(exp(recruit_out$b15.0 [x]+ 
                                                             recruit_out$b15.12[x]*ele.r + 
                                                             recruit_out$b15.13[x]*mean(recruit_data$head) + 
                                                             recruit_out$b15.14[x]*mean(recruit_data$main)+
                                                             recruit_out$b15.10[x]*mean(recruit_data$sa_density) +
                                                             recruit_out$b15.11[x]*(recruit_out$b11.0[x] + 
                                                                                      recruit_out$b11.12[x]*ele.r + 
                                                                                      recruit_out$b11.14[x]*mean(recruit_data$main,na.rm = T) + 
                                                                                      recruit_out$b11.2[x]*mean(recruit_data$mean_water_temp,na.rm = T))))) 
  
  # Bind predictions from each sample into one dataframe
  recruit.predict.ele.all <- rbind(recruit.predict.ele.all,recruit.predict.ele)
  
}

### Mean and credible intervals  ###
recruit.ele.md <- apply(recruit.predict.ele.all, 2, quantile, probs=0.5)
recruit.ele.up <- apply(recruit.predict.ele.all, 2, quantile, probs=0.975)
recruit.ele.lo <- apply(recruit.predict.ele.all, 2, quantile, probs=0.025)

### Head  ######################################################################

### Create vector of values across the range of predictor  ###
head.r <- seq(min(recruit_data$head,na.rm = T), 
              max(recruit_data$head,na.rm = T),
              length.out = 1000)

### Empty Vector to contain predictions  ##
recruit.predict.head.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(recruit_out),1)
  recruit.predict.head  <- exp(recruit_out$b0.0[x] + 
                                 recruit_out$b0.2[x]*mean(recruit_data$mean_water_temp,na.rm = T) +
                                 recruit_out$b0.5[x]*mean(recruit_data$mean_sal,na.rm = T) + 
                                 recruit_out$b0.10[x]*(exp(recruit_out$b10.0[x] + 
                                                             recruit_out$b10.13[x]*head.r + 
                                                             recruit_out$b10.14[x]*mean(recruit_data$main))) +
                                 recruit_out$b0.11[x]*mean(recruit_data$sa_height,na.rm = T) + 
                                 recruit_out$b0.14[x]*mean(recruit_data$main,na.rm = T) +
                                 recruit_out$b0.15[x]*(exp(recruit_out$b15.0[x] + 
                                                             recruit_out$b15.12[x]*mean(recruit_data$ele) + 
                                                             recruit_out$b15.13[x]*head.r + 
                                                             recruit_out$b15.14[x]*mean(recruit_data$main)+
                                                             recruit_out$b15.10[x]*(exp(recruit_out$b10.0[x] + 
                                                                                          recruit_out$b10.13[x]*head.r + 
                                                                                          recruit_out$b10.14[x]*mean(recruit_data$main))) +
                                                             recruit_out$b15.11[x]*mean(recruit_data$sa_height))))
  # Bind predictions from each sample into one dataframe
  recruit.predict.head.all <- rbind(recruit.predict.head.all,recruit.predict.head)
  
}

### Mean and credible intervals  ###
recruit.head.md <- apply(recruit.predict.head.all, 2, quantile, probs=0.5)
recruit.head.up <- apply(recruit.predict.head.all, 2, quantile, probs=0.975)
recruit.head.lo <- apply(recruit.predict.head.all, 2, quantile, probs=0.025)

### Main #######################################################################

### Create vector of values across the range of predictor  ###
main.r <- seq(min(recruit_data$main,na.rm = T), 
              max(recruit_data$main,na.rm = T),
              length.out = 1000)

### Empty Vector to contain predictions  ##
recruit.predict.main.all <- NULL

### Sample 1000 predictions of posterior distributions  ###
# predictor of interests varies across range of values, all others held constant
for (i in 1:1000)  {
  x <- sample(1:nrow(recruit_out),1)
  recruit.predict.main  <- exp(recruit_out$b0.0[x] + 
                                 recruit_out$b0.2[x]*mean(recruit_data$mean_water_temp,na.rm = T) +
                                 recruit_out$b0.5[x]*(recruit_out$b5.0[x] + 
                                                        recruit_out$b5.2[x]*mean(recruit_data$mean_water_temp,na.rm = T) + 
                                                        recruit_out$b5.2q[x]*mean(recruit_data$mean_water_temp,na.rm = T)^2+
                                                        recruit_out$b5.14[x]*main.r) + 
                                 recruit_out$b0.10[x]*(exp(recruit_out$b10.0[x] + 
                                                             recruit_out$b10.13[x]*mean(recruit_data$head) + 
                                                             recruit_out$b10.14[x]*main.r)) +
                                 recruit_out$b0.11[x]*(recruit_out$b11.0[x] + 
                                                         recruit_out$b11.12[x]*mean(recruit_data$ele,na.rm = T) + 
                                                         recruit_out$b11.14[x]*main.r + 
                                                         recruit_out$b11.2[x]*mean(recruit_data$mean_water_temp,na.rm = T)) + 
                                 recruit_out$b0.14[x]*main.r +
                                 recruit_out$b0.15[x]*(exp(recruit_out$b15.0[x] + 
                                                             recruit_out$b15.12[x]*mean(recruit_data$ele) + 
                                                             recruit_out$b15.13[x]*mean(recruit_data$head) + 
                                                             recruit_out$b15.14[x]*main.r+
                                                             recruit_out$b15.10[x]*(exp(recruit_out$b10.0[x] + 
                                                                                          recruit_out$b10.13[x]*mean(recruit_data$head) + 
                                                                                          recruit_out$b10.14[x]*main.r)) +
                                                             recruit_out$b15.11[x]*(recruit_out$b11.0[x] + 
                                                                                      recruit_out$b11.12[x]*mean(recruit_data$ele,na.rm = T) + 
                                                                                      recruit_out$b11.14[x]*main.r + 
                                                                                      recruit_out$b11.2[x]*mean(recruit_data$mean_water_temp,na.rm = T)))))
  # Bind predictions from each sample into one dataframe
  recruit.predict.main.all <- rbind(recruit.predict.main.all,recruit.predict.main)
  
}

### Mean and credible intervals  ###
recruit.main.md <- apply(recruit.predict.main.all, 2, quantile, probs=0.5)
recruit.main.up <- apply(recruit.predict.main.all, 2, quantile, probs=0.975)
recruit.main.lo <- apply(recruit.predict.main.all, 2, quantile, probs=0.025)


################################################################################
### Plotting features  ###
################################################################################

### y axis labels  ###
y.growth <- "Growth (mm/mo.)"
y.recruit <- "Recruits"
y.pred <- "Prop. depredated"

### x axis labels  ###
x.temp.g <- "Maximum mussel temperature (C)"
x.temp.pr <- "Water temperature (C)"
x.temp <- "Temperature (c)"
x.ele <- "Elevation (m)"
x.head <- "Distance (m) to head"
x.main <- "Distance (m) to Subtidal"

### Colors ###
temp.col.d <- "#5E163A"
temp.col.l <- "#882255"

ele.col.d <- "#A76E00"
ele.col.l <- "#E69F00"

main.col.d <- "#2E7A6B"
main.col.l <- "#44AA99"

head.col.d <- "#762F6E"
head.col.l <- "#AA4499"

### y axis limits ###
y.lim.g <- c(0,3)
y.lim.r <- c(0,30)
y.lim.p <- c(0,0.4)

### Sizes ###
size <- 150
point.size <- 5

### Exporting  ###
dpi <- 600
dim.x <- 2000
dim.y <- 1600


################################################################################
### Create plots  ###
################################################################################

scatter.size = 12

### Growth  ####################################################################

### Temp ###
growth.temp.plot <- ggplot(data = pred.main, aes(temp,growth.temp.md*30))+
  geom_smooth(size = point.size, col = temp.col.d)+
  geom_ribbon(aes(ymin = growth.temp.lo*30,ymax =growth.temp.up*30),
              alpha=0.2,fill = temp.col.l)+
  geom_point(data = growth_data, aes(x =mean_water_temp, y=growth_rate*30),
             size = scatter.size,col=temp.col.d) +
  xlab("") +
  ylab(y.growth) +
  scale_y_continuous(limits = y.lim.g,
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Elevation ###
growth.ele.plot <- ggplot(data = pred.main, aes(ele,growth.ele.md*30))+
  geom_smooth(size = point.size, col = ele.col.d)+
  geom_ribbon(aes(ymin = growth.ele.lo*30,ymax =growth.ele.up*30),
              alpha=0.2,fill = ele.col.l)+
  geom_point(data = growth_data, aes(x =ele, y=growth_rate*30),
             size = scatter.size,col=ele.col.d) +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = y.lim.g,
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Head ###
growth.head.plot <- ggplot(data = pred.main, aes(head,growth.head.md*30))+
  geom_smooth(size = point.size, col = head.col.d)+
  geom_ribbon(aes(ymin = growth.head.lo*30,ymax =growth.head.up*30),
              alpha=0.2, fill=head.col.l)+
  geom_point(data = growth_data, aes(x =head, y=growth_rate*30),
             size = scatter.size,col=head.col.d) +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = y.lim.g,
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Main ###
growth.main.plot <- ggplot(data = pred.main, aes(main,growth.main.md*30))+
  geom_smooth(size = point.size, col = main.col.d)+
  geom_ribbon(aes(ymin = growth.main.lo*30,ymax =growth.main.up*30),alpha=0.2,
              fill = main.col.l)+
  geom_point(data = growth_data, aes(x =main, y=growth_rate*30),
             size = scatter.size,col=main.col.d) +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = y.lim.g,
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Predation  #################################################################

### Temp ###
pred.temp.plot <- ggplot(data = pred.main, aes(temp.p,pred.temp.md))+
  geom_smooth(size = point.size, col = temp.col.d)+
  geom_ribbon(aes(ymin = pred.temp.lo,ymax =pred.temp.up),
              alpha=0.2, fill = temp.col.l)+
  geom_point(data = pred_data, aes(x =temp, y=prop_pred),
             size = scatter.size,col=temp.col.d) +
  xlab(x.temp) +
  ylab(y.pred) +
  scale_y_continuous(limits = c(0,1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Elevation  ###
pred.ele.plot <- ggplot(data = pred.main, aes(ele.p,pred.ele.md))+
  geom_smooth(size = point.size, col = ele.col.d)+
  geom_ribbon(aes(ymin = pred.ele.lo,ymax =pred.ele.up),
              alpha=0.2, fill = ele.col.l)+
  geom_point(data = pred_data, aes(x =ele, y=prop_pred),
             size = scatter.size,col=ele.col.d) +
  xlab(x.ele) +
  ylab("") +
  scale_y_continuous(limits = c(0,1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Head ###
pred.head.plot <- ggplot(data = pred.main, aes(head.p,pred.head.md))+
  geom_smooth(size = point.size, col = head.col.d)+
  geom_ribbon(aes(ymin = pred.head.lo,ymax =pred.head.up),
              alpha=0.2, fill = head.col.d)+
  geom_point(data = pred_data, aes(x =head, y=prop_pred),
             size = scatter.size,col=head.col.d) +
  xlab(x.head) +
  ylab("") +
  scale_y_continuous(limits = c(0,1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Main ###
pred.main.plot <- ggplot(data = pred.main, aes(main.p,pred.main.md))+
  geom_smooth(size = point.size, col = main.col.d)+
  geom_ribbon(aes(ymin = pred.main.lo,ymax =pred.main.up),
              alpha=0.2,fill = main.col.l)+
  geom_point(data = pred_data, aes(x =main, y=prop_pred),
             size = scatter.size,col=main.col.d) +
  xlab(x.main) +
  ylab(y.pred) +
  ylab("") +
  scale_y_continuous(limits = c(0,1),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(limits = c(0,175))+
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Recruitment  ###############################################################

### Temp ###
recruit.temp.plot <- ggplot(data = pred.main, aes(temp.r,recruit.temp.md))+
  geom_smooth(size = point.size, col = temp.col.d)+
  geom_ribbon(aes(ymin = recruit.temp.lo,ymax =recruit.temp.up),
              alpha=0.2,fill = temp.col.l)+
  geom_point(data = recruit_data, aes(x =mean_water_temp, y=recruit),
             size = scatter.size,col=temp.col.d) +
  xlab("") +
  ylab(y.recruit) +
  scale_y_continuous(limits = c(0,55)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Elevation ###
recruit.ele.plot <- ggplot(data = pred.main, aes(ele.r,recruit.ele.md))+
  geom_smooth(size = point.size, col = ele.col.d)+
  geom_ribbon(aes(ymin = recruit.ele.lo,ymax =recruit.ele.up),
              alpha=0.2,fill = ele.col.l)+
  geom_point(data = recruit_data, aes(x =ele, y=recruit),
             size = scatter.size,col=ele.col.d) +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(0,55)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Head ###
recruit.head.plot <- ggplot(data = pred.main, aes(head.r,recruit.head.md))+
  geom_smooth(size = point.size, col = head.col.d)+
  geom_ribbon(aes(ymin = recruit.head.lo,ymax =recruit.head.up),
              alpha=0.2, fill = head.col.l)+
  geom_point(data = recruit_data, aes(x =head, y=recruit),
             size = scatter.size,col=head.col.d) +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(0,55)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Main ###
recruit.main.plot <- ggplot(data = pred.main, aes(main.r,recruit.main.md))+
  geom_smooth(size = point.size, col = main.col.d)+
  geom_ribbon(aes(ymin = recruit.main.lo,ymax =recruit.main.up),
              alpha=0.2,fill = main.col.l)+
  geom_point(data = recruit_data, aes(x =main, y=recruit),
             size = scatter.size,col=main.col.d) +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(0,55)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))


################################################################################
### Export plots  ###
################################################################################

### Growth  ####################################################################

### Temp ###
png("plots/total_effect/growth_temp_plot_v2.png",
    width = dim.x,
    height = dim.y,
    pointsize = dpi)
print(growth.temp.plot)
dev.off()

### Elevation ###
png("plots/total_effect/growth_ele_plot_v2.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(growth.ele.plot)
dev.off()

### Head ###
png("plots/total_effect/growth_head_plot_v2.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(growth.head.plot)
dev.off()

### Main ###
png("plots/total_effect/growth_main_plot_v2.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(growth.main.plot)
dev.off()

### Predation ##################################################################

### Temp  ###
png("plots/total_effect/pred_temp_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(pred.temp.plot)
dev.off()

### Elevation  ###
png("plots/total_effect/pred_ele_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(pred.ele.plot)
dev.off()

### Head  ###
png("plots/total_effect/pred_head_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(pred.head.plot)
dev.off()

### Main ###
png("plots/total_effect/pred_main_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(pred.main.plot)
dev.off()

### Recruitment ################################################################

### Temp  ###
png("plots/total_effect/recruit_temp_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(recruit.temp.plot)
dev.off()

### Elevation  ###
png("plots/total_effect/recruit_ele_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(recruit.ele.plot)
dev.off()

### Head  ###
png("plots/total_effect/recruit_head_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(recruit.head.plot)
dev.off()

### Main ###
png("plots/total_effect/recruit_main_plot.png",
    width = dim.x, 
    height = dim.y,
    pointsize = dpi)
print(recruit.main.plot)
dev.off()


