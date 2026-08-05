################################################################################
# 
#  Summary plots and tables                            
#
#  by Anonymized                                  
#
################################################################################

# This script creates summary plots for the population vital rate experiments.
# For the growth experiment, growth rates, condition index, and mortality plots
# are created. For the predation experiment box plots for the total number of 
# mussels predated during a trial are created. For the recruitment experiment,
# the raw monthly mussel recruit counts are plotted. In addition, the mean
# and standard deviation for the predictor variables are calculated at the
# trial level for each vital rate and are exported as tables.


################################################################################
##  House Keeping  ##
################################################################################


### Clear workspace  ###
rm(list = ls())

### Load in packages  ###
library(ggplot2)
library(dplyr)
library(readxl)


################################################################################
##  Load in data  ##
################################################################################

### Vital rate and predictor data  ###
growth <-read_excel("mussel_population_vital_data.xlsx",sheet = 1)
predate <- read_excel("mussel_population_vital_data.xlsx",sheet = 2)
recruit <-read_excel("mussel_population_vital_data.xlsx",sheet = 3)

### Monthly recruitment data  ###
recruit_mon<-read_excel("mussel_population_vital_data.xlsx",sheet = 6)  #per day

### Plot parameter  ###
size = 130

################################################################################
### Growth experiment plots  ###
################################################################################

### Format data  ###############################################################
# change season names to sentence case
substr(growth$season,1,1) <- toupper(substr(growth$season,1,1))

# Change season to a factor order chronologically
growth$season <- factor(growth$season, levels = c("Spring","Summer", "Fall"))

# CHange growth rate to monthly scale for better visualization
growth$growth_rate_mon <-growth$growth_rate*30

### Create vital rate plots  ###################################################

### Growth  ###
g_plot <- ggplot(growth, aes(x = season, y = growth_rate_mon))+
  #geom_violin()
  geom_boxplot(fill ="grey",outlier.size=16,lwd =5) +
  xlab("Trial") +
  ylab("Growth rate (mm/day)") +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Condition  ###
c_plot <- ggplot(growth, aes(x = season, y = condition))+
  #geom_violin()
  geom_boxplot(fill ="grey",outlier.size=16,lwd =5) +
  xlab("Trial") +
  ylab("Mussel condition") +
  #scale_y_continuous(breaks = seq(0, 15, by = 3)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Mortality ###
m_plot <- ggplot(
  growth %>% mutate(prop_mort = mortality/n), 
  aes(x = season, y = prop_mort)
  )+
  geom_boxplot(fill ="grey", outlier.size=16,lwd =5) +
  xlab("Trial") +
  ylab("Proportion desiccated") +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Export plots  ##############################################################

### Growth  ###
png("plots/summary/growth_box.png",width = 2000, height = 1600,
    pointsize = 600)
print(g_plot)
dev.off()

### Condition  ###
png("plots/summary/cond_box.png",width = 2000, height = 1600,
    pointsize = 600)
print(c_plot)
dev.off()

### Mortality  ###
png("plots/summary/mort_box.png",width = 2000, height = 1600,
    pointsize = 600)
print(m_plot)
dev.off()


################################################################################
### Predation experiment ###
###############################################################################

### Format data  ###############################################################

# Create a vector for trial
pred_trial <- 1:5

# name each vector after the month
names(pred_trial) <- c("Jun", "Jul", "Aug", "Oct", "Dec")

# Match trial to main data to replace experiment
predate$trial <- names(pred_trial)[match(predate$experiment,pred_trial)]

# CHange trial into factor
predate$trial <- factor(predate$trial, levels = names(pred_trial))

### Create plot  ###############################################################
p_plot <- ggplot(predate, aes(x = trial, y = prop_pred))+
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),size = 1.5) +
  geom_boxplot(fill = "grey", outlier.size=16,lwd =5) + ##FEE79E
  xlab("Trial") +
  ylab("Proportion depredated") +
  scale_y_continuous(limits = c(0,1),
                     labels = scales::number_format(accuracy = 0.1))+
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))


### Export  ####################################################################
png("plots/summary/pred_box.png",width = 2000, height = 1600,
    pointsize = 600)
print(p_plot)
dev.off()


################################################################################
### Recruitment experiment plots ###
################################################################################

### Format data ###############################################################

# CHange month to factor
r.month = c("Feb", "Mar", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec")
recruit_mon$month <- factor(recruit_mon$month, levels = r.month)
recruit_mon$recruit_day <- as.numeric(recruit_mon$recruit_day)

### Creat plots  ###############################################################
rm_plot <- ggplot(recruit_mon, aes(x = month, y = recruit_day))+
  geom_boxplot(fill ="grey", outlier.size=16,lwd =5) +
  xlab("Month") +
  ylab("Recruits per day") +
  #scale_y_continuous(breaks = seq(0, 60, by = 10)) +
  theme_classic(base_size = size) +
  theme(legend.position = "none",
        panel.border =  element_rect(color = "black", fill = NA, size = 5))

### Export plot  ###############################################################
png("plots/summary/recruit_mon_box.png",width = 4000, height = 1600,
    pointsize = 600)
print(rm_plot)
dev.off()


################################################################################
### Population predictor summaries  ###
################################################################################
# Find the trial mean and standard deviation of each vital rates predictors. 
# These are used of table S1 in Appendix S6

### Format data ################################################################

### Growth ###
growth_format <- growth %>% 
  mutate(season = case_when(
    # change season names to trials used in manuscript
    season == "Spring" ~ "Trial 1",
    season == "Summer" ~ "Trial 2",
    season == "Fall" ~ "Trial 3",
    TRUE ~ season
  )) %>% 
  # retain only columns used in hyp sem
  select(season, mussel_size,ele,head,main,mean_water_temp,
         mean_sal,mean_bact,mean_chla,mean_om,sub_prop,flow_rate,
         gd_density_max_sum,sa_density,sa_height)

### Predation ###

pred_format <- predate %>% 
  mutate(season = paste("Trial",experiment))%>% 
  
  # retain only columns used in hyp sem
  select(season,mean_size_pred_exp,ele,head,main,temp,
         sal,sub_time,burrow,
         gd_density_max_sum,sa_density,sa_height)

### Recruitment
recruit_format <- recruit %>% 
  mutate(season = case_when(
    # change season names to trials used in manuscript
    season == "spring" ~ "Trial 1",
    season == "summer" ~ "Trial 2",
    season == "fall" ~ "Trial 3",
    TRUE ~ season
  )) %>% 
  
  # retain only columns used in hyp sem
  select(season,ele,head,main,mean_water_temp,
         mean_sal,sub_prop,flow_rate,burrow,
         gd_density_max_sum,sa_density,sa_height)

### Summarize data #############################################################

# Function to summarize each variable at each trial. Produces long format
# data frame with trial as column, and variable as row
sumFun <- function(data) {
  
  # Function to calculate mean and sd, and format into one columns i.e., MEAN (SD)
  mean.sd.fun <- function(x) {
    mean = round(mean(x,na.rm = T),2)
    sd = round(sd(x, na.rm = T),2)
    mean_sd = paste(mean," (",sd,")", sep = "")
    return(mean_sd)
  }  # end of mean.sd function
  
  # calculate the mean and sd of each var at each trial  
  sum <- data %>% 
    group_by(season) %>% 
    summarise(across(.cols = where(is.numeric),
                     .fns = mean.sd.fun))
  
  # Transpose and format
  sum.t <- t(sum)  
  colnames(sum.t) <- sum.t[1,]
  sum.t <- sum.t[-1,]
  sum.t <- as.data.frame(sum.t)
  
  # Return transposed summaries
  return(sum.t)
}  # end of sumFun

### Calculate for each vital rate  #############################################
growth_sum <- sumFun(growth_format)
recruit_sum <- sumFun(recruit_format)
pred_sum <- sumFun(pred_format)

### export #####################################################################
write.csv(growth_sum,"plots/summary/growth_driver_means.csv", row.names = T)
write.csv(recruit_sum,"plots/summary/recruit_driver_means.csv", row.names = T)
write.csv(pred_sum,"plots/summary/pred_driver_means.csv", row.names = T)
