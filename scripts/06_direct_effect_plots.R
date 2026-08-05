################################################################################
# 
#  Direct effect plots for all SEM pathways            
#
#  by Anonymized                                
#
################################################################################

# This script creates predictions for every response variable in each structural
# equation models. These are used to create direct effect plots for each 
# SEM pathway. Each plot is color coded based on predictor variable. A function
# is created in this script that creates the predictions and plots for every
# pathway, and then exports each to a folder. This plots were used to create
# figure 3 of the main text using adobe illustrated and the standard effect
# queries calculated in the SEM scripts.


################################################################################
##  House Keeping  ##
################################################################################


### Clear workspace  ###
rm(list = ls())

### Load in packages  ###
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)


################################################################################
##  Load in data  ##
################################################################################


### growth ###
growth_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 1) %>% 
  rename(sub = sub_prop,
         temp = mean_water_temp,
         sal = mean_sal) %>% 
  select(growth_rate,mussel_size,sub,flow_rate,sal,temp,ele,main,head)

growth_out <- readRDS("model_out/growth_SEM_final_output.rds")
growth_out <- as.data.frame(growth_out$sims.list)

### Recruit ###
recruit_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 3) %>% 
  rename(temp = mean_water_temp,
         sal = mean_sal) %>% 
  select(recruit,sa_height,sa_density,burrow,sal,temp,ele,main,head)
recruit_out <- readRDS("model_out/recruitment_SEM_final_output.rds")
recruit_out <- as.data.frame(recruit_out$sims.list)

### Predation ###
pred_data <- read_excel("mussel_population_vital_data.xlsx",sheet = 2) %>% 
  rename(sub = sub_time) %>% 
  select(prop_pred,sub,sa_height,temp,sal,ele,main,head)
pred_out <- readRDS("model_out/predation_SEM_final_output.rds")
pred_out <- as.data.frame(pred_out$sims.list)


################################################################################
### Prep data ####
################################################################################

### Growth  ####################################################################

# Load in variable code key
growth_codes <- read.csv("plot_parameters/growth_codes.csv")

# Extract all combination of response and predictor codes
g_out_codes <- colnames(growth_out)[substr(colnames(growth_out),1,1) == "b"]

# Create data frame
g_out_codes_df <- data.frame(y_code = as.numeric(substr(g_out_codes,2,2)),
                       x_code = as.numeric(substr(g_out_codes,4,5)))

## Merge in variable names
growth_code_df <- g_out_codes_df %>% 
  # Remove na and zero codes
  filter(!is.na(x_code),  
         x_code !=0) %>% 
  
  # Merge in response codes
  left_join(growth_codes, by =join_by(y_code == code)) %>% 
  rename(var_y = var) %>% 
  
  # Merge in predictor codes
  left_join(growth_codes, by =join_by(x_code == code)) %>% 
  rename(var_x = var) %>% 
  
  # Add vital code
  mutate(vital = "g") %>% 
  select(-y_code,-x_code)

### Predation ##################################################################

# Load in variable code key
pred_codes <- read.csv("plot_parameters/pred_codes.csv")

# Extract all combination of response and predictor codes
p_out_codes <- substr(colnames(pred_out)[substr(colnames(pred_out),1,1) == "b"],
                      2,6)

# Create data frame
p_out_codes_df <- data.frame(y_code = as.numeric(str_extract(p_out_codes, "^[^\\.]+")),
                             x_code = as.numeric(str_extract(p_out_codes, "(?<=\\.)[^\\.]*$")))



## Merge in variable names
pred_code_df <- p_out_codes_df %>% 
  # Remove na and zero codes
  filter(!is.na(x_code),  
         x_code !=0) %>% 
  
  # Merge in response codes
  left_join(pred_codes, by =join_by(y_code == code)) %>% 
  rename(var_y = var) %>% 
  
  # Merge in predictor codes
  left_join(pred_codes, by =join_by(x_code == code)) %>% 
  rename(var_x = var) %>% 
  
  # Add vital code
  mutate(vital = "p") %>% 
  select(-y_code,-x_code)

### Recruitment ################################################################

# Load in variable code key
recruit_codes <- read.csv("plot_parameters/recruit_codes.csv")

# Extract all combination of response and predictor codes
r_out_codes <- substr(colnames(recruit_out)[substr(colnames(recruit_out),1,1) == "b"],
                      2,6)

# Create data frame
r_out_codes_df <- data.frame(y_code = as.numeric(str_extract(r_out_codes, "^[^\\.]+")),
                             x_code = as.numeric(str_extract(r_out_codes, "(?<=\\.)[^\\.]*$")))



## Merge in variable names
recruit_code_df <- r_out_codes_df %>% 
  # Remove na and zero codes
  filter(!is.na(x_code),  
         x_code !=0) %>% 
  
  # Merge in response codes
  left_join(recruit_codes, by =join_by(y_code == code)) %>% 
  rename(var_y = var) %>% 
  
  # Merge in predictor codes
  left_join(recruit_codes, by =join_by(x_code == code)) %>% 
  rename(var_x = var) %>% 
  
  # Add vital code
  mutate(vital = "r") %>% 
  select(-y_code,-x_code)

## Merge all data into one df ##################################################
vital_df <- do.call(bind_rows,list(growth_code_df,recruit_code_df,pred_code_df))

# Merge in color codes
color_codes <- read.csv("plot_parameters/color_palette.csv")
vital_df <- vital_df %>% 
  left_join(color_codes)

# Merge in manually selected y axis limits
ymax_df <- read.csv("plot_parameters/ymax.csv")  # ylimits for some variables

vital_df <- vital_df %>% 
  left_join(ymax_df) #%>% 
  # filter(!is.na(ymax))

################################################################################
### Create plotting function ###
################################################################################

# This function creates predictions for each combination of predictor and 
# response variables given. Data is extracted fro the plotting df that selects
# the x and y variables of interest, which vital rate SEM they come from, 
# y-limits for plots, and color palette of the predictor variables. Predictions  
# are estimated across the range of the predictor variable of interest, with all 
# others held constant. 

# Plots are created using these predictions with parameters defined by the user.
# Plots are colore coded based on predictor of interest and then are exported
# to a directory folder.

directPlotter <- function(df, 
                           point.size = 8, 
                           line.size = 5, 
                           text.size =60, 
                           rect.size = 5,
                           dpi = 600,
                           width = 10,
                           height = 10) {
  
  ##################
  ### Format data ##############################################################
  ##################
  
  df = as.data.frame(df)
  
  ## Extract necessary info from data frame ##
  var_x = df$var_x
  var_y = df$var_y
  vital = df$vital
  dark = df$col_dark
  light = df$col_light
  
  ymin = df$ymin
  ymax = df$ymax
  
  ## Select data for specified vital rate ##
  if (vital == "g") {
    data = growth_data
    out = growth_out
  }
  
  if (vital == "r") {
    data = recruit_data
    out = recruit_out
  }
  
  if (vital == "p") {
    data = pred_data
    out = pred_out
  }
  
  # Create list of mean values for drivers
  xvar_list <- as.list(colMeans(data,na.rm = T))
  
  # replace the mean with range of values of predictor variable
  xvar_list[[var_x]] <- seq(min(data[var_x],na.rm = T), 
                            max(data[var_x],na.rm = T),
                            length.out = 1000)
  
  ################################################
  ### Calculate response variable predictions ##################################
  ################################################
  
  # response variable-specific prediction equations for each population vitals
  
  ### Growth rate ##############################################################
  
  ## growth ##
  if (vital == "g" & var_y == "growth_rate"){
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b0.0[x] + 
        out$b0.2[x]*xvar_list[["temp"]] +
        out$b0.3[x]*xvar_list[["sub"]] + 
        out$b0.4[x]*xvar_list[["flow_rate"]] +
        out$b0.5[x]*xvar_list[["sal"]] + 
        out$b0.9[x]*xvar_list[["mussel_size"]] 
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
      
    }  # end loop
  }  # end if statement
  
  ## submergance ##
  if (vital == "g" & var_y == "sub"){
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b3.0[x] + 
        out$b3.12[x]*xvar_list[["ele"]] +
        out$b3.13[x]*xvar_list[["head"]]  
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
    }  # end loop
  }  # end if statement
  
  ## Flow rate ##
  if (vital == "g" & var_y == "flow_rate"){
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b4.0[x] + 
        out$b4.3[x]*xvar_list[["sub"]] +
        out$b4.12[x]*xvar_list[["ele"]] +
        out$b4.13[x]*xvar_list[["head"]] +
        out$b4.14[x]*xvar_list[["main"]] 
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
    }  # end loop
  }  #end if statement
  
  ## Salinity ##
  if (vital == "g" & var_y == "sal"){
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b5.0[x] + 
        out$b5.2[x]*xvar_list[["temp"]] + out$b5.2q[x]*xvar_list[["temp"]]^2 +
        out$b5.14[x]*xvar_list[["main"]] 
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
    }  # end loop
  }  #end if statement
  
  ### Recruitment ##############################################################
  
  ## recruitment ## 
  if(vital == "r" & var_y == "recruit") {
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- exp(out$b0.0[x] + 
                              out$b0.2[x]*xvar_list[["temp"]] +
                              out$b0.5[x]*xvar_list[["sal"]] +
                              out$b0.10[x]*xvar_list[["sa_density"]] +
                              out$b0.11[x]*xvar_list[["sa_height"]] +
                              out$b0.14[x]*xvar_list[["main"]] +
                              out$b0.15[x]*xvar_list[["burrow"]])
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
      
    }  # end loop
  }  # end if statement
  
  ## burrow ## 
  if(vital == "r" & var_y == "burrow") {
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- exp(out$b15.0[x] + 
                              out$b15.12[x]*xvar_list[["ele"]] +
                              out$b15.13[x]*xvar_list[["head"]] +
                              out$b15.14[x]*xvar_list[["main"]] +
                              out$b15.10[x]*xvar_list[["sa_density"]] +
                              out$b15.11[x]*xvar_list[["sa_height"]])
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
      
    }  # end loop
  }  # end if statement
  
  ## sa_height ## 
  if(vital == "r" & var_y == "sa_height") {
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b11.0[x] + 
        out$b11.12[x]*xvar_list[["ele"]] +
        out$b11.14[x]*xvar_list[["main"]] +
        out$b11.2[x]*xvar_list[["temp"]]
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
      
    }  # end loop
  }  # end if statement
  
  ## sa_density ## 
  if(vital == "r" & var_y == "sa_density") {
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- exp(out$b10.0[x] + 
                              out$b10.13[x]*xvar_list[["head"]] +
                              out$b10.14[x]*xvar_list[["main"]] )
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
      
    }  # end loop
  }  # end if statement
  
  ## salinity ## 
  if(vital == "r" & var_y == "sal") {
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b5.0[x] + 
        out$b5.14[x]*xvar_list[["main"]] +
        out$b5.2[x]*xvar_list[["temp"]] + out$b5.2q[x]*xvar_list[["temp"]]^2
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
      
    }  # end loop
  }  # end if statement
  
  ### Predation ################################################################
  
  ## prop_pred ##
  if(vital == "p" & var_y == "prop_pred"){  
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b0.0[x]+
        out$b0.2[x]*xvar_list[["temp"]] +
        out$b0.6[x]*xvar_list[["sal"]] +
        out$b0.10[x]*xvar_list[["sa_height"]] +
        out$b0.11[x]*xvar_list[["sub"]]	  
      predict.vector <- (1/(1+1/(exp(predict.vector))))
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
      
    }
  }
  
  ## sub ##
  if(vital == "p" & var_y == "sub"){  
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b11.0[x]+
        out$b11.2[x]*xvar_list[["temp"]] + out$b11.2q[x]*xvar_list[["temp"]]^2 +
        out$b11.3[x]*xvar_list[["ele"]] +
        out$b11.4[x]*xvar_list[["head"]]
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
    }
  }
  
  ## sa_height ##
  if(vital == "p" & var_y == "sa_height"){  
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b10.0[x]+
        out$b10.5[x]*xvar_list[["main"]] 
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
    }
  }
  
  ## sal ##
  if(vital == "p" & var_y == "sal"){ 
    predict.vector.all <- NULL
    for (i in 1:1000)  {
      x <- sample(1:nrow(out),1)
      predict.vector <- out$b6.0[x]+
        out$b6.2[x]*xvar_list[["temp"]] + out$b6.2q[x]*xvar_list[["temp"]]^2 +
        out$b6.5[x]*xvar_list[["main"]]
      
      predict.vector.all <- rbind(predict.vector.all,predict.vector)
    }
  }
  
  
  ### Calculate mean and credible intervals of predictions #####
  predict.md <- apply(predict.vector.all, 2, quantile, probs=0.5)
  predict.up <- apply(predict.vector.all, 2, quantile, probs=0.975)
  predict.lo <- apply(predict.vector.all, 2, quantile, probs=0.025)
  
  predict.df <- data.frame(x_pred = xvar_list[[var_x]],
                           y_pred = predict.md,
                           up = predict.up,
                           lo = predict.lo)
  
  #####################
  ### Create plots #############################################################
  #####################
  
  # data for plotting
  plot_data <- data[,c(var_y,var_x)]
  colnames(plot_data) <- c("y","x")
  
  
  ### Set plotting parameters ##################################################
  
  # For growth rate, alter the scale to growth per month, not days
  if(var_y == "growth_rate"){
    plot_data$y <- plot_data$y*30
    predict.df$y_pred <- predict.df$y_pred*30
    predict.df$up <- predict.df$up*30
    predict.df$lo <- predict.df$lo*30
  }

  ### Export directory ###
  folder.path <- "plots/direct/"
  dir <- paste(folder.path,vital,"/",var_y,"_",var_x,".png",sep = "")
 
  ### Create plot ##############################################################
  
  # If specific y limits are specified, use these in plotting
   if(!is.na(ymax)){
     plot_scatter <- ggplot(data = predict.df,aes(x_pred,y_pred)) +
       geom_point(data = plot_data, aes(x,y), color = dark,
                  size = point.size)+
       geom_smooth(data = predict.df,aes(x = x_pred,y=y_pred),
                   color = dark, size = line.size) +
       geom_ribbon(aes(ymin = lo,ymax =up), alpha=0.2,fill = light) +
       scale_y_continuous(limits = c(ymin,ymax)) +
       xlab("") +
       ylab("")+
       theme_classic(base_size = text.size)+
       theme(panel.border =  element_rect(color = "black", fill = NA, linewidth = rect.size))
     
     # otherwise use defult values
   } else {
     plot_scatter <- ggplot(data = predict.df,aes(x_pred,y_pred)) +
       geom_point(data = plot_data, aes(x,y), color = dark,
                  size = point.size)+
       geom_smooth(data = predict.df,aes(x = x_pred,y=y_pred),
                   color = dark, size = line.size) +
       geom_ribbon(aes(ymin = lo,ymax =up), alpha=0.2,fill = light) +
       xlab("") +
       ylab("")+
       theme_classic(base_size = text.size)+
       theme(panel.border =  element_rect(color = "black", fill = NA, linewidth = rect.size))
   }

  ###################
  ## Export plot  ##############################################################
  ###################
  
  ggsave(dir,
         plot = plot_scatter,
         dpi = dpi,
         width = width,
         height = height)
}  # end of function

################################################################################
### Loop function through all variables  ###
################################################################################


for (i in 1:nrow(vital_df)) {
  directPlotter(vital_df[i,])
}  # end of loop


