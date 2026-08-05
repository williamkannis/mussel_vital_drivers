################################################################################
#
#  RIBBED MUSSEL PREDATION STRUCTURAL EQUATION MODELS     
#  
#  Created by Anonymized                               
#  Modified from Grace et al., 2012 - Ecosphere           
#
################################################################################

# This code performs the structural equation modelling (SEM) of the hierarchical
# drivers of ribbed mussel depredation rates across the salt marsh landscape. 
# Here, we create a global Bayesian piecewise SEM based on prior knowledge from 
# literature. We then perform model reparameterization  based on the methods of
# Grace et al. 2012. In this manner, nonsignificant pathways are pruned and 
# missing pathways are added based on relationship among residuals and prior
# ecological knowledge. Finally, standardized queries for direct effects and 
# total effects are calculated to compare effect sizes.

################################################################################
### HOUSE KEEPING  ###
################################################################################

# Clear work space
rm(list = ls())

# Load in packages
library(dplyr)
library(rjags)
library(jagsUI)
library(rstatix)
library(readxl)

################################################################################
### LOAD DATA  ###
################################################################################

## Predation and explanatory data  ##
response <- read_excel("mussel_population_vital_data.xlsx",sheet = 2)

## Check for correlation among explanatory variables  ##
response_cov <- response[,3:20]
response_TF_cov<-abs(cor(response_cov, use = "complete.obs"))>.60

## Check for missing data ##

# what columns have missing data?
colnames(response)[colSums(is.na(response)) > 0]

# Number of missing data per column
length(response$sub_time[is.na(response$sub_time)])
length(response$burrow[is.na(response$burrow)])
length(response$sa_density[is.na(response$sa_height)])
length(response$sa_height[is.na(response$sa_height)])
length(response$temp[is.na(response$temp)])
length(response$sal[is.na(response$sal)])

################################################################################
### VARIABLE CODES FOR MODELS ###
################################################################################

#1. mussel_size
#2. water_temp
#3. ele
#4. head
#5. main
#6. sal
#7. gd_density
#8. burrow
#9. sa_density
#10. sa_height
#11. sub

################################################################################
### Global SEM MODEL  ###
################################################################################

sink("pred.global.model.txt")		
cat("
model {

### LIKELIHOODS  ###############################################################
	for(i in 1:N) {

## Predation ##
#0.Pred	 
	  pred[i] ~ dbin(p[i],n[i])                           # 
		logit(p[i]) <- b0.0 +
	  	b0.1*mussel_size[i] +
	  	b0.2*water_temp[i] + 
	  	b0.6*sal[i] +
	  	b0.7*gd_density[i] + 
	  	b0.8*burrow[i] + 
	    b0.9*sa_density[i] + 
	  	b0.10*sa_height[i] +
	  	b0.11*sub_prop[i]  

## Experiment level variables ##
#2. Water temp
water_temp[i] ~dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0

#6. Salinity
  sal[i] ~ dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <- b6.0 +
    b6.2*water_temp[i] + b6.2q*water_temp[i]*water_temp[i] + 
    b6.5*main[i] 

#7.Mussel density
  gd_density[i] ~ dpois(gd.hat[i])
  log(gd.hat[i]) <- b7.0 +
    b7.3*ele[i]+
    b7.4*head[i]+
    b7.5*main[i] +
    b7.9*sa_density[i]

#8. burrows
burrow[i] ~ dpois(burrow.hat[i])
  log(burrow.hat[i]) <- b8.0 +   
    b8.3*ele[i]+
    b8.4*head[i]+
   b8.5*main[i]

#9. Sa_denisty
  sa_density[i] ~ dpois(sad.hat[i])
  log(sad.hat[i]) <- b9.0 +
    b9.3*ele[i]+
    b9.4*head[i]

#10. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b10.0 +  
    b10.2*water_temp[i] +
    b10.3*ele[i]+
    b10.4*head[i]+
    b10.5*main[i] 
     
#11. #Submergence Time
  sub_prop[i] ~ dnorm(sub.hat[i], tau.sub)
  sub.hat[i] <- b11.0 + 
   b11.3*ele[i] + 
   b11.4*head[i]  

}

### PRECISION VARIABLES  #######################################################
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.height <- 1/(sigma.height*sigma.height)
  tau.temp <- 1/(sigma.temp*sigma.temp)

### PRIORS #####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.1 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001) 
	b0.6 ~ dnorm(0,0.00001); b0.7 ~ dnorm(0,0.00001); b0.8 ~ dnorm(0,0.00001) 
	b0.9 ~ dnorm(0,0.00001); b0.10 ~ dnorm(0,0.00001); b0.11 ~ dnorm(0,0.00001)

	b2.0~ dnorm(0,0.00001); sigma.temp ~ dunif(0,100)
	
	b6.0 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); b6.2q ~ dnorm(0,0.00001);
	b6.5 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)

	b7.0 ~ dnorm(0,0.00001); b7.3 ~ dnorm(0,0.00001); b7.4 ~ dnorm(0,0.00001) 
	b7.5 ~ dnorm(0,0.00001);b7.9 ~ dnorm(0,0.00001)
	
  b8.0 ~ dnorm(0,0.00001); b8.3 ~ dnorm(0,0.00001); b8.4 ~ dnorm(0,0.00001) 
	b8.5 ~ dnorm(0,0.00001)
	
	b9.0 ~ dnorm(0,0.00001); b9.3 ~ dnorm(0,0.00001); b9.4 ~ dnorm(0,0.00001) 
	
	b10.0 ~ dnorm(0,0.00001); b10.2 ~ dnorm(0,0.00001); b10.3 ~ dnorm(0,0.00001); 
	b10.4 ~ dnorm(0,0.00001); b10.5 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100)
  
  b11.0 ~ dnorm(0,0.00001); b11.3 ~ dnorm(0,0.00001); b11.4 ~ dnorm(0,0.00001)  
  sigma.sub ~ dunif(0,100)
}
    ",fill=TRUE)
sink()
# end of JAGS code creation

### CREATE OBJECTS TO HAND OFF TO JAGS  #########################################

N=length(response$plot) # number plots x exp

data = list(N = N,
            pred = as.numeric(response$n_predated),
            n = as.numeric(response$tot_mussels),
            mussel_size = response$mean_size_pred_exp,
            water_temp = response$temp, 
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            sal = response$sal, 
            gd_density = response$gd_density_max_sum,
            burrow = response$burrow, 
            sa_density = response$sa_density,
            sa_height = response$sa_height, 
            sub_prop = response$sub_time)   

parameters <- c("b0.0", "b0.1","b0.2","b0.6", "b0.7","b0.8", "b0.9","b0.10","b0.11",
                
                "b2.0","sigma.temp",
                
                "b6.0","b6.2","b6.2q", "b6.5", "sigma.sal",
                
                "b7.0", "b7.3", "b7.4", "b7.5", "b7.9",
                
                "b8.0","b8.3","b8.4","b8.5",
                
                "b9.0", "b9.3", "b9.4", 
                
                "b10.0", "b10.3", "b10.4", "b10.5","sigma.height", "b10.2",
                
                "b11.0", "b11.3", "b11.4","sigma.sub" )

inits <- function(){list(b0.0 = 0, b0.1 = 0, b0.2 = 0,b0.6 = 0,b0.7 = 0, 
                         b0.8 = 0, b0.9 = 0,b0.10 = 0, b0.11 = 0 ,b0.4=0,
                         
                         b2.0 = 0, sigma.temp = 100,
                         
                         b6.0 = 0, b6.2 = 0,b6.2q = 0,b6.5 = 0, sigma.sal = 100,
                         
                         b7.0 = 0, b7.3 = 0, b7.4 = 0, b7.5 = 0, b7.9 = 0,
                         
                         b8.0 = 0, b8.3 = 0, b8.4 = 0, b8.5 = 0,
                         b9.0 = 0, b9.3 = 0, b9.4 = 0, 
                         
                         b10.0 = 0, b10.2 = 0, b10.3 = 0, b10.4 = 0, b10.5 = 0, 
                         sigma.height = 100,
                         
                         b11.0 = 0, b11.3 = 0, b11.4 = 0,sigma.sub = 100)}


# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain     
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

### RUN MODEL IN JAGS  ##########################################################
global.out <- jags(data=data, 
                            inits=inits, 
                            parameters.to.save=parameters, 
                            model.file="pred.global.model.txt", 
                            n.chains=nc, 
                            n.iter=ni, 
                            n.burnin=nb, 
                            n.thin=nt,  
                            DIC=TRUE,
                            parallel = T)

# Print some basic results
print(global.out,digits=4)  

# Examine traceplots
traceplot(global.out)

# ## Export results for table ##
# hyp_out <- as.data.frame(global.out$summary)  # result summary
# hyp_out <- round(hyp_out,digits = 3)  # round
# hyp_out$median_95 <- paste(hyp_out$`50%`," (",hyp_out$`2.5%`,",",hyp_out$`97.5%`,")",sep ="")  # create column with medians and 95CI
# hyp_out <- hyp_out %>% select(median_95,Rhat,overlap0)
# write.csv(hyp_out,"~/Documents/School/Master's/Thesis documents/ch2 pub/results/model outputs/pred_hyp_out_noblue.csv")

### PREDICTION EQUATIONS ########################################################

# Predation
pred.hat <- global.out$mean$b0.0 +
  global.out$mean$b0.1*response$mean_size_pred_exp+
  global.out$mean$b0.2*response$temp + 
  global.out$mean$b0.6*response$sal +
  global.out$mean$b0.7*response$gd_density_max_sum +
  global.out$mean$b0.8*response$burrow +
  global.out$mean$b0.9*response$sa_density +
  global.out$mean$b0.10*response$sa_height +
  global.out$mean$b0.11*response$sub_time 
pred.hat_exp <- (1/(1+1/(exp(pred.hat))))

# GD density
gd.hat <- exp(global.out$mean$b7.0 +
                global.out$mean$b7.3*response$ele + 
                global.out$mean$b7.4*response$head +
                global.out$mean$b7.5*response$main +
                global.out$mean$b7.9*response$sa_density)

# Submergance
sub.hat <- global.out$mean$b11.0 +
  global.out$mean$b11.3*response$ele +
  global.out$mean$b11.4*response$head

# Salinity
sal.hat <- global.out$mean$b6.0 +
  global.out$mean$b6.2*response$temp + 
  global.out$mean$b6.2q*response$temp^2+
  global.out$mean$b6.5*response$main

# Burrows
bur.hat <- exp(global.out$mean$b8.0 +
                global.out$mean$b8.3*response$ele + 
                global.out$mean$b8.4*response$head +
                global.out$mean$b8.5*response$main)

# SA density
sad.hat <- exp(global.out$mean$b9.0 +
                 global.out$mean$b9.3*response$ele + 
                 global.out$mean$b9.4*response$head)

# SA HEIGHT
sah.hat <- global.out$mean$b10.0 +
  global.out$mean$b10.3*response$ele + 
  global.out$mean$b10.4*response$head +
  global.out$mean$b10.5*response$main +
  global.out$mean$b10.2*response$temp

### MODEL R2 VALUES #############################################################
summary(lm(response$prop_pred ~ pred.hat_exp))
summary(lm(response$gd_density_max_sum~gd.hat))
summary(lm(response$sub_time ~ sub.hat))
summary(lm(response$sal ~ sal.hat))
summary(lm(response$burrow ~ bur.hat))
summary(lm(response$sa_density ~ sad.hat))
summary(lm(response$sa_height ~ sah.hat))


################################################################################
### GLOBAL MODEL AFTER PRUNING NS PATHWAYS  ###
################################################################################

	
sink("pred.global.pruned.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {

## Predation ##	
#0.Pred	 
	  pred[i] ~ dbin(p[i],n[i])                            
		logit(p[i]) <- b0.0 +
	  	b0.2*water_temp[i] + 
	  	b0.6*sal[i] +
	  	b0.10*sa_height[i] +
	  	b0.11*sub_prop[i]  

## Experiment level variables ##	
#2. Water temp
water_temp[i] ~dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0

#6. Salinity
  sal[i] ~ dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <- b6.0 +
    b6.2*water_temp[i] + b6.2q*water_temp[i]*water_temp[i] + 
    b6.5*main[i] 

#10. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b10.0 +   
    b10.5*main[i]

#11. #SUbmergence Time
  sub_prop[i] ~ dnorm(sub.hat[i], tau.sub)
  sub.hat[i] <- b11.0 + 
   b11.3*ele[i] + 
   b11.4*head[i]  
}

### PRECISION VARIABLES  #######################################################
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.height <- 1/(sigma.height*sigma.height)
  tau.temp <- 1/(sigma.temp*sigma.temp)

### PRIORS #####################################################################	
  b0.0 ~ dnorm(0,0.00001);  b0.2 ~ dnorm(0,0.00001) ;b0.6 ~ dnorm(0,0.00001)
	b0.10 ~ dnorm(0,0.00001); b0.11 ~ dnorm(0,0.00001)

	b2.0~ dnorm(0,0.00001);sigma.temp~ dunif(0,100)
	
	b6.0 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); b6.2q ~ dnorm(0,0.00001);
  b6.5 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)

	b10.0 ~ dnorm(0,0.00001); b10.5 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100)
  
  b11.0 ~ dnorm(0,0.00001); b11.3 ~ dnorm(0,0.00001); b11.4 ~ dnorm(0,0.00001)  
  sigma.sub ~ dunif(0,100)
}
    ",fill=TRUE)
sink()
# end of JAGS code creation

### CREATE OBJECTS TO HAND OFF TO JAGS  #########################################

N=length(response$plot) # number of plot x experiment

data = list(N = N,
            pred = as.numeric(response$n_predated),
            n = as.numeric(response$tot_mussels),
            water_temp = response$temp, 
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            sal = response$sal, 
            sa_height = response$sa_height, 
            sub_prop = response$sub_time) 

parameters <- c("b0.0", "b0.2","b0.6","b0.10","b0.11",
                
                "b2.0","sigma.temp",
                
                "b6.0","b6.2","b6.2q", "b6.5", "sigma.sal",
                
                "b10.0","b10.5","sigma.height",
                
                "b11.0", "b11.3", "b11.4","sigma.sub")

inits <- function(){list(b0.0 = 0, b0.2 = 0,b0.6 = 0, b0.10 = 0, b0.11 = 0,
                         
                         b2.0 = 0, sigma.temp = 100,
                         
                         b6.0 = 0, b6.2 = 0,b6.2q = 0,b6.5 = 0, sigma.sal = 100,
                         
                         b10.0 = 0, b10.5 = 0, sigma.height = 100,
                         
                         b11.0 = 0, b11.3 = 0, b11.4 = 0,sigma.sub = 100)}

# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain      
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

### RUN MODEL IN JAGS  ##########################################################
out_ns <- jags(data=data, 
               inits=inits, 
               parameters.to.save=parameters, 
               model.file="pred.global.pruned.model.txt",
               n.chains=nc, 
               n.iter=ni, 
               n.burnin=nb, 
               n.thin=nt,  
               DIC=TRUE,
               parallel = T)


# Print some basic results
print(out_ns,digits=4)  

###  PREDICTION EQUATIONS ######################################################

# Predation
pred.hat <- out_ns$mean$b0.0+
  out_ns$mean$b0.2*response$temp + 
  out_ns$mean$b0.6*response$sal +
  out_ns$mean$b0.10*response$sa_height +
  out_ns$mean$b0.11*response$sub_time
pred.hat_exp <- (1/(1+1/(exp(pred.hat))))

# Salinity
sal.hat <- out_ns$mean$b6.0 +
  out_ns$mean$b6.2*response$temp + 
  out_ns$mean$b6.2q*response$temp^2+
  out_ns$mean$b6.5*response$main

# SA Height
sah.hat <- out_ns$mean$b10.0 +
  out_ns$mean$b10.5*response$main

# Submergance
sub.hat <- out_ns$mean$b11.0 +
  out_ns$mean$b11.3*response$ele +
  out_ns$mean$b11.4*response$head
  
### RESIDUALS ################################################################## 
pred.res <- response$prop_pred-	pred.hat_exp	 
sal.res <- response$sal-sal.hat
sah.res <- response$sa_height-sah.hat
sub.res <- response$sub_time-sub.hat
ele.res <- response$ele
main.res <- response$main
head.res <- response$head
temp.res <- response$temp

# Residual matrix
resid.mat <- data.frame(pred.res,sal.res,sah.res,sub.res,ele.res,main.res,head.res,temp.res)

### MODEL R2 ###################################################################
summary(lm(response$prop_pred ~pred.hat_exp))
summary(lm(response$sal ~sal.hat))
summary(lm(response$sa_height ~sah.hat))
summary(lm(response$sub_time ~sub.hat))

### EXAMINING RESIDUALS FOR MISSING PATHWAYS  ##################################

## Pearson's  ##

# Matrix of pearson correlation p values
pers <- as.data.frame(rstatix::cor_pmat(resid.mat))
rownames(pers)<-pers$rowname  
pers <- pers[,-1]

# Select only comparisons with p value less than 0.05
sig_mask <- pers < 0.05
upper_tri <- upper.tri(sig_mask)
indices <- which(sig_mask & upper_tri, arr.ind = TRUE)

# Create data frame with variable whose residuals are significantly ( p < 0.05)
# correlated
pers_vars <- data.frame(var1 =rownames(pers)[indices[,1]],
                        var2 = colnames(pers)[indices[,2]])

## Spearman's ##

# Matrix of Spearman's correlation p values
spear <- as.data.frame(rstatix::cor_pmat(resid.mat,method = "spearman"))
rownames(spear)<-spear$rowname
spear <- spear[,-1]

# Select only comparisons with p value less than 0.05
sig_mask <- spear < 0.05
upper_tri <- upper.tri(sig_mask)
indices <- which(sig_mask & upper_tri, arr.ind = TRUE)

# Create data frame with variable whose residuals are significantly ( p < 0.05)
# correlated
spear_vars <- data.frame(var1 =rownames(spear)[indices[,1]],
                         var2 = colnames(spear)[indices[,2]])

## Combine both correlation test and find variables that are significantly
# correlated in either
vars <- rbind(pers_vars,spear_vars)
var_combo <- paste(vars$var1,vars$var2,sep = "-")  # combine both variables into one column
var_combo <- var_combo[!duplicated(var_combo)]  # remove duplicates
var_combo

# "sal.res-sub.res"   "ele.res-head.res"  "main.res-head.res" "sub.res-temp.res"  
# "ele.res-main.res"   "sal.res-sub.res"    "sub.res-temp.res"   

################################################################################
### SEM MODEL WITH ADDED PATHWAYS  ###
################################################################################

sink("pred.added.paths.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Predation ##	
#0.Pred	 
	  pred[i] ~ dbin(p[i],n[i])                          
		logit(p[i]) <- b0.0 +
	  	b0.2*water_temp[i] + 
	  	b0.6*sal[i] +
	  	b0.10*sa_height[i] +
	  	b0.11*sub_prop[i]  

## Experiment level variables ##
#2. Water temp
water_temp[i] ~dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0

#6. Salinity
  sal[i] ~ dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <- b6.0 +
    b6.2*water_temp[i] + b6.2q*water_temp[i]*water_temp[i] + 
    b6.5*main[i] +
    b6.11*sub_prop[i] # ADDED

#10. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b10.0 +   
    b10.5*main[i]

#11. #Submergence Time
  sub_prop[i] ~ dnorm(sub.hat[i], tau.sub)
  sub.hat[i] <- b11.0 + 
  b11.2*water_temp[i] + # ADDED
  b11.3*ele[i] + 
  b11.4*head[i]  
}

### PRECISION VARIABLES  #######################################################
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.height <- 1/(sigma.height*sigma.height)
  tau.temp <- 1/(sigma.temp*sigma.temp)

### PRIORS  ####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001) ;b0.6 ~ dnorm(0,0.00001)
	b0.10 ~ dnorm(0,0.00001); b0.11 ~ dnorm(0,0.00001); 

	b2.0~ dnorm(0,0.00001);sigma.temp~ dunif(0,100)
	
  b6.0 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); b6.2q ~ dnorm(0,0.00001)
  b6.5 ~ dnorm(0,0.00001); b6.11~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)
	
	b10.0 ~ dnorm(0,0.00001); b10.5 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100)
  
  b11.0 ~ dnorm(0,0.00001); b11.2 ~ dnorm(0,0.00001); b11.3 ~ dnorm(0,0.00001)  
  b11.4 ~ dnorm(0,0.00001); sigma.sub ~ dunif(0,100)
}
    ",fill=TRUE)
sink()
# end of JAGS code creation

### CREATING OBJECTS TO HAND TO JAGS  ###########################################

N=length(response$plot) # number of plots x experiment

data = list(N = N,
            pred = as.numeric(response$n_predated),
            n = as.numeric(response$tot_mussels),
            water_temp = response$temp, 
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            sal = response$sal, 
            sa_height = response$sa_height, 
            sub_prop = response$sub_time)   

parameters <- c("b0.0", "b0.2","b0.6","b0.10","b0.11",
                
                "b2.0","sigma.temp",
                
                "b6.0","b6.2","b6.2q", "b6.5","b6.11", "sigma.sal",
                
                "b10.0","b10.5","sigma.height",
                
                "b11.0", "b11.2", "b11.3", "b11.4","sigma.sub")

inits <- function(){list(b0.0 = 0, b0.2 = 0, b0.6 = 0,b0.10 = 0, b0.11 = 0, 
                         
                         b2.0 = 0, sigma.temp = 100,
                         
                         b6.0 = 0, b6.2 = 0,b6.2q = 0,b6.5 = 0, b6.11=0, sigma.sal = 100,
                         
                         b10.0 = 0, b10.5 = 0, sigma.height = 100,
                         
                         b11.0 = 0, b11.2=0, b11.3 = 0, b11.4 = 0,sigma.sub = 100)}

# MCMC settings to hand to winbugs
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain      
nb <- 5000   # number of samples to discard for burn in
nt <- 5  

## RUN MODEL IN JAGS ###########################################################
added_paths_out <- jags(data=data, 
                             inits=inits, 
                             parameters.to.save=parameters, 
                             model.file="pred.added.paths.model.txt", 
                             n.chains=nc, 
                             n.iter=ni, 
                             n.burnin=nb, 
                             n.thin=nt,  
                             DIC=TRUE, 
                             parallel = T)

# Print some basic results
print(added_paths_out,digits=4)  

################################################################################
### FINAL PRUNING OF NS PARHWAYS ###
################################################################################

sink("pred.final.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {

## Predation  ##
#0.Pred	 
	  pred[i] ~ dbin(p[i],n[i])                            
		logit(p[i]) <- b0.0 +
	  	b0.2*water_temp[i] + 
	  	b0.6*sal[i] +
	  	b0.10*sa_height[i] +
	  	b0.11*sub_prop[i] 

## Experiment level varaibles ##	
#2. Water temp
water_temp[i] ~dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0

#6. Salinity
  sal[i] ~ dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <- b6.0 +
    b6.2*water_temp[i] + b6.2q*water_temp[i]*water_temp[i] + 
    b6.5*main[i] 

#10. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b10.0 +
    b10.5*main[i]

#11. #Submergence Time
  sub_prop[i] ~ dnorm(sub.hat[i], tau.sub)
  sub.hat[i] <- b11.0 + 
  b11.2*water_temp[i] + b11.2q*water_temp[i]*water_temp[i] +# ADDED
   b11.3*ele[i] + 
   b11.4*head[i]  
}

### PRECISION VARIABLES  #######################################################
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.height <- 1/(sigma.height*sigma.height)
  tau.temp <- 1/(sigma.temp*sigma.temp)

### PRIORS  ####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001) ;b0.6 ~ dnorm(0,0.00001)
	b0.10 ~ dnorm(0,0.00001); b0.11 ~ dnorm(0,0.00001)

	b2.0~ dnorm(0,0.00001);sigma.temp~ dunif(0,100)

  b6.0 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); b6.2q ~ dnorm(0,0.00001)
  b6.5 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)

	b10.0 ~ dnorm(0,0.00001); b10.5 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100)
  
  b11.0 ~ dnorm(0,0.00001);  b11.2 ~ dnorm(0,0.00001);b11.2q ~ dnorm(0,0.00001)
  b11.3 ~ dnorm(0,0.00001); b11.4 ~ dnorm(0,0.00001); sigma.sub ~ dunif(0,100);
}
    ",fill=TRUE)
sink()
# end of JAGS code creation

### CREATING OBJECTS TO HAND TO JAGS  ###########################################

N=length(response$plot) # number of plots x exp

data = list(N = N,
            pred = as.numeric(response$n_predated),
            n = as.numeric(response$tot_mussels),
            water_temp = response$temp, 
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            sal = response$sal, 
            sa_height = response$sa_height, 
            sub_prop = response$sub_time) 

parameters <- c("b0.0", "b0.2","b0.6","b0.10","b0.11",
                
                "b2.0","sigma.temp",
                
                "b6.0","b6.2","b6.2q", "b6.5", "sigma.sal",
                
                "b10.0","b10.5","sigma.height",
                
                "b11.0", "b11.2","b11.2q", "b11.3", "b11.4","sigma.sub")

inits <- function(){list(b0.0 = 0, b0.2 = 0, b0.6 = 0, b0.10 = 0, b0.11 = 0, 
                         
                         b2.0 = 0, sigma.temp = 100,
                         
                         b6.0 = 0, b6.2 = 0,b6.2q = 0,b6.5 = 0, sigma.sal = 100,
                         
                         b10.0 = 0, b10.3 = 0, b10.4 = 0, b10.5 = 0, 
                         sigma.height = 100,
                         
                         b11.0 = 0,b11.2=0, b11.2q=0, b11.3 = 0, b11.4 = 0,
                         sigma.sub = 100)}

# MCMC settings to hand to winbugs
nc <- 3      # number of chains
ni <- 500000  # number of samples for each chain      
nb <- 50000   # number of samples to discard for burn in
nt <- 5  

## RUN MODEL IN JAGS ###########################################################
out_final <- jags(data=data, 
                  inits=inits, 
                  parameters.to.save=parameters, 
                  model.file="pred.final.model.txt", 
                  n.chains=nc, 
                  n.iter=ni, 
                  n.burnin=nb, 
                  n.thin=nt,  
                  DIC=TRUE, 
                  parallel = T)

# Print some basic results
print(out_final,digits=4)  

## Export results for prediction maps ##
# saveRDS(out_final,"model_out/predation_SEM_final_output.rds")

# ## Export results for table ##
# hyp_out <- as.data.frame(out_final$summary)  # result summary
# hyp_out <- round(hyp_out,digits = 3)  # round
# hyp_out$median_95 <- paste(hyp_out$`50%`," (",hyp_out$`2.5%`,",",hyp_out$`97.5%`,")",sep ="")# create column for median and 95CI
# hyp_out <- hyp_out %>% select(median_95,Rhat,overlap0)
# write.csv(hyp_out,"~/Documents/School/Master's/Thesis documents/ch2 pub/results/model outputs/pred_final_out.csv")

### PREDICTION EQUATIONS #######################################################

# Predation
pred.hat <- out_final$mean$b0.0+
  out_final$mean$b0.2*response$temp + 
  out_final$mean$b0.6*response$sal +
  out_final$mean$b0.10*response$sa_height +
  out_final$mean$b0.11*response$sub_time
pred.hat_exp <- (1/(1+1/(exp(pred.hat))))

# Salinity
sal.hat <- out_final$mean$b6.0 +
  out_final$mean$b6.2*response$temp + 
  out_final$mean$b6.2q*response$temp^2+
  out_final$mean$b6.5*response$main

# SA Height
sah.hat <- out_final$mean$b10.0 +
  out_final$mean$b10.5*response$main

# Submergance
sub.hat <- out_final$mean$b11.0 +
  out_final$mean$b11.2*response$temp +
  out_final$mean$b11.2q*response$temp^2 +
  out_final$mean$b11.3*response$ele +
  out_final$mean$b11.4*response$head


### RESIDUALS ###################################################################
pred.res <- response$prop_pred-	pred.hat_exp	 
sal.res <- response$sal-sal.hat
sah.res <- response$sa_height-sah.hat
sub.res <- response$sub_time-sub.hat
ele.res <- response$ele
main.res <- response$main
head.res <- response$head
temp.res <- response$temp

# Residual matrix
final.res <- data.frame(pred.res,sal.res,sah.res,sub.res,ele.res,main.res,head.res,temp.res)

# Export growth residuals
resid_out <- data.frame(
  plot = response$plot,
  season = response$experiment,
  resid = pred.res
)
saveRDS(resid_out,"model_out/pred_resid.rds")

###  FINAL MODEL R2 ############################################################
summary(lm(response$prop_pred ~pred.hat_exp))
summary(lm(response$sal ~sal.hat))
summary(lm(response$sa_height ~sah.hat))
summary(lm(response$sub_time ~sub.hat))

###  EXAMINING RESIDUALS FOR MISSING PATHWAYS  #################################

## Pearson's  ##

# Matrix of pearson correlation p values
pers <- as.data.frame(rstatix::cor_pmat(final.res))
rownames(pers)<-pers$rowname  
pers <- pers[,-1]

# Select only comparisons with p value less than 0.05
sig_mask <- pers < 0.05
upper_tri <- upper.tri(sig_mask)
indices <- which(sig_mask & upper_tri, arr.ind = TRUE)

# Create data frame with variable whose residuals are significantly ( p < 0.05)
# correlated
pers_vars <- data.frame(var1 =rownames(pers)[indices[,1]],
                        var2 = colnames(pers)[indices[,2]])

## Spearman's ##

# Matrix of Spearman's correlation p values
spear <- as.data.frame(rstatix::cor_pmat(final.res,method = "spearman"))
rownames(spear)<-spear$rowname
spear <- spear[,-1]

# Select only comparisons with p value less than 0.05
sig_mask <- spear < 0.05
upper_tri <- upper.tri(sig_mask)
indices <- which(sig_mask & upper_tri, arr.ind = TRUE)

# Create data frame with variable whose residuals are significantly ( p < 0.05)
# correlated
spear_vars <- data.frame(var1 =rownames(spear)[indices[,1]],
                         var2 = colnames(spear)[indices[,2]])

## Combine both correlation test and find variables that are significant
# correlated in either
vars <- rbind(pers_vars,spear_vars)
var_combo <- paste(vars$var1,vars$var2,sep = "-")  # combine both variables into one column
var_combo <- var_combo[!duplicated(var_combo)]  # remove duplicates
var_combo

## all pathways added or not significant, final model ##


################################################################################
### QUERIES IN LIEU OF STANDARDIZED COEFS  ###
################################################################################

# Standardized queries allow for the comparison of effect sizes across non-
# standardized variables. Calculated using the methods of Grace et al., 2012,
# where the change in a response variable between the maximum and minimum values 
# of an explanatory variable, with all other variables held constant. This 
# difference was standardized by the maximum range of the response variable

### FOR PREDATION ##############################################################
  
### TEMP ###

# min
pred.predict.2 <- out_final$mean$b0.0+
  out_final$mean$b0.2*min(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T) +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T) 			      	

# max
pred.predict.3 <- out_final$mean$b0.0+
  out_final$mean$b0.2*max(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T) +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T) 			        

# using straight probabilities
pred.predict.2pr.raw <- (1/(1+1/(exp(pred.predict.2))))  		
pred.predict.3pr.raw <- (1/(1+1/(exp(pred.predict.3))))  		

# QUERY
(pred.predict.3pr.raw-pred.predict.2pr.raw)/(max(response$prop_pred)-min(response$prop_pred))

### SALINITY ###

# min
pred.predict.4 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*min(response$sal,na.rm = T) +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T) 			      

#max
pred.predict.5 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*max(response$sal,na.rm = T) +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T) 			      	 

# using straight probabilities
pred.predict.4pr.raw <- (1/(1+1/(exp(pred.predict.4))))
pred.predict.5pr.raw <- (1/(1+1/(exp(pred.predict.5))))

# QUERY
(pred.predict.5pr.raw-pred.predict.4pr.raw)/(max(response$prop_pred)-min(response$prop_pred))

### SA HEIGHT ###

# min
pred.predict.6 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T) +
  out_final$mean$b0.10*min(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T) 			      

# max
pred.predict.7 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T) +
  out_final$mean$b0.10*max(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T) 			      	 

# using straight probabilities
pred.predict.6pr.raw <- (1/(1+1/(exp(pred.predict.6))))
pred.predict.7pr.raw <- (1/(1+1/(exp(pred.predict.7))))

# QUERY
(pred.predict.7pr.raw-pred.predict.6pr.raw)/(max(response$prop_pred)-min(response$prop_pred))

### SUBMERGANCE ###

# min
pred.predict.8 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T) +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*min(response$sub_time,na.rm = T) 

# max
pred.predict.9 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T) +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*max(response$sub_time,na.rm = T)			      	 

# using straight probabilities
pred.predict.8pr.raw <- (1/(1+1/(exp(pred.predict.8))))
pred.predict.9pr.raw <- (1/(1+1/(exp(pred.predict.9))))

# QUERY
(pred.predict.9pr.raw-pred.predict.8pr.raw)/(max(response$prop_pred)-min(response$prop_pred))


### FOR SALINITY ###############################################################

### TEMP ###

# min
sal.predict.2 <- out_final$mean$b6.0 +
  out_final$mean$b6.2*min(response$temp,na.rm = T) + 
  out_final$mean$b6.2q*min(response$temp,na.rm = T)^2+
  out_final$mean$b6.5*mean(response$main)

# max
sal.predict.3 <- out_final$mean$b6.0 +
  out_final$mean$b6.2*max(response$temp,na.rm = T) + 
  out_final$mean$b6.2q*max(response$temp,na.rm = T)^2+
  out_final$mean$b6.5*mean(response$main)

# QUERY
(sal.predict.3-sal.predict.2)/(max(response$sal,na.rm = T)-min(response$sal,na.rm = T))

### Main ###

# min
sal.predict.6 <- out_final$mean$b6.0 +
  out_final$mean$b6.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b6.2q*mean(response$temp,na.rm = T)^2+
  out_final$mean$b6.5*min(response$main)

# max
sal.predict.7 <- out_final$mean$b6.0 +
  out_final$mean$b6.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b6.2q*mean(response$temp,na.rm = T)^2+
  out_final$mean$b6.5*max(response$main)

# QUERY
(sal.predict.7-sal.predict.6)/(max(response$sal,na.rm = T)-min(response$sal,na.rm = T))

### FOR SA HEIGHT  #############################################################

### MAIN  ####

# min
sah.predict.2 <- out_final$mean$b10.0 +
  out_final$mean$b10.5*min(response$main)

# max
sah.predict.3 <- out_final$mean$b10.0 +
  out_final$mean$b10.5*max(response$main)

# QUERY
(sah.predict.3-sah.predict.2)/(max(response$sa_height,na.rm = T)-min(response$sa_height,na.rm = T))

### FOR SUB ####################################################################

### TEMP ###

# min
sub.predict.2 <- out_final$mean$b11.0 +
  out_final$mean$b11.2*min(response$temp,na.rm = T) +
  out_final$mean$b11.2q*min(response$temp,na.rm = T)^2 +
  out_final$mean$b11.3*mean(response$ele) +
  out_final$mean$b11.4*mean(response$head)

# max
sub.predict.3 <- out_final$mean$b11.0 +
  out_final$mean$b11.2*max(response$temp,na.rm = T) +
  out_final$mean$b11.2q*max(response$temp,na.rm = T)^2 +
  out_final$mean$b11.3*mean(response$ele) +
  out_final$mean$b11.4*mean(response$head) 

# QUERY
(sub.predict.3-sub.predict.2)/(max(response$sub_time,na.rm = T)-min(response$sub_time,na.rm = T))

### ELEVATION ###

# min
sub.predict.6 <- out_final$mean$b11.0 +
  out_final$mean$b11.2*mean(response$temp,na.rm = T) +
  out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
  out_final$mean$b11.3*min(response$ele) +
  out_final$mean$b11.4*mean(response$head) 

# max
sub.predict.7 <- out_final$mean$b11.0 +
  out_final$mean$b11.2*mean(response$temp,na.rm = T) +
  out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
  out_final$mean$b11.3*max(response$ele) +
  out_final$mean$b11.4*mean(response$head) 

# QUERY
(sub.predict.7-sub.predict.6)/(max(response$sub_time,na.rm = T)-min(response$sub_time,na.rm = T))

### HEAD ###

# min
sub.predict.8 <- out_final$mean$b11.0 +
  out_final$mean$b11.2*mean(response$temp,na.rm = T) +
  out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
  out_final$mean$b11.3*mean(response$ele) +
  out_final$mean$b11.4*min(response$head) 

# max
sub.predict.9 <- out_final$mean$b11.0 +
  out_final$mean$b11.2*mean(response$temp,na.rm = T) +
  out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
  out_final$mean$b11.3*mean(response$ele) +
  out_final$mean$b11.4*max(response$head) 

# QUERY
(sub.predict.9-sub.predict.8)/(max(response$sub_time,na.rm = T)-min(response$sub_time,na.rm = T))

################################################################################
## TOTAL EFFECTS QUERIES ##
################################################################################

## TEMP ###

# Min
pred.predict.tot2  <- out_final$mean$b0.0+
  out_final$mean$b0.2*min(response$temp,na.rm = T) + 
  out_final$mean$b0.6*(out_final$mean$b6.0 +
     out_final$mean$b6.2*min(response$temp,na.rm = T) + 
     out_final$mean$b6.2q*min(response$temp,na.rm = T)^2+
     out_final$mean$b6.5*mean(response$main))  +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 +
    out_final$mean$b11.2*min(response$temp,na.rm = T) +
    out_final$mean$b11.2q*min(response$temp,na.rm = T)^2 +
    out_final$mean$b11.3*mean(response$ele) +
    out_final$mean$b11.4*mean(response$head))	  
pred.predict.tot2pr.raw <- (1/(1+1/(exp(pred.predict.tot2))))

# Max
pred.predict.tot3 <- pred.hat <- out_final$mean$b0.0+
  out_final$mean$b0.2*max(response$temp,na.rm = T) + 
  out_final$mean$b0.6*(out_final$mean$b6.0 +
     out_final$mean$b6.2*max(response$temp,na.rm = T) + 
     out_final$mean$b6.2q*max(response$temp,na.rm = T)^2+
     out_final$mean$b6.5*mean(response$main))  +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 +
      out_final$mean$b11.2*max(response$temp,na.rm = T) +
      out_final$mean$b11.2q*max(response$temp,na.rm = T)^2 +
      out_final$mean$b11.3*mean(response$ele) +
      out_final$mean$b11.4*mean(response$head))
pred.predict.tot3pr.raw <- (1/(1+1/(exp(pred.predict.tot3))));

# QUERY
(pred.predict.tot3pr.raw-pred.predict.tot2pr.raw)/(max(response$prop_pred)-min(response$prop_pred))


### ELEVATION ###

# min
pred.predict.tot4 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T)  +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 +
      out_final$mean$b11.2*mean(response$temp,na.rm = T) +
      out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
      out_final$mean$b11.3*min(response$ele) +
      out_final$mean$b11.4*mean(response$head))
pred.predict.tot4pr.raw <- (1/(1+1/(exp(pred.predict.tot4))))

# max
pred.predict.tot5 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T)  +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 +
      out_final$mean$b11.2*mean(response$temp,na.rm = T) +
      out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
      out_final$mean$b11.3*max(response$ele) +
      out_final$mean$b11.4*mean(response$head)) 
pred.predict.tot5pr.raw <- (1/(1+1/(exp(pred.predict.tot5))))

# QUERY
(pred.predict.tot5pr.raw-pred.predict.tot4pr.raw)/(max(response$prop_pred)-min(response$prop_pred))


### HEAD ###

# min
pred.predict.tot4 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T)  +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 +
      out_final$mean$b11.2*mean(response$temp,na.rm = T) +
      out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
      out_final$mean$b11.3*mean(response$ele) +
      out_final$mean$b11.4*min(response$head))	  
pred.predict.tot4pr.raw <- (1/(1+1/(exp(pred.predict.tot4))))

#max
pred.predict.tot5 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*mean(response$sal,na.rm = T)  +
  out_final$mean$b0.10*mean(response$sa_height,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 +
      out_final$mean$b11.2*mean(response$temp,na.rm = T) +
      out_final$mean$b11.2q*mean(response$temp,na.rm = T)^2 +
      out_final$mean$b11.3*mean(response$ele) +
      out_final$mean$b11.4*max(response$head))  
pred.predict.tot5pr.raw <- (1/(1+1/(exp(pred.predict.tot5))))

# QUERY
(pred.predict.tot5pr.raw-pred.predict.tot4pr.raw)/(max(response$prop_pred)-min(response$prop_pred))


### MAIN ###

# min
pred.predict.tot4 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*(out_final$mean$b6.0 +
       out_final$mean$b6.2*mean(response$temp,na.rm = T) + 
       out_final$mean$b6.2q*mean(response$temp,na.rm = T)^2+
       out_final$mean$b6.5*min(response$main))  +
  out_final$mean$b0.10*(out_final$mean$b10.0 +
      out_final$mean$b10.5*min(response$main)) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T)
pred.predict.tot4pr.raw <- (1/(1+1/(exp(pred.predict.tot4))))

#max
pred.predict.tot5 <- out_final$mean$b0.0+
  out_final$mean$b0.2*mean(response$temp,na.rm = T) + 
  out_final$mean$b0.6*(out_final$mean$b6.0 +
     out_final$mean$b6.2*mean(response$temp,na.rm = T) + 
     out_final$mean$b6.2q*mean(response$temp,na.rm = T)^2+
     out_final$mean$b6.5*max(response$main))  +
  out_final$mean$b0.10*(out_final$mean$b10.0 +
    out_final$mean$b10.5*max(response$main)) +
  out_final$mean$b0.11*mean(response$sub_time,na.rm = T)  
pred.predict.tot5pr.raw <- (1/(1+1/(exp(pred.predict.tot5))))

# QUERY
(pred.predict.tot5pr.raw-pred.predict.tot4pr.raw)/(max(response$prop_pred)-min(response$prop_pred))

