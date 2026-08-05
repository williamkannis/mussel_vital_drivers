################################################################################
#
#  RIBBED MUSSEL RECRUITMENT STRUCTURAL EQUATION MODELS   
#  
#  Created by Anonymized                               
#  Modified from Grace et al., 2012 - Ecosphere           
#
################################################################################

# This code performs the structural equation modelling (SEM) of the hierarchical
# drivers of ribbed mussel recruitment across the salt marsh landscape. Here,
# we create a global Bayesian piecewise SEM based on prior knowledge from 
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

## recruitment and explanatory data  ##
response <- read_excel("mussel_population_vital_data.xlsx",sheet = 3)

#Order by experiment
response <- response[order(response$season),]

## Check for correlation among explanatory variables  ##
response_cov <- response[,3:22]
response_TF_cov<-abs(cor(response_cov, use = "complete.obs"))>.60

## Check for missing data ##

# what columns have missing data?
colnames(response)[colSums(is.na(response)) > 0]  

# fix the one NA value in crab burrows by assigning its yearly average value
response$burrow[is.na(response$burrow)] <- 
  as.integer(response$burrow_m[is.na(response$burrow)])  

################################################################################
### VARIABLE CODES FOR MODELS ###
################################################################################

#1. gd_density
#2. temp
#3. sub_prop
#4. flow_rate
#5. mean_sal
#10. sa_density
#11. sa_height
#12. ele
#13. head
#14. main
#15. burrow


################################################################################
### Global SEM MODEL  ###
################################################################################

sink("recruit.global.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Recruitment ## 
	  recruit[i] ~ dpois(mu.recruit[i])                         
		log(mu.recruit[i]) <- b0.0 + 
	  	b0.1*gd_density[i] + 
	  	b0.2*temp[i] + 
	  	b0.3*sub_prop[i] + 
	  	b0.4*flow_rate[i] + 
	  	b0.5*mean_sal[i] + 
	  	b0.10*sa_density[i] +
	  	b0.11*sa_height[i] +
	  	b0.14*main[i] +
	  	b0.15*burrow[i] 

## Experiment level variables ##
#1. Mussel density
  gd_density[i] ~ dpois(gd.hat[i])
  log(gd.hat[i]) <- b1.0 +
    b1.10*sa_density[i] +
    b1.12*ele[i]+
    b1.13*head[i]+
    b1.14*main[i]	  	

#2. Temp
  temp[i] ~ dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0  + 
    b2.3*sub_prop[i] + 
      b2.10*sa_density[i] +
      b2.11*sa_height[i] 

#5. Salinity
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]

#10. Sa_denisty
  sa_density[i] ~dpois(sad.hat[i])
  log(sad.hat[i]) <- b10.0 +
    b10.12*ele[i]+
    b10.13*head[i]

#11. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b11.0 +   
    b11.12*ele[i]+
    b11.13*head[i]+
    b11.14*main[i]

#15. Crab burrow
burrow[i] ~ dpois(burrow.hat[i])
  log(burrow.hat[i]) <- b15.0 +   
    b15.12*ele[i]+
    b15.13*head[i]+
    b15.14*main[i]
	}
	
## Plot level variables ##	
for(j in 1:30) {
#3. SUbmergence Time
  sub_prop[j] ~ dnorm(sub.hat[j], tau.sub)
  sub.hat[j] <- b3.0 + 
   b3.12*ele[j] + 
   b3.13*head[j]

#4. Flow
  flow_rate[j] ~ dnorm(flow.hat[j], tau.flow)
    flow.hat[j] <- b4.0 + 
      b4.10*sa_density_m[j] +
      b4.11*sa_height_m[j] +
      b4.12*ele[j]+
      b4.13*head[j]
	}
	
### PRECISION VARIABLES  #######################################################
	tau.growth <- 1/(sigma.growth*sigma.growth)
	tau.gd <- 1/(sigma.gd*sigma.gd)
	tau.temp <- 1/(sigma.temp*sigma.temp)
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.flow <- 1/(sigma.flow*sigma.flow)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.sad <- 1/(sigma.sad*sigma.sad)
	tau.height <- 1/(sigma.height*sigma.height)
	
### PRIORS #####################################################################
  b0.0 ~ dnorm(0,0.00001); b0.1 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001) 
	b0.3 ~ dnorm(0,0.00001); b0.4 ~ dnorm(0,0.00001); b0.5 ~ dnorm(0,0.00001)
	b0.10 ~ dnorm(0,0.00001);b0.11 ~ dnorm(0,0.00001); 	b0.14 ~ dnorm(0,0.00001)
	b0.15 ~ dnorm(0,0.00001);sigma.growth ~ dunif(0,100);
	
	b1.0 ~ dnorm(0,0.00001); b1.10 ~ dnorm(0,0.00001); b1.12 ~ dnorm(0,0.00001) 
	b1.13 ~ dnorm(0,0.00001); b1.14 ~ dnorm(0,0.00001); sigma.gd ~ dunif(0,100);
	
	b2.0 ~ dnorm(0,0.00001); b2.3 ~ dnorm(0,0.00001); b2.10 ~ dnorm(0,0.00001) 
	b2.11 ~ dnorm(0,0.00001);b2.12 ~ dnorm(0,0.00001); sigma.temp ~ dunif(0,100)
	
  b3.0 ~ dnorm(0,0.00001); b3.12 ~ dnorm(0,0.00001); b3.13 ~ dnorm(0,0.00001)
  sigma.sub ~ dunif(0,100)
  
  b4.0 ~ dnorm(0,0.00001); b4.10 ~ dnorm(0,0.00001); b4.11 ~ dnorm(0,0.00001) 
	b4.12 ~ dnorm(0,0.00001); b4.13 ~ dnorm(0,0.00001); sigma.flow ~ dunif(0,100)

  b5.0 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);b5.2q ~ dnorm(0,0.00001)
  b5.14 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)
  
  b10.0 ~ dnorm(0,0.00001);  b10.12 ~ dnorm(0,0.00001)
  b10.13 ~ dnorm(0,0.00001); b10.14 ~ dnorm(0,0.00001); sigma.sad ~ dunif(0,100);
  
  b11.0 ~ dnorm(0,0.00001);  b11.12 ~ dnorm(0,0.00001)
  b11.13 ~ dnorm(0,0.00001); b11.14 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100);

  b15.0 ~ dnorm(0,0.00001);  b15.12 ~ dnorm(0,0.00001)
  b15.13 ~ dnorm(0,0.00001); b15.14 ~ dnorm(0,0.00001)
}
",fill=TRUE)
sink()
# end of JAGS code creation

### CREATE OBJECTS TO HAND OFF TO JAGS  #########################################

N=90 # Number of plots*seasons

data = list(N = N,
            recruit = response$recruit,
            gd_density = response$gd_density_max_sum,
            temp = response$mean_water_temp, 
            sub_prop = response$sub_prop, 
            flow_rate = response$flow_rate, 
            mean_sal = response$mean_sal,
            sa_density = response$sa_density,
            sa_density_m = response$sa_denisty_m,
            sa_height = response$sa_height,
            sa_height_m = response$sa_height_m,
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            burrow=response$burrow)

parameters <- c("b0.0","b0.1","b0.2","b0.3","b0.4","b0.5","b0.10","b0.11",
                "b0.15","b0.14",
                
                "b1.0","b1.10","b1.12","b1.13","b1.14","sigma.gd",
                
                "b2.0","b2.3","b2.10","b2.11","b2.12","sigma.temp",
                
                "b3.0", "b3.12","b3.13","sigma.sub",
                
                "b4.0","b4.10","b4.11","b4.12","b4.13","sigma.flow",
                
                "b5.0","b5.2","b5.2q","b5.14","sigma.sal",

                "b10.0","b10.12","b10.13","b10.14","sigma.sad",
                
                "b11.0","b11.10","b11.12","b11.13","b11.14","sigma.height",
                
                "b15.0","b15.12","b15.13","b15.14")

inits <- function(){list(b0.0=0,b0.1=0,b0.2=0,b0.3=0,b0.4=0,b0.5=0,                                
                         b0.10=0,b0.11=0,b0.14=0,sigma.growth=100, 
                         
                         b1.0=0,b1.10=0,b1.12=0,b1.13=0,b1.14=0,sigma.gd=100,
                         
                         b2.0=0, b2.3=0,b2.10=0,b2.11=0,b2.12=0, sigma.temp=100,
                         
                         b3.0=0, b3.12=0, b3.13=0, sigma.sub=100,
                         
                         b4.0=0,b4.10=0,b4.11=0,b4.12=0,b4.13=0,sigma.flow=100,
                         
                         b5.0=0,b5.2=0,b5.2q=0,b5.14=0,sigma.sal=100,
                         
                         b10.0=0,b10.12=0,b10.13=0,b10.14=0,sigma.sad=100,
                         
                         b11.0=0,b11.12=0,b11.13=0,b11.14=0,sigma.height=100,
                         
                         b15.0=0,b15.12=0,b15.13=0,b15.14=0)}

# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain      
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

### RUN MODEL IN JAGS  ##########################################################
global.out <- jags(data=data, 
                            inits=inits, 
                            parameters.to.save=parameters, 
                            model.file="recruit.global.model.txt", 
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

## Export results for table ##
# hyp_out <- as.data.frame(global.out$summary)  # result summary
# hyp_out <- round(hyp_out,digits = 3)  # round results
# hyp_out$median_95 <- paste(hyp_out$`50%`," (",hyp_out$`2.5%`,",",hyp_out$`97.5%`,")",sep ="")  # create column with medians and 95CI
# hyp_out <- hyp_out %>% select(median_95,Rhat,overlap0)
# write.csv(hyp_out,"~/Documents/School/Master's/Thesis documents/ch2 pub/results/model outputs/recruit_hyp_out.csv")

### PREDICTION EQUATIONS ########################################################

# Recruitment
recruit.hat <- exp(global.out$mean$b0.0 + 
  global.out$mean$b0.1*response$gd_density_max_sum +
  global.out$mean$b0.2*response$mean_water_temp + 
  global.out$mean$b0.3*response$sub_prop +
  global.out$mean$b0.4*response$flow_rate +
  global.out$mean$b0.5*response$mean_sal + 
  global.out$mean$b0.10*response$sa_density +
  global.out$mean$b0.11*response$sa_height + 
  global.out$mean$b0.14*response$main +
  global.out$mean$b0.15*response$burrow)

# GD density
gd.hat <- exp(global.out$mean$b1.0 +
                global.out$mean$b1.10*response$sa_density +
                global.out$mean$b1.12*response$ele +
                global.out$mean$b1.13*response$head +
                global.out$mean$b1.14*response$main)

# Temp
temp.hat <- global.out$mean$b2.0 + 
  global.out$mean$b2.3*response$sub_prop +
  global.out$mean$b2.10*response$sa_density +
  global.out$mean$b2.11*response$sa_height 

# Sumbergance
sub.hat <- global.out$mean$b3.0 +
  global.out$mean$b3.12*response$ele +
  global.out$mean$b3.13*response$head 

# Flow
flow.hat <- global.out$mean$b4.0 + 
  global.out$mean$b4.10*response$sa_density +
  global.out$mean$b4.11*response$sa_height +
  global.out$mean$b4.12*response$ele +
  global.out$mean$b4.13*response$head

# Salinity
sal.hat <- global.out$mean$b5.0 + 
  global.out$mean$b5.2*response$mean_water_temp+ 
  global.out$mean$b5.2q*response$mean_water_temp^2+
  global.out$mean$b5.14*response$main

# SA density
sad.hat <- exp(global.out$mean$b10.0 + 
  global.out$mean$b10.12*response$ele +
  global.out$mean$b10.13*response$head)

# SA Heigth
sah.hat <- global.out$mean$b11.0 + 
  global.out$mean$b11.12*response$ele + 
  global.out$mean$b10.13*response$head +
  global.out$mean$b11.14*response$main

# Burrows
crab.hat <- exp(global.out$mean$b15.0 + 
  global.out$mean$b15.12*response$ele + 
  global.out$mean$b15.13*response$head + 
  global.out$mean$b15.14*response$main)

### MODEL R2 VALUES #############################################################
summary(lm(response$recruit~recruit.hat))
summary(lm(response$gd_density_max_sum~gd.hat))
summary(lm(response$mean_water_temp~temp.hat))
summary(lm(response$sub_prop~sub.hat))
summary(lm(response$flow_rate~flow.hat))
summary(lm(response$mean_sal~sal.hat))
summary(lm(response$sa_density~sad.hat))
summary(lm(response$sa_height~sah.hat))
summary(lm(response$burrow~crab.hat))

################################################################################
### GLOBAL MODEL AFTER PRUNING NS PATHWAYS  ###
################################################################################

sink("recruit.global.pruned.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Recruitment ##	 
	  recruit[i] ~ dpois(mu.recruit[i])                           # 
		log(mu.recruit[i]) <- b0.0 + 
	  	b0.2*temp[i] + 
	  	b0.3*sub_prop[i] +
	  	b0.5*mean_sal[i] + 
	  	b0.10*sa_density[i] +
	  	b0.11*sa_height[i] +
	  	b0.14*main[i] +
	  	b0.15*burrow[i]
	  	
## Experiment level variables ##	 
#2. Temp
  temp[i] ~ dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0  + 
      b2.11*sa_height[i] 

#5. sal
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]

#10. Sa_denisty
  sa_density[i] ~dpois(sad.hat[i])
  log(sad.hat[i]) <- b10.0 +
    b10.13*head[i]

#11. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b11.0 +   
    b11.12*ele[i]+
    b11.14*main[i]

#15. Crab burrow
burrow[i] ~ dpois(burrow.hat[i])
  log(burrow.hat[i]) <- b15.0 +   
    b15.12*ele[i]+
    b15.13*head[i]+
    b15.14*main[i]
}

## Plot level variables  ##
for(j in 1:30) {
#3. Submergence Time
  sub_prop[j] ~ dnorm(sub.hat[j], tau.sub)
  sub.hat[j] <- b3.0 + 
   b3.12*ele[j] + 
   b3.13*head[j]
	}
	
### PRECISION VARIABLES  #######################################################
	tau.growth <- 1/(sigma.growth*sigma.growth)
	tau.temp <- 1/(sigma.temp*sigma.temp)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.sad <- 1/(sigma.sad*sigma.sad)
	tau.height <- 1/(sigma.height*sigma.height)
	tau.sub <- 1/(sigma.sub*sigma.sub)
	
	
### PRIORS #####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001);b0.3 ~ dnorm(0,0.00001)
	b0.5 ~ dnorm(0,0.00001); b0.10 ~ dnorm(0,0.00001); b0.11 ~ dnorm(0,0.00001)	
	b0.14 ~ dnorm(0,0.00001);b0.15 ~ dnorm(0,0.00001);sigma.growth ~ dunif(0,100)
	
	b2.0 ~ dnorm(0,0.00001); b2.11 ~ dnorm(0,0.00001);sigma.temp ~ dunif(0,100)
	 
	b3.0 ~ dnorm(0,0.00001); b3.12 ~ dnorm(0,0.00001); b3.13 ~ dnorm(0,0.00001); sigma.sub ~ dunif(0,100)

  b5.0 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);b5.2q ~ dnorm(0,0.00001)
  b5.14 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)
  
  b10.0 ~ dnorm(0,0.00001); b10.13 ~ dnorm(0,0.00001);sigma.sad ~ dunif(0,100);
  
  b11.0 ~ dnorm(0,0.00001);  b11.12 ~ dnorm(0,0.00001);
  b11.14 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100);

  b15.0 ~ dnorm(0,0.00001);  b15.12 ~ dnorm(0,0.00001)
  b15.13 ~ dnorm(0,0.00001); b15.14 ~ dnorm(0,0.00001)
}
",fill=TRUE)
sink()
# end of JAGS code creation

### CREATE OBJECTS TO HAND OFF TO JAGS  #########################################

N=90 # plots x seasons

data = list(N = N,
            recruit = response$recruit,
            gd_density = response$gd_density_max_sum,
            temp = response$mean_water_temp, 
            sub_prop = response$sub_prop, 
            mean_sal = response$mean_sal,
            sa_density = response$sa_density,
            sa_height = response$sa_height,
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            burrow=response$burrow)

parameters <- c("b0.0","b0.2","b0.3","b0.5","b0.10","b0.11","b0.14", "b0.15",
                
                "b2.0", "b2.11","sigma.temp",
                
                "b3.0", "b3.12","b3.13","sigma.sub",
                
                "b5.0","b5.14","b5.2","b5.2q","sigma.sal",
                
                "b10.0","b10.13","sigma.sad",
                
                "b11.0","b11.12","b11.14","sigma.height",
                
                "b15.0","b15.12","b15.13","b15.14")

inits <- function(){list(b0.0=0,b0.2=0,b0.3=0,b0.5=0,b0.10=0,b0.11=0,b0.14=0,b0.15=0,
                         sigma.growth=100, 
   
                         b2.0=0, b2.11=0,sigma.temp=100,
                         
                         b5.0=0,b5.2=0,b5.2q=0,b5.14=0,sigma.sal=100,
                         
                         b10.0=0,b10.13=0,sigma.sad=100,
                         
                         b11.0=0,b11.12=0,b11.14=0,sigma.height=100,
                         
                         b15.0=0,b15.12=0,b15.13=0,b15.14=0)}

# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain     
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

### RUN MODEL IN JAGS  ##########################################################
out_ns <- jags(data=data, 
               inits=inits, 
               parameters.to.save=parameters, 
               model.file="recruit.global.pruned.model.txt", 
               n.chains=nc, 
               n.iter=ni, 
               n.burnin=nb, 
               n.thin=nt,  
               DIC=TRUE, 
               parallel = T)

# Print some basic results
print(out_ns,digits=4)  

###  PREDICTION EQUATIONS ######################################################

# Recruitment
recruit.hat <- out_ns$mean$b0.0 + 
  out_ns$mean$b0.2*response$mean_water_temp + 
  out_ns$mean$b0.3*response$sub_prop+
  out_ns$mean$b0.5*response$mean_sal + 
  out_ns$mean$b0.10*response$sa_density +
  out_ns$mean$b0.11*response$sa_height + 
  out_ns$mean$b0.14*response$main +
  out_ns$mean$b0.15*response$burrow
recruit.hat_exp <- exp(recruit.hat)  

# Temp
temp.hat <- out_ns$mean$b2.0 + 
  out_ns$mean$b2.11*response$sa_height 

# Submergance
sub.hat <- out_ns$mean$b3.0 +
  out_ns$mean$b3.12*response$ele +
  out_ns$mean$b3.13*response$head 

# Salinity
sal.hat <- out_ns$mean$b5.0 + 
  out_ns$mean$b5.2*response$mean_water_temp + 
  out_ns$mean$b5.2q*response$mean_water_temp^2 +
  out_ns$mean$b5.14*response$main

# SA density
sad.hat <- exp(out_ns$mean$b10.0 + 
  out_ns$mean$b10.13*response$head)

# SA Height
sah.hat <- out_ns$mean$b11.0 + 
  out_ns$mean$b11.12*response$ele + 
  out_ns$mean$b11.14*response$main

# Burrows
crab.hat <- exp(out_ns$mean$b15.0 + 
  out_ns$mean$b15.12*response$ele + 
  out_ns$mean$b15.13*response$head + 
  out_ns$mean$b15.14*response$main)



### RESIDUALS ################################################################## 
recruit.res <- response$recruit-recruit.hat_exp
temp.res <- response$mean_water_temp-temp.hat 	
sub.res <- response$sub_prop-sub.hat
sal.res <- response$mean_sal-sal.hat 	 
sad.res <- response$sa_density-sad.hat		
sah.res <- response$sa_height-sah.hat 	  
crab.res <- response$burrow-crab.hat	  
ele.res <- response$ele 	 
main.res <- response$main  
head.res <- response$head	  

# Residual matrix
resid.mat <- data.frame(recruit.res,temp.res,sad.res,sah.res,crab.res,sal.res,ele.res,main.res,head.res,sub.res)


### MODEL R2 ################################################################### 
summary(lm(response$recruit ~recruit.hat_exp))
summary(lm(response$mean_water_temp~temp.hat))
summary(lm(response$sub_prop~sub.hat))
summary(lm(response$mean_sal~sal.hat))
summary(lm(response$sa_density~sad.hat))
summary(lm(response$sa_height~sah.hat))
summary(lm(response$burrow~crab.hat))

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

# sig combos
# "sad.res-crab.res"     "sah.res-crab.res"     "ele.res-head.res"     
# "main.res-head.res"    "sah.res-sub.res"      "recruit.res-crab.res"
# "ele.res-main.res"  

# meaningful
# [1] "sad.res-crab.res"     "sah.res-crab.res"            
#          

## SIg path ways to add
# main->sad #ADDED
# sad->crab # ADDED
# sah->crab  #ADDED



################################################################################
### SEM MODEL WITH ADDED PATHWAYS  ###
################################################################################

sink("recruit.added.paths.model.txt")		
cat("
model {
	# Likelihoods
	for(i in 1:N) {
### LIKELIHOODS  ###############################################################	

## Recruitment ##	 
	  recruit[i] ~ dpois(mu.recruit[i])                            
		log(mu.recruit[i]) <- b0.0 + 
	  	b0.2*temp[i] + 
	  	b0.3*sub_prop[i] +
	  	b0.5*mean_sal[i] + 
	  	b0.10*sa_density[i] +
	  	b0.11*sa_height[i] +
	  	b0.14*main[i] +
	  	b0.15*burrow[i] 

## Experiment level variables ##	 
#5. sal
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]

#10. Sa_denisty
  sa_density[i] ~dpois(sad.hat[i])
  log(sad.hat[i]) <- b10.0 +
    b10.13*head[i] +
    b10.14*main[i]  # ADDED

#11. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b11.0 +   
    b11.12*ele[i]+
    b11.14*main[i]+
    b11.2*temp[i] #ADDED

#15. Crab burrow
burrow[i] ~ dpois(burrow.hat[i])
  log(burrow.hat[i]) <- b15.0 +   
    b15.12*ele[i]+
    b15.13*head[i]+
    b15.14*main[i] +
    b15.10*sa_density[i] +  #ADDED
	  b15.11*sa_height[i]  # ADDED
}

## Plot level variables ##
for(j in 1:30) {
#3. SUbmergence Time
  sub_prop[j] ~ dnorm(sub.hat[j], tau.sub)
  sub.hat[j] <- b3.0 + 
   b3.12*ele[j] + 
   b3.13*head[j]
	}
	
### PRECISION VARIABLES  #######################################################
	tau.growth <- 1/(sigma.growth*sigma.growth)
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.sad <- 1/(sigma.sad*sigma.sad)
	tau.height <- 1/(sigma.height*sigma.height)

### PRIORS  ####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001);b0.3 ~ dnorm(0,0.00001)
	b0.5 ~ dnorm(0,0.00001); b0.10 ~ dnorm(0,0.00001); b0.11 ~ dnorm(0,0.00001) 
	b0.14 ~ dnorm(0,0.00001); b0.15 ~ dnorm(0,0.00001);sigma.growth ~ dunif(0,100)
	
	b3.0 ~ dnorm(0,0.00001); b3.12 ~ dnorm(0,0.00001); b3.13 ~ dnorm(0,0.00001)
	sigma.sub ~ dunif(0,100)

  b5.0 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);b5.2q ~ dnorm(0,0.00001)
  b5.14 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)
  
  b10.0 ~ dnorm(0,0.00001); b10.13 ~ dnorm(0,0.00001);b10.14 ~ dnorm(0,0.00001) 
  sigma.sad ~ dunif(0,100)
  
  b11.0 ~ dnorm(0,0.00001); b11.2 ~ dnorm(0,0.00001); b11.12 ~ dnorm(0,0.00001) 
  b11.14 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100)

  b15.0 ~ dnorm(0,0.00001); b15.10 ~ dnorm(0,0.00001); b15.11 ~ dnorm(0,0.00001)  
  b15.12 ~ dnorm(0,0.00001);b15.13 ~ dnorm(0,0.00001); b15.14 ~ dnorm(0,0.00001)
}
",fill=TRUE)
sink()
# end of JAGS code creation

### CREATING OBJECTS TO HAND TO JAGS  ###########################################
N=90 # number of plot x season

data = list(N = N,
            recruit = response$recruit,
            gd_density = response$gd_density_max_sum,
            temp = response$mean_water_temp, 
            sub_prop = response$sub_prop, 
            mean_sal = response$mean_sal,
            sa_density = response$sa_density,
            sa_height = response$sa_height,
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            burrow=response$burrow)

parameters <- c("b0.0","b0.2","b0.3","b0.5","b0.10","b0.11","b0.14","b0.15", 

                "b3.0", "b3.12","b3.13","sigma.sub",
                
                "b5.0","b5.2","b5.2q","b5.14","sigma.sal",
                
                "b10.0","b10.13","b10.14","sigma.sad",
                
                "b11.0","b11.2","b11.12","b11.14","sigma.height",
                
                "b15.0","b15.10","b15.11","b15.12","b15.13","b15.14")

inits <- function(){list(b0.0=0,b0.2=0,b0.3=0,b0.5=0,b0.10=0,b0.11=0,b0.14=0,

                         b5.0=0,b5.2=0,b5.2q=0,b5.14.0=0,sigma.sal=100,
                         
                         b10.0=0,b10.13=0,b10.14=0,sigma.sad=100,
                         
                         b11.0=0,b11.2=0,b11.12=0,b11.14=0,sigma.height=100,
                         
                         b15.0=0,b15.10=0,b15.11=0,b15.12=0,b15.13=0,b15.14=0)}

# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

## RUN MODEL IN JAGS ###########################################################
added_paths_out <- jags(data=data, 
                        inits=inits, 
                        parameters.to.save=parameters, 
                        model.file="recruit.added.paths.model.txt", 
                        n.chains=nc, 
                        n.iter=ni, 
                        n.burnin=nb, 
                        n.thin=nt,  
                        DIC=TRUE, 
                        parallel = T)

print(added_paths_out,digits=4) 


################################################################################
### FINAL PRUNING OF NS PARHWAYS ###
################################################################################

sink("recruit.final.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Recruitment ##	 
	  recruit[i] ~ dpois(mu.recruit[i])                          
		log(mu.recruit[i]) <- b0.0 + 
	  	b0.2*temp[i] + 
	  	b0.5*mean_sal[i] + 
	  	b0.10*sa_density[i] +
	  	b0.11*sa_height[i] +
	  	b0.14*main[i] +
	  	b0.15*burrow[i] 
	 
## Experiment level varaibles ##
#5. sal
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]

#10. Sa_denisty
  sa_density[i] ~dpois(sad.hat[i])
  log(sad.hat[i]) <- b10.0 +
    b10.13*head[i] +
    b10.14*main[i]  # ADDED

#11. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b11.0 +   
    b11.12*ele[i]+
    b11.14*main[i]+  
    b11.2*temp[i] #ADDED

#15. Crab burrow
burrow[i] ~ dpois(burrow.hat[i])
  log(burrow.hat[i]) <- b15.0 +   
    b15.12*ele[i]+
    b15.13*head[i]+
    b15.14*main[i] +
    b15.10*sa_density[i] +  #ADDED
	  b15.11*sa_height[i]  # ADDED
}

	
### PRECISION VARIABLES  #######################################################
	tau.growth <- 1/(sigma.growth*sigma.growth)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.sad <- 1/(sigma.sad*sigma.sad)
	tau.height <- 1/(sigma.height*sigma.height)
	
### PRIORS  ####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001);b0.5 ~ dnorm(0,0.00001)
	b0.10 ~ dnorm(0,0.00001); b0.11 ~ dnorm(0,0.00001); 	b0.14 ~ dnorm(0,0.00001)
	b0.15 ~ dnorm(0,0.00001);sigma.growth ~ dunif(0,100)

  b5.0 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);b5.2q ~ dnorm(0,0.00001)
  b5.14 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)
  
  b10.0 ~ dnorm(0,0.00001); b10.13 ~ dnorm(0,0.00001);b10.14 ~ dnorm(0,0.00001); 
  sigma.sad ~ dunif(0,100);
  
  b11.0 ~ dnorm(0,0.00001); b11.2 ~ dnorm(0,0.00001); b11.12 ~ dnorm(0,0.00001)
  b11.14 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100)

  b15.0 ~ dnorm(0,0.00001); b15.10 ~ dnorm(0,0.00001); b15.11 ~ dnorm(0,0.00001)  
  b15.12 ~ dnorm(0,0.00001); b15.13 ~ dnorm(0,0.00001); b15.14 ~ dnorm(0,0.00001)
}
",fill=TRUE)
sink()

### CREATING OBJECTS TO HAND TO JAGS  ###########################################

N=90  # number of seasons x plots

data = list(N = N,
            recruit = response$recruit,
            gd_density = response$gd_density_max_sum,
            temp = response$mean_water_temp, 
            mean_sal = response$mean_sal,
            sa_density = response$sa_density,
            sa_height = response$sa_height,
            ele = response$ele, 
            head = response$head, 
            main = response$main, 
            burrow=response$burrow)

parameters <- c("b0.0","b0.1","b0.2","b0.5","b0.10","b0.11","b0.14","b0.15", 

                "b5.0","b5.2","b5.2q","b5.14","sigma.sal",
                
                "b10.0","b10.13","b10.14","sigma.sad",
                
                "b11.0","b11.2","b11.12","b11.14","sigma.height",
                
                "b15.0","b15.10","b15.11","b15.12","b15.13","b15.14")

inits <- function(){list(b0.0=0,b0.2=0,b0.5=0,b0.10=0,b0.11=0,b0.14=0,
         
                         b5.0=0,b5.2=0,b5.2q=0,b5.14.0=0,sigma.sal=100,
                         
                         b10.0=0,b10.13=0,b10.14=0,sigma.sad=100,
                         
                         b11.0=0,b11.2=0,b11.12=0,b11.14=0,sigma.height=100,
                         
                         b15.0=0,b15.12=0,b15.13=0,b15.14=0,b15.10=0,b15.11=0)}

# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

## RUN MODEL IN JAGS ###########################################################
out_final <- jags(data=data, 
                  inits=inits, 
                  parameters.to.save=parameters, 
                  model.file="recruit.final.model.txt",      
                  n.chains=nc, 
                  n.iter=ni, 
                  n.burnin=nb, 
                  n.thin=nt,  
                  DIC=TRUE,
                  parallel = T)

# Print some basic results
print(out_final,digits=4)  

## Export results for prediction maps ##
# saveRDS(out_final,"model_out/recuit_SEM_final_output.rds")

## Export results for table ##
# hyp_out_final <- as.data.frame(out_final$summary)  # result summary
# hyp_out_final <- round(hyp_out_final,digits = 3)  # round results
# hyp_out_final$median_95 <- paste(hyp_out_final$`50%`," (",hyp_out_final$`2.5%`,",",hyp_out_final$`97.5%`,")",sep ="")  # create column for median and 95CI
# hyp_out_final <- hyp_out_final %>% select(median_95,Rhat,overlap0)
# write.csv(hyp_out,"~/Documents/School/Master's/Thesis documents/ch2 pub/results/model outputs/recruit_final_out.csv")

### PREDICTION EQUATIONS #######################################################

# Recruitment
recruit.hat <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*response$mean_water_temp +
  out_final$mean$b0.5*response$mean_sal + 
  out_final$mean$b0.10*response$sa_density +
  out_final$mean$b0.11*response$sa_height + 
  out_final$mean$b0.14*response$main +
  out_final$mean$b0.15*response$burrow
recruit.hat_exp <- exp(recruit.hat)  

# Salinity
sal.hat <- out_final$mean$b5.0 + 
  out_final$mean$b5.2*response$mean_water_temp +
  out_final$mean$b5.2q*response$mean_water_temp^2+
  out_final$mean$b5.14*response$main

# SA density
sad.hat <- exp(out_final$mean$b10.0 + 
                 out_final$mean$b10.13*response$head + 
                 out_final$mean$b10.14*response$main)

# SA height
sah.hat <- out_final$mean$b11.0 + 
  out_final$mean$b11.12*response$ele + 
  out_final$mean$b11.14*response$main + 
  out_final$mean$b11.2*response$mean_water_temp

# Burrow
crab.hat <- exp(out_final$mean$b15.0 + 
  out_final$mean$b15.12*response$ele + 
  out_final$mean$b15.13*response$head + 
  out_final$mean$b15.14*response$main+
  out_final$mean$b15.10*response$sa_density +
  out_final$mean$b15.11*response$sa_height)


### RESIDUALS ###################################################################
recruit.res <- response$recruit-recruit.hat_exp
temp.res <- response$mean_water_temp	  
sad.res <- response$sa_density-sad.hat		  
sah.res <- response$sa_height-sah.hat 	  
crab.res <- response$burrow-crab.hat 	  
sal.res <- response$mean_sal - sal.hat 	 
ele.res <- response$ele 	  
main.res <- response$main  
head.res <- response$head	 

# Residual matrix
final.res <- data.frame(recruit.res,temp.res,sad.res,sah.res,crab.res,sal.res,ele.res,main.res,head.res)

# Export growth residuals
resid_out <- data.frame(
  plot = response$plot,
  season = response$season,
  resid = recruit.res
)
saveRDS(resid_out,"model_out/recruit_resid.rds")

###  FINAL MODEL R2 ############################################################
summary(lm(response$recruit ~recruit.hat_exp))
summary(lm(response$mean_sal~sal.hat))
summary(lm(response$sa_density~sad.hat))
summary(lm(response$sa_height~sah.hat))
summary(lm(response$burrow~crab.hat))

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

## Combine both correlation test and find variables that are significantly
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

### FOR RECRUITMENT ############################################################

### TEMP ###

# min
recruit.predict.2 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*min(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow)

# max
recruit.predict.3 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*max(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow)

# using straight probabilities
recruit.predict.2pr.raw <- exp(recruit.predict.2)
recruit.predict.3pr.raw <- exp(recruit.predict.3)

# QUERY
(recruit.predict.3pr.raw-recruit.predict.2pr.raw)/(max(response$recruit)-min(response$recruit))

### SALINITY ###

# min
recruit.predict.4 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*min(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow)

# max
recruit.predict.5 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*max(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow)

# using straight probabilities
recruit.predict.4pr.raw <- exp(recruit.predict.4)
recruit.predict.5pr.raw <- exp(recruit.predict.5)

# QUERY
(recruit.predict.5pr.raw-recruit.predict.4pr.raw)/(max(response$recruit)-min(response$recruit))

### SA DENSITY ###

# min
recruit.predict.6 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*min(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow)

# max
recruit.predict.7 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*max(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow); 			      

# using straight probabilities
recruit.predict.6pr.raw <- exp(recruit.predict.6)
recruit.predict.7pr.raw <- exp(recruit.predict.7)
# QUERY
(recruit.predict.7pr.raw-recruit.predict.6pr.raw)/(max(response$recruit)-min(response$recruit))

### SA HEIGHT ###

# min
recruit.predict.8 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*min(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow)

# max
recruit.predict.9 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*max(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*mean(response$burrow); 			      

# using straight probabilities
recruit.predict.8pr.raw <- exp(recruit.predict.8)
recruit.predict.9pr.raw <- exp(recruit.predict.9)

# QUERY
(recruit.predict.9pr.raw-recruit.predict.8pr.raw)/(max(response$recruit)-min(response$recruit))

### MAIN ###

# min
recruit.predict.10 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*min(response$main) +
  out_final$mean$b0.15*mean(response$burrow)

# max
recruit.predict.11 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*max(response$main) +
  out_final$mean$b0.15*mean(response$burrow); 			      

# using straight probabilities
recruit.predict.10pr.raw <- exp(recruit.predict.10)
recruit.predict.11pr.raw <- exp(recruit.predict.11)

# QUERY
(recruit.predict.11pr.raw-recruit.predict.10pr.raw)/(max(response$recruit)-min(response$recruit))

### BURROWS ###

# min
recruit.predict.12 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*min(response$burrow)

# max
recruit.predict.13 <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp) +
  out_final$mean$b0.5*mean(response$mean_sal) + 
  out_final$mean$b0.10*mean(response$sa_density) +
  out_final$mean$b0.11*mean(response$sa_height) + 
  out_final$mean$b0.14*mean(response$main) +
  out_final$mean$b0.15*max(response$burrow); 			      

# using straight probabilities
recruit.predict.12pr.raw <- exp(recruit.predict.12)
recruit.predict.13pr.raw <- exp(recruit.predict.13)

# QUERY
(recruit.predict.13pr.raw-recruit.predict.12pr.raw)/(max(response$recruit)-min(response$recruit))

### FOR SALINITY ###############################################################

### TEMP ###

# min
sal.predict.2 <- out_final$mean$b5.0 + 
  out_final$mean$b5.2*min(response$mean_water_temp,na.rm = T) + 
  out_final$mean$b5.2q*min(response$mean_water_temp,na.rm = T)^2+
  out_final$mean$b5.14*mean(response$main)

# max
sal.predict.3 <- out_final$mean$b5.0 + 
  out_final$mean$b5.2*max(response$mean_water_temp,na.rm = T) + 
  out_final$mean$b5.2q*max(response$mean_water_temp,na.rm = T)^2+
  out_final$mean$b5.14*mean(response$main)

# QUERY
(sal.predict.3-sal.predict.2)/(max(response$mean_sal)-min(response$mean_sal))

### Main ###

# min
sal.predict.6 <- out_final$mean$b5.0 + 
  out_final$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
  out_final$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
  out_final$mean$b5.14*min(response$main)

# max
sal.predict.7 <- out_final$mean$b5.0 + 
  out_final$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
  out_final$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
  out_final$mean$b5.14*max(response$main)

# QUERY
(sal.predict.7-sal.predict.6)/(max(response$mean_sal)-min(response$mean_sal))

### FOR SA DENSITY #############################################################

### HEAD ###

# min
sad.predict.2 <- out_final$mean$b10.0 + 
  out_final$mean$b10.13*min(response$head) + 
  out_final$mean$b10.14*mean(response$main)		      

# max
sad.predict.3 <- out_final$mean$b10.0 + 
  out_final$mean$b10.13*max(response$head) + 
  out_final$mean$b10.14*mean(response$main)			      

# using straight probabilities
sad.predict.2pr.raw <- exp(sad.predict.2)
sad.predict.3pr.raw <- exp(sad.predict.3)

# QUERY
(sad.predict.3pr.raw-sad.predict.2pr.raw)/(max(response$sa_density)-min(response$sa_density))

### MAIN ###

# Min
sad.predict.4 <- out_final$mean$b10.0 + 
  out_final$mean$b10.13*mean(response$head) + 
  out_final$mean$b10.14*min(response$main)		      

# Max
sad.predict.5 <- out_final$mean$b10.0 + 
  out_final$mean$b10.13*mean(response$head) + 
  out_final$mean$b10.14*max(response$main) 			      

# using straight probabilities
sad.predict.4pr.raw <- exp(sad.predict.4)
sad.predict.5pr.raw <- exp(sad.predict.5)

# QUERY
(sad.predict.5pr.raw-sad.predict.4pr.raw)/(max(response$sa_density)-min(response$sa_density))

### FOR SA HEIGHT ##############################################################

### ELE ###

# min
sah.predict.2 <- out_final$mean$b11.0 + 
  out_final$mean$b11.12*min(response$ele) + 
  out_final$mean$b11.14*mean(response$main) + 
  out_final$mean$b11.2*mean(response$mean_water_temp)

# max
sah.predict.3 <- out_final$mean$b11.0 + 
  out_final$mean$b11.12*max(response$ele) + 
  out_final$mean$b11.14*mean(response$main) + 
  out_final$mean$b11.2*mean(response$mean_water_temp)

# QUERY
(sah.predict.3-sah.predict.2)/(max(response$sa_height)-min(response$sa_height))

### MAIN ###

# min
sah.predict.4 <- out_final$mean$b11.0 + 
  out_final$mean$b11.12*mean(response$ele) + 
  out_final$mean$b11.14*min(response$main) + 
  out_final$mean$b11.2*mean(response$mean_water_temp)		       

# max
sah.predict.5 <- out_final$mean$b11.0 + 
  out_final$mean$b11.12*mean(response$ele) + 
  out_final$mean$b11.14*max(response$main) + 
  out_final$mean$b11.2*mean(response$mean_water_temp)   	

# QUERY
(sah.predict.5-sah.predict.4)/(max(response$sa_height)-min(response$sa_height))

### TEMP ###

# min
sah.predict.6 <- out_final$mean$b11.0 + 
  out_final$mean$b11.12*mean(response$ele) + 
  out_final$mean$b11.14*mean(response$main) + 
  out_final$mean$b11.2*min(response$mean_water_temp)

# max
sah.predict.7 <- out_final$mean$b11.0 + 
  out_final$mean$b11.12*mean(response$ele) + 
  out_final$mean$b11.14*mean(response$main) + 
  out_final$mean$b11.2*max(response$mean_water_temp)

# QUERY
(sah.predict.7-sah.predict.6)/(max(response$sa_height)-min(response$sa_height))

### FOR BURROW ##################################################################

### ELE ###

# min
crab.predict.2 <- exp(out_final$mean$b15.0 + 
  out_final$mean$b15.12*min(response$ele) + 
  out_final$mean$b15.13*mean(response$head) + 
  out_final$mean$b15.14*mean(response$main)+
  out_final$mean$b15.10*mean(response$sa_density) +
  out_final$mean$b15.11*mean(response$sa_height)) 

# max
crab.predict.3 <- exp(out_final$mean$b15.0 + 
  out_final$mean$b15.12*max(response$ele) + 
  out_final$mean$b15.13*mean(response$head) + 
  out_final$mean$b15.14*mean(response$main)+
  out_final$mean$b15.10*mean(response$sa_density) +
  out_final$mean$b15.11*mean(response$sa_height)) 			      

# QUERY
(crab.predict.3-crab.predict.2)/(max(response$burrow)-min(response$burrow))

### HEAD ###

# min
crab.predict.4 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*min(response$head) + 
                        out_final$mean$b15.14*mean(response$main)+
                        out_final$mean$b15.10*mean(response$sa_density) +
                        out_final$mean$b15.11*mean(response$sa_height)) 	

# max
crab.predict.5 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*max(response$head) + 
                        out_final$mean$b15.14*mean(response$main)+
                        out_final$mean$b15.10*mean(response$sa_density) +
                        out_final$mean$b15.11*mean(response$sa_height))  			      

# QUERY
(crab.predict.5-crab.predict.4)/(max(response$burrow)-min(response$burrow))

### MAIN ###

# min
crab.predict.6 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*mean(response$head) + 
                        out_final$mean$b15.14*min(response$main)+
                        out_final$mean$b15.10*mean(response$sa_density) +
                        out_final$mean$b15.11*mean(response$sa_height))  

# max
crab.predict.7 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*mean(response$head) + 
                        out_final$mean$b15.14*max(response$main)+
                        out_final$mean$b15.10*mean(response$sa_density) +
                        out_final$mean$b15.11*mean(response$sa_height)) 			      

# QUERY
(crab.predict.7-crab.predict.6)/(max(response$burrow)-min(response$burrow))

### SA DENSITY ###

# min
crab.predict.8 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*mean(response$head) + 
                        out_final$mean$b15.14*mean(response$main)+
                        out_final$mean$b15.10*min(response$sa_density) +
                        out_final$mean$b15.11*mean(response$sa_height)) 	

# max
crab.predict.9 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*mean(response$head) + 
                        out_final$mean$b15.14*mean(response$main)+
                        out_final$mean$b15.10*max(response$sa_density) +
                        out_final$mean$b15.11*mean(response$sa_height)) 			      
# QUERY
(crab.predict.9-crab.predict.8)/(max(response$burrow)-min(response$burrow))

### SA HEIGHT ###

# min
crab.predict.10 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*mean(response$head) + 
                        out_final$mean$b15.14*mean(response$main)+
                        out_final$mean$b15.10*mean(response$sa_density) +
                        out_final$mean$b15.11*min(response$sa_height)) 

# max
crab.predict.11 <- exp(out_final$mean$b15.0 + 
                        out_final$mean$b15.12*mean(response$ele) + 
                        out_final$mean$b15.13*mean(response$head) + 
                        out_final$mean$b15.14*mean(response$main)+
                        out_final$mean$b15.10*mean(response$sa_density) +
                        out_final$mean$b15.11*max(response$sa_height))  			      

# QUERY
(crab.predict.11-crab.predict.10)/(max(response$burrow)-min(response$burrow))

################################################################################
## TOTAL EFFECTS QUERIES ##
################################################################################

### TEMP ###

# min
recruit.predict.tot2  <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*min(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*(out_final$mean$b5.0 + 
                         out_final$mean$b5.2*min(response$mean_water_temp,na.rm = T) + 
                         out_final$mean$b5.2q*min(response$mean_water_temp,na.rm = T)^2+
                         out_final$mean$b5.14*mean(response$main,na.rm = T)) + 
  out_final$mean$b0.10*mean(response$sa_density,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 + 
                          out_final$mean$b11.12*mean(response$ele,na.rm = T) + 
                          out_final$mean$b11.14*mean(response$main,na.rm = T) + 
                          out_final$mean$b11.2*min(response$mean_water_temp,na.rm = T)) + 
  out_final$mean$b0.14*mean(response$main,na.rm = T) +
  out_final$mean$b0.15*mean(response$burrow,na.rm = T) 
recruit.predict.tot2pr.raw <- exp(recruit.predict.tot2)

# max
recruit.predict.tot3  <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*max(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*(out_final$mean$b5.0 + 
                         out_final$mean$b5.2*max(response$mean_water_temp,na.rm = T) + 
                         out_final$mean$b5.2q*max(response$mean_water_temp,na.rm = T)^2+
                         out_final$mean$b5.14*mean(response$main,na.rm = T)) + 
  out_final$mean$b0.10*mean(response$sa_density,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 + 
                          out_final$mean$b11.12*mean(response$ele,na.rm = T) + 
                          out_final$mean$b11.14*mean(response$main,na.rm = T) + 
                          out_final$mean$b11.2*max(response$mean_water_temp,na.rm = T)) + 
  out_final$mean$b0.14*mean(response$main,na.rm = T) +
  out_final$mean$b0.15*mean(response$burrow,na.rm = T) 
recruit.predict.tot3pr.raw <- exp(recruit.predict.tot3)

# QUERY
(recruit.predict.tot3pr.raw-recruit.predict.tot2pr.raw)/(max(response$recruit)-min(response$recruit))

## ELE ###

# min
recruit.predict.tot2 <- recruit.hat <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*mean(response$mean_sal,na.rm = T) + 
  out_final$mean$b0.10*mean(response$sa_density,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 + 
      out_final$mean$b11.12*min(response$ele,na.rm = T) + 
      out_final$mean$b11.14*mean(response$main,na.rm = T) + 
      out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T)) + 
  out_final$mean$b0.14*mean(response$main,na.rm = T) +
  out_final$mean$b0.15*(exp(out_final$mean$b15.0 + 
      out_final$mean$b15.12*min(response$ele) + 
      out_final$mean$b15.13*mean(response$head) + 
      out_final$mean$b15.14*mean(response$main)+
      out_final$mean$b15.10*mean(response$sa_density) +
      out_final$mean$b15.11*(out_final$mean$b11.0 + 
             out_final$mean$b11.12*min(response$ele,na.rm = T) + 
             out_final$mean$b11.14*mean(response$main,na.rm = T) + 
             out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T))))
recruit.predict.tot2pr.raw <- exp(recruit.predict.tot2)

# max
recruit.predict.tot3 <- recruit.hat <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*mean(response$mean_sal,na.rm = T) + 
  out_final$mean$b0.10*mean(response$sa_density,na.rm = T) +
  out_final$mean$b0.11*(out_final$mean$b11.0 + 
      out_final$mean$b11.12*max(response$ele,na.rm = T) + 
      out_final$mean$b11.14*mean(response$main,na.rm = T) + 
      out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T)) + 
  out_final$mean$b0.14*mean(response$main,na.rm = T) +
  out_final$mean$b0.15*(exp(out_final$mean$b15.0 + 
      out_final$mean$b15.12*max(response$ele) + 
      out_final$mean$b15.13*mean(response$head) + 
      out_final$mean$b15.14*mean(response$main)+
      out_final$mean$b15.10*mean(response$sa_density) +
      out_final$mean$b15.11*(out_final$mean$b11.0 + 
           out_final$mean$b11.12*max(response$ele,na.rm = T) + 
           out_final$mean$b11.14*mean(response$main,na.rm = T) + 
           out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T)))) 
recruit.predict.tot3pr.raw <- exp(recruit.predict.tot3)

# QUERY
(recruit.predict.tot3pr.raw-recruit.predict.tot2pr.raw)/(max(response$recruit)-min(response$recruit))

### HEAD ###

# min
recruit.predict.tot2 <- recruit.hat <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*mean(response$mean_sal,na.rm = T) + 
  out_final$mean$b0.10*(exp(out_final$mean$b10.0 + 
      out_final$mean$b10.13*min(response$head) + 
      out_final$mean$b10.14*mean(response$main))) +
  out_final$mean$b0.11*mean(response$sa_height,na.rm = T) + 
  out_final$mean$b0.14*mean(response$main,na.rm = T) +
  out_final$mean$b0.15*(exp(out_final$mean$b15.0 + 
      out_final$mean$b15.12*mean(response$ele) + 
      out_final$mean$b15.13*min(response$head) + 
      out_final$mean$b15.14*mean(response$main)+
      out_final$mean$b15.10*(exp(out_final$mean$b10.0 + 
         out_final$mean$b10.13*min(response$head) + 
         out_final$mean$b10.14*mean(response$main))) +
      out_final$mean$b15.11*mean(response$sa_height))) 
recruit.predict.tot2pr.raw <- exp(recruit.predict.tot2)

# max
recruit.predict.tot3 <- recruit.hat <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*mean(response$mean_sal,na.rm = T) + 
  out_final$mean$b0.10*(exp(out_final$mean$b10.0 + 
      out_final$mean$b10.13*max(response$head) + 
      out_final$mean$b10.14*mean(response$main))) +
  out_final$mean$b0.11*mean(response$sa_height,na.rm = T) + 
  out_final$mean$b0.14*mean(response$main,na.rm = T) +
  out_final$mean$b0.15*(exp(out_final$mean$b15.0 + 
      out_final$mean$b15.12*mean(response$ele) + 
      out_final$mean$b15.13*max(response$head) + 
      out_final$mean$b15.14*mean(response$main)+
      out_final$mean$b15.10*(exp(out_final$mean$b10.0 + 
         out_final$mean$b10.13*max(response$head) + 
         out_final$mean$b10.14*mean(response$main))) +
      out_final$mean$b15.11*mean(response$sa_height))) 
recruit.predict.tot3pr.raw <- exp(recruit.predict.tot3)

# QUERY
(recruit.predict.tot3pr.raw-recruit.predict.tot2pr.raw)/(max(response$recruit)-min(response$recruit))

### Main ###

# min
recruit.predict.tot2 <- recruit.hat <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*(out_final$mean$b5.0 + 
       out_final$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
       out_final$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
       out_final$mean$b5.14*min(response$main,na.rm = T)) + 
  out_final$mean$b0.10*(exp(out_final$mean$b10.0 + 
      out_final$mean$b10.13*mean(response$head) + 
      out_final$mean$b10.14*min(response$main))) +
  out_final$mean$b0.11*(out_final$mean$b11.0 + 
      out_final$mean$b11.12*mean(response$ele,na.rm = T) + 
      out_final$mean$b11.14*min(response$main,na.rm = T) + 
      out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T)) + 
  out_final$mean$b0.14*min(response$main,na.rm = T) +
  out_final$mean$b0.15*(exp(out_final$mean$b15.0 + 
      out_final$mean$b15.12*mean(response$ele) + 
      out_final$mean$b15.13*mean(response$head) + 
      out_final$mean$b15.14*min(response$main)+
      out_final$mean$b15.10*(exp(out_final$mean$b10.0 + 
           out_final$mean$b10.13*mean(response$head) + 
           out_final$mean$b10.14*min(response$main))) +
      out_final$mean$b15.11*(out_final$mean$b11.0 + 
           out_final$mean$b11.12*mean(response$ele,na.rm = T) + 
           out_final$mean$b11.14*min(response$main,na.rm = T) + 
           out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T)))) 
recruit.predict.tot2pr.raw <- exp(recruit.predict.tot2)

# max
recruit.predict.tot3 <- recruit.hat <- out_final$mean$b0.0 + 
  out_final$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  out_final$mean$b0.5*(out_final$mean$b5.0 + 
       out_final$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
       out_final$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
       out_final$mean$b5.14*max(response$main,na.rm = T)) + 
  out_final$mean$b0.10*(exp(out_final$mean$b10.0 + 
       out_final$mean$b10.13*mean(response$head) + 
       out_final$mean$b10.14*max(response$main))) +
  out_final$mean$b0.11*(out_final$mean$b11.0 + 
        out_final$mean$b11.12*mean(response$ele,na.rm = T) + 
        out_final$mean$b11.14*max(response$main,na.rm = T) + 
        out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T)) + 
  out_final$mean$b0.14*max(response$main,na.rm = T) +
  out_final$mean$b0.15*(exp(out_final$mean$b15.0 + 
        out_final$mean$b15.12*mean(response$ele) + 
        out_final$mean$b15.13*mean(response$head) + 
        out_final$mean$b15.14*max(response$main)+
        out_final$mean$b15.10*(exp(out_final$mean$b10.0 + 
             out_final$mean$b10.13*mean(response$head) + 
             out_final$mean$b10.14*max(response$main))) +
        out_final$mean$b15.11*(out_final$mean$b11.0 + 
             out_final$mean$b11.12*mean(response$ele,na.rm = T) + 
             out_final$mean$b11.14*max(response$main,na.rm = T) + 
             out_final$mean$b11.2*mean(response$mean_water_temp,na.rm = T))))
recruit.predict.tot3pr.raw <- exp(recruit.predict.tot3)

# QUERY
(recruit.predict.tot3pr.raw-recruit.predict.tot2pr.raw)/(max(response$recruit)-min(response$recruit))

