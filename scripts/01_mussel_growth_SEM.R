################################################################################
#
#  RIBBED MUSSEL GROWTH STRUCTURAL EQUATION MODELS        
#  
#  Created by Anonymized                              
#  Modified from Grace et al., 2012 - Ecosphere           
#
################################################################################

# This code performs the structural equation modelling (SEM) of the hierarchical
# drivers of ribbed mussel growth rates across the salt marsh landscape. Here,
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

## Growth and explanatory data  ##
response <- read_excel("mussel_population_vital_data.xlsx",sheet = 1)

#Order by experiment
response <- response[order(response$season),]

## Check for correlation among explanatory variables  ##
response_cov <- response[,3:27]
response_TF_cov <- abs(cor(response_cov, use = "complete.obs"))>.60

## Check for missing data ##

# what columns have missing data?
colnames(response)[colSums(is.na(response)) > 0]  

# Number of missing data per column
length(response$MaxTemp_ibut[is.na(response$MaxTemp_ibut)])
length(response$mean_bact[is.na(response$mean_bact)])
length(response$mean_chla[is.na(response$mean_chla)])
length(response$mean_om[is.na(response$mean_om)])
length(response$burrow[is.na(response$burrow)])

################################################################################
### VARIABLE CODES FOR MODELS ###
################################################################################

#1. gd_density
#2. temp
#3. sub_prop
#4. flow_rate
#5. mean_sal
#6. mean_bact
#7. mean_chla
#8. mean_om
#9. mussel_size
#10. sa_density
#11. sa_height
#12. ele
#13. head
#14. main

################################################################################
### Global SEM MODEL  ###
################################################################################

sink("growth.global.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Growth	## 
	  growth[i] ~ dnorm(growth.hat[i],tau.growth)                            
		growth.hat[i] <- b0.0 + 
	  	b0.1*gd_density[i] + 
	  	b0.2*temp[i] + b0.2q*temp[i]*temp[i] +
	  	b0.3*sub_prop[i] + 
	  	b0.4*flow_rate[i] + 
	  	b0.5*mean_sal[i] + 
	  	b0.6*mean_bact[i] + 
	  	b0.7*mean_chla[i] + 
	  	b0.8*mean_om[i]+
	  	b0.9*mussel_size[i] 

## Experiment level variables ##
#1. Mussel density
  gd_density[i] ~ dpois(gd.hat[i])
  log(gd.hat[i])<- b1.0 +
    b1.10*sa_density[i] +
    b1.12*ele[i]+
    b1.13*head[i]+
    b1.14*main[i]

#2. Temp
  temp[i] ~ dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0  + 
    b2.3*sub_prop[i] + 
      b2.10*sa_density[i] +
      b2.11*sa_height[i] +
      b2.12*ele[i]

#5. Salinity
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]

#6. bact 
  mean_bact[i] ~dnorm(bact.hat[i],tau.bact)
    bact.hat[i] <-b6.0 +
    b6.4*flow_rate[i] +
    b6.12*ele[i]+
    b6.13*head[i] + 
    b6.14*main[i] 

#7. chla 
  mean_chla[i] ~dnorm(chla.hat[i],tau.chla)
    chla.hat[i] <-b7.0 +
    b7.4*flow_rate[i] +
    b7.12*ele[i] +
    b7.13*head[i] +
    b7.14*main[i] 

#8. om 
  mean_om[i] ~dnorm(om.hat[i],tau.om)
    om.hat[i] <-b8.0 +
    b8.4*flow_rate[i] +
    b8.12*ele[i] +
    b8.13*head[i] +
    b8.14*main[i] 

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
	}

## Plot level variables ##	
for(j in 1:30) {
#3. Submergence Time
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
	tau.temp <- 1/(sigma.temp*sigma.temp)
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.flow <- 1/(sigma.flow*sigma.flow)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.height <- 1/(sigma.height*sigma.height)
	tau.bact <- 1/(sigma.bact*sigma.bact)
	tau.chla <- 1/(sigma.chla*sigma.chla)
	tau.om <- 1/(sigma.om*sigma.om)

### PRIORS #####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.1 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001) 
	b0.2q ~ dnorm(0,0.00001); b0.3 ~ dnorm(0,0.00001); b0.4 ~ dnorm(0,0.00001)
	b0.5 ~ dnorm(0,0.00001); b0.6 ~ dnorm(0,0.00001); b0.7 ~ dnorm(0,0.00001) 
	b0.8 ~ dnorm(0,0.00001);b0.9 ~ dnorm(0,0.00001); sigma.growth ~ dunif(0,100)

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
  
  b6.0 ~ dnorm(0,0.00001); b6.4 ~ dnorm(0,0.00001);b6.12 ~ dnorm(0,0.00001)
  b6.13 ~ dnorm(0,0.00001);b6.14 ~ dnorm(0,0.00001); sigma.bact ~ dunif(0,100)
  
  b7.0 ~ dnorm(0,0.00001); b7.4 ~ dnorm(0,0.00001);b7.12 ~ dnorm(0,0.00001)
  b7.13 ~ dnorm(0,0.00001);b7.14 ~ dnorm(0,0.00001); sigma.chla ~ dunif(0,100)
  
  b8.0 ~ dnorm(0,0.00001); b8.4 ~ dnorm(0,0.00001);b8.12 ~ dnorm(0,0.00001)
  b8.13 ~ dnorm(0,0.00001);b8.14 ~ dnorm(0,0.00001); sigma.om ~ dunif(0,100)
  
  b10.0 ~ dnorm(0,0.00001);  b10.12 ~ dnorm(0,0.00001);b10.13 ~ dnorm(0,0.00001) 
  b10.14 ~ dnorm(0,0.00001); sigma.sad ~ dunif(0,100);
  
  b11.0 ~ dnorm(0,0.00001);  b11.12 ~ dnorm(0,0.00001);b11.13 ~ dnorm(0,0.00001)
  b11.14 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100);
}
    ",fill=TRUE)
sink()
# end of JAGS code creation

### CREATE OBJECTS TO HAND OFF TO JAGS  #########################################

N=90 # Number of plots*seasons

data = list(N = N,
            growth = response$growth_rate,
            gd_density = response$gd_density_max_sum,
            temp = response$mean_water_temp, 
            sub_prop = response$sub_prop, 
            flow_rate = response$flow_rate, 
            mean_sal = response$mean_sal,
            mean_bact = response$mean_bact, 
            mean_chla = response$mean_chla, 
            mean_om = response$mean_om,
            mussel_size = response$mussel_size, 
            sa_density = response$sa_density,
            sa_density_m = response$sa_denisty_m,
            sa_height = response$sa_height,
            sa_height_m = response$sa_height_m,
            ele = response$ele, 
            head = response$head, 
            main = response$main) 

parameters <- c("b0.0","b0.1","b0.2","b0.2q","b0.3","b0.4","b0.5","b0.6",                          
                "b0.7","b0.8","b0.9","sigma.growth",
                
                "b1.0","b1.10","b1.12","b1.13","b1.14",
                
                "b2.0","b2.3","b2.10","b2.11","b2.12", "sigma.temp",
                
                "b3.0", "b3.12","b3.13","sigma.sub",
                
                "b4.0","b4.10","b4.11","b4.12","b4.13","sigma.flow",
                
                "b5.0","b5.2","b5.2q","b5.14","sigma.sal",
                
                "b6.0", "b6.4","b6.12","b6.13", "b6.14", "sigma.bact",
                
                "b7.0", "b7.4","b7.12","b7.13", "b7.14", "sigma.chla",
                
                "b8.0", "b8.4","b8.12","b8.13", "b8.14", "sigma.om",
                
                "b10.0","b10.12","b10.13","b10.14",
                
                "b11.0","b11.10","b11.12","b11.13","b11.14","sigma.height")

inits <- function(){list(b0.0=0,b0.1=0,b0.2=0,b0.2q=0,b0.3=0,b0.4=0,b0.5=0,                                
                         b0.6=0,b0.7=0,b0.8=0,b0.9=0,sigma.growth=100, 
                         
                         b1.0=0,b1.10=0,b1.12=0,b1.13=0,b1.14=0,sigma.gd=100,
                      
                         b2.0=0, b2.3=0,b2.10=0,b2.11=0,b2.12=0, sigma.temp=100,
                         
                         b3.0=0, b3.12=0, b3.13=0, sigma.sub=100,
                         
                         b4.0=0,b4.10=0,b4.11=0,b4.12=0,b4.13=0,sigma.flow=100,
                         
                         b5.0=0,b5.2=0,b5.2q=0,b5.14=0,sigma.sal=100,
                         
                         b6.0=0,b6.4=0,b6.12=0,b6.13=0,b6.14=0,sigma.bact=100,
                         
                         b7.0=0,b7.4=0,b7.12=0,b7.13=0,b7.14=0,sigma.chla=100,
                         
                         b8.0=0,b8.4=0,b8.12=0,b8.13=0,b8.14=0,sigma.om=100,
                         
                         b10.0=0,b10.12=0,b10.13=0,b10.14=0,sigma.sad=100,
                         
                         b11.0=0,b11.12=0,b11.13=0,b11.14=0,sigma.height=100,
                         
                         b15.0=0,b15.3=0,b15.13=0,b15.14=0,sigma.pH=100)}

# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain   
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

### RUN MODEL IN JAGS  ##########################################################
global.out <- jags(data=data, 
                            inits=inits, 
                            parameters.to.save=parameters, 
                            model.file="growth.global.model.txt", 
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
# hyp_out <- round(hyp_out,digits = 3)  # round results
# hyp_out$median_95 <- paste(hyp_out$`50%`," (",hyp_out$`2.5%`,",",hyp_out$`97.5%`,")",sep ="")  # create column with medians and 95CI
# hyp_out <- hyp_out %>% select(median_95,Rhat,overlap0)
# write.csv(hyp_out,"~/Documents/School/Master's/Thesis documents/ch2 pub/results/model outputs/growth_hyp_out_water_temp.csv")

### PREDICTION EQUATIONS ########################################################

# Growth
growth.hat <- global.out$mean$b0.0 + 
  global.out$mean$b0.1*response$gd_density_max_sum +
  global.out$mean$b0.2*response$mean_water_temp + global.out$mean$b0.2q*response$mean_water_temp^2 +
  global.out$mean$b0.3*response$sub_prop + 
  global.out$mean$b0.4*response$flow_rate +
  global.out$mean$b0.5*response$mean_sal + 
  global.out$mean$b0.6*response$mean_bact +
  global.out$mean$b0.7*response$mean_chla +
  global.out$mean$b0.8*response$mean_om +
  global.out$mean$b0.9*response$mussel_size

# GD density
gd.hat <- exp(global.out$mean$b1.0 +
  global.out$mean$b1.10*response$sa_density +
  global.out$mean$b1.12*response$ele +
  global.out$mean$b1.13*response$head +
  global.out$mean$b1.14*response$main)

# Temperature
temp.hat <- global.out$mean$b2.0 +
  global.out$mean$b2.3*response$sub_prop +
  global.out$mean$b2.10*response$sa_density +
  global.out$mean$b2.11*response$sa_height +
  global.out$mean$b2.12*response$ele

# Submergance
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
  global.out$mean$b5.2*response$mean_water_temp+ global.out$mean$b5.2q*response$mean_water_temp^2+
  global.out$mean$b5.14*response$main

# Bacteria
bact.hat <- global.out$mean$b6.0 +
  global.out$mean$b6.4*response$flow_rate + 
  global.out$mean$b6.12*response$ele + 
  global.out$mean$b6.13*response$head +
  global.out$mean$b6.14*response$main

# Chla
chla.hat <- global.out$mean$b7.0 +
  global.out$mean$b7.4*response$flow_rate + 
  global.out$mean$b7.12*response$ele + 
  global.out$mean$b7.13*response$head +
  global.out$mean$b7.14*response$main

# POM
pom.hat <- global.out$mean$b8.0 +
  global.out$mean$b8.4*response$flow_rate + 
  global.out$mean$b8.12*response$ele + 
  global.out$mean$b8.13*response$head +
  global.out$mean$b8.14*response$main

# SA Density
sad.hat <- exp(global.out$mean$b10.0 + 
  global.out$mean$b10.12*response$ele + 
  global.out$mean$b10.13*response$head)

# SA Height
sah.hat <- global.out$mean$b11.0 + 
  global.out$mean$b11.12*response$ele + 
  global.out$mean$b11.13*response$head +
  global.out$mean$b11.14*response$main

### MODEL R2 VALUES #############################################################
summary(lm(response$growth_rate~growth.hat))
summary(lm(response$gd_density_max_sum~gd.hat))
summary(lm(response$mean_water_temp~temp.hat))
summary(lm(response$sub_prop~sub.hat))
summary(lm(response$flow~flow.hat))
summary(lm(response$mean_sal~sal.hat))
summary(lm(response$mean_bact~bact.hat))
summary(lm(response$mean_chla~chla.hat))
summary(lm(response$mean_om~pom.hat))
summary(lm(response$sa_density~sad.hat))
summary(lm(response$sa_height~sah.hat))

################################################################################
### GLOBAL MODEL AFTER PRUNING NS PATHWAYS  ###
################################################################################

sink("growth.global.pruned.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Growth	 ##
	  growth[i] ~ dnorm(growth.hat[i],tau.growth)                            
		growth.hat[i] <- b0.0 + 
	  	b0.2*temp[i] +
	  	b0.3*sub_prop[i] + 
	  	b0.4*flow_rate[i] + 
	  	b0.5*mean_sal[i] + 
	  	b0.9*mussel_size[i] 

## Experiment level variables ##	  	
#2. Temp
  temp[i] ~ dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0  

#5. Salinity
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]

#11. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b11.0 +   
    b11.12*ele[i]+
    b11.14*main[i]

	}
	
## Plot level variables  ##
for(j in 1:30) {
#3. Submergence Time
  sub_prop[j] ~ dnorm(sub.hat[j], tau.sub)
  sub.hat[j] <- b3.0 + 
   b3.12*ele[j] + 
   b3.13*head[j]

#4. Flow
  flow_rate[j] ~ dnorm(flow.hat[j], tau.flow)
    flow.hat[j] <- b4.0 + 
      b4.11*sa_height_m[j]

	}
	
### PRECISION VARIABLES  #######################################################
	tau.growth <- 1/(sigma.growth*sigma.growth)
	tau.temp <- 1/(sigma.temp*sigma.temp)
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.flow <- 1/(sigma.flow*sigma.flow)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.height <- 1/(sigma.height*sigma.height)

### PRIORS #####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001); b0.3 ~ dnorm(0,0.00001)
	b0.4 ~ dnorm(0,0.00001); b0.5 ~ dnorm(0,0.00001); b0.9 ~ dnorm(0,0.00001) 
	sigma.growth ~ dunif(0,100)

	b2.0 ~ dnorm(0,0.00001); sigma.temp ~ dunif(0,100)
	
  b3.0 ~ dnorm(0,0.00001); b3.12 ~ dnorm(0,0.00001); b3.13 ~ dnorm(0,0.00001) 
  sigma.sub ~ dunif(0,100)
  
  b4.0 ~ dnorm(0,0.00001); b4.11 ~ dnorm(0,0.00001) ; sigma.flow ~ dunif(0,100)

  b5.0 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);b5.2q ~ dnorm(0,0.00001)
  b5.14 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)
  
  b11.0 ~ dnorm(0,0.00001);  b11.12 ~ dnorm(0,0.00001); b11.14 ~ dnorm(0,0.00001) 
  sigma.height ~ dunif(0,100);

}
    ",fill=TRUE)
sink()

# end of jags code creation

### CREATE OBJECTS TO HAND OFF TO JAGS  #########################################

N=90 # plots x seasons
  
data = list(N = N,
            growth = response$growth_rate,
            temp = response$mean_water_temp, 
            sub_prop = response$sub_prop, 
            flow_rate = response$flow_rate, 
            mean_sal = response$mean_sal,
            mussel_size = response$mussel_size, 
            sa_height = response$sa_height,
            sa_height_m = response$sa_height_m,
            ele = response$ele, 
            head = response$head, 
            main = response$main) 

parameters <- c("b0.0","b0.2","b0.3","b0.4","b0.5","b0.9","sigma.growth",
                
                "b2.0",
                
                "b3.0", "b3.12","b3.13","sigma.sub",
                
                "b4.0","b4.11","sigma.flow",
                
                "b5.0","b5.2","b5.2q","b5.14","sigma.sal",
                
                "b11.0","b11.10","b11.12","b11.14","sigma.height")

inits <- function(){list(b0.0=0,b0.1=0,b0.2=0,b0.2q=0,b0.3=0,b0.4=0,b0.5=0,                                
                         b0.9=0,sigma.growth=100, 
                         
                         b2.0=0, sigma.temp=100,
                         
                         b3.0=0, b3.12=0, b3.13=0, sigma.sub=100,
                         
                         b4.0=0,b4.11=0,sigma.flow=100,
                         
                         b5.0=0,b5.2=0,b5.2q=0,b5.14=0,sigma.sal=100,
                         
                         b11.0=0,b11.12=0,b11.13=0,b11.14=0,sigma.height=100)}


# MCMC settings to hand to Jags
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain   
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

### RUN MODEL IN JAGS  ##########################################################
out_ns <- jags(data=data, inits=inits, 
               parameters.to.save=parameters,
               model.file="growth.global.pruned.model.txt", 
               n.chains=nc, 
               n.iter=ni, 
               n.burnin=nb, 
               n.thin=nt, 
               DIC=TRUE, 
               parallel = T)

# Print some basic results
print(out_ns,digits=4)  

###  PREDICTION EQUATIONS ######################################################

# Growth
growth.hat <- out_ns$mean$b0.0 + 
  out_ns$mean$b0.2*response$mean_water_temp +
  out_ns$mean$b0.3*response$sub_prop + 
  out_ns$mean$b0.4*response$flow_rate +
  out_ns$mean$b0.5*response$mean_sal + 
  out_ns$mean$b0.9*response$mussel_size

# Submergence
sub.hat <- out_ns$mean$b3.0 + 
  out_ns$mean$b3.12*response$ele + 
  out_ns$mean$b3.13*response$head

# Flow
flow.hat <- out_ns$mean$b4.0 + 
  out_ns$mean$b4.11*response$sa_height

# SA height
sah.hat <- out_ns$mean$b11.0 + 
  out_ns$mean$b11.12*response$ele + 
  out_ns$mean$b11.14*response$main

# Salinity
sal.hat <- out_ns$mean$b5.0 + 
  out_ns$mean$b5.2*response$mean_water_temp+ out_ns$mean$b5.2q*response$mean_water_temp^2+
  out_ns$mean$b5.14*response$main

### RESIDUALS ################################################################## 
temp.res <- response$mean_water_temp		  
size.res <- response$mussel_size
ele.res <- response$ele
head.res <- response$head
main.res <- response$main
growth.res <- response$growth_rate-growth.hat 	
sub.res <- response$sub_prop-sub.hat
flow.res <- response$flow_rate -flow.hat
sah.res <- response$sa_height-sah.hat
sal.res <- response$mean_sal-sal.hat

# Residual matrix
resid.mat <- data.frame(temp.res,flow.res,sal.res,size.res,growth.res,sub.res,ele.res,sah.res,head.res,main.res)


### MODEL R2 ################################################################### 
summary(lm(response$growth_rate ~ growth.hat))
summary(lm(response$sub_prop~sub.hat))
summary(lm(response$mean_sal~sal.hat))
summary(lm(response$flow_rate~flow.hat))
summary(lm(response$sa_height~sah.hat))
summary(lm(response$mean_sal~sal.hat))

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

# All combos significant
# "temp.res-sal.res"    "temp.res-size.res"   "flow.res-size.res"   
# "sal.res-size.res"    "flow.res-sub.res"    "flow.res-ele.res"   
# "ele.res-head.res"    "flow.res-main.res"   "growth.res-main.res" 
# "head.res-main.res"   "temp.res-flow.res"   "temp.res-main.res"  
# "ele.res-main.res" 

# Ecologically meaningful relationships
# "flow.res-sub.res"    "flow.res-ele.res"   
# "flow.res-main.res"   "growth.res-main.res" 
  


# Growth ~
#     main  # not significant
#  flow~
#     sub
#     main
#     ele


################################################################################
### SEM MODEL WITH ADDED PATHWAYS  ###
################################################################################

sink("growth.added.paths.model.txt")		
cat("
model {

### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Growth	### 
	  growth[i] ~ dnorm(growth.hat[i],tau.growth)                            
		growth.hat[i] <- b0.0 + 
	  	b0.2*temp[i] +
	  	b0.3*sub_prop[i] + 
	  	b0.4*flow_rate[i] + 
	  	b0.5*mean_sal[i] + 
	  	b0.9*mussel_size[i]

## Experiment level variables ##	  	
#2. Temp
  temp[i] ~ dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0  

#5. Salinity
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]

#11. Sa_height
 sa_height[i] ~dnorm(height.hat[i],tau.height)
  height.hat[i] <- b11.0 +   
    b11.12*ele[i]+
    b11.14*main[i]

	}

## Plot level variables ##
for(j in 1:30) {
#3. Submergence Time
  sub_prop[j] ~ dnorm(sub.hat[j], tau.sub)
  sub.hat[j] <- b3.0 + 
   b3.12*ele[j] + 
   b3.13*head[j]

#4. Flow
  flow_rate[j] ~ dnorm(flow.hat[j], tau.flow)
    flow.hat[j] <- b4.0 + 
      b4.3*sub_prop[j] + # ADDED
      b4.11*sa_height_m[j] +
      b4.12*ele[j] +  # ADDED
      b4.14*main[j]  #ADDED
	}
	
### PRECISION VARIABLES  #######################################################
	tau.growth <- 1/(sigma.growth*sigma.growth)
	tau.temp <- 1/(sigma.temp*sigma.temp)
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.flow <- 1/(sigma.flow*sigma.flow)
	tau.sal <- 1/(sigma.sal*sigma.sal)
	tau.height <- 1/(sigma.height*sigma.height)

### PRIORS  ####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001); b0.3 ~ dnorm(0,0.00001); 
	b0.4 ~ dnorm(0,0.00001); b0.5 ~ dnorm(0,0.00001); b0.9 ~ dnorm(0,0.00001); 
	sigma.growth ~ dunif(0,100)

	b2.0 ~ dnorm(0,0.00001); sigma.temp ~ dunif(0,100)
	
  b3.0 ~ dnorm(0,0.00001); b3.12 ~ dnorm(0,0.00001); b3.13 ~ dnorm(0,0.00001) 
  sigma.sub ~ dunif(0,100)
  
  b4.0 ~ dnorm(0,0.00001); b4.3 ~ dnorm(0,0.00001); b4.11 ~ dnorm(0,0.00001) 
	b4.12 ~ dnorm(0,0.00001); b4.14 ~ dnorm(0,0.00001); sigma.flow ~ dunif(0,100);

  b5.0 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);b5.2q ~ dnorm(0,0.00001)
  b5.14 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)
  
  b11.0 ~ dnorm(0,0.00001);  b11.12 ~ dnorm(0,0.00001)
  b11.14 ~ dnorm(0,0.00001); sigma.height ~ dunif(0,100);

}
    ",fill=TRUE)
sink()
# end of JAGS code creation

### CREATING OBJECTS TO HAND TO JAGS  ###########################################
data = list(N = N,
            growth = response$growth_rate,
            temp = response$mean_water_temp, 
            sub_prop = response$sub_prop, 
            flow_rate = response$flow_rate, 
            mean_sal = response$mean_sal,
            mussel_size = response$mussel_size, 
            sa_height = response$sa_height,
            sa_height_m = response$sa_height_m,
            ele = response$ele, 
            head = response$head, 
            main = response$main)   

parameters <- c("b0.0","b0.2","b0.3","b0.4","b0.5","b0.9","b0.14","sigma.growth",
                
                "b2.0","b2.4","b2.14","sigma.temp",
                
                "b3.0", "b3.12","b3.13","sigma.sub",
                
                "b4.0","b4.3","b4.11","b4.12","b4.14","sigma.flow",
                
                "b5.0","b5.2","b5.2q","b5.14","sigma.sal",
                
                "b11.0","b11.10","b11.12","b11.14","sigma.height")

inits <- function(){list(b0.0=0,b0.2=0,b0.3=0,b0.4=0,b0.5=0,b0.9=0,sigma.growth=100, 
                         
                         b2.0=0,sigma.temp=100,
                         
                         b3.0=0, b3.12=0, b3.13=0, sigma.sub=100,
                         
                         b4.0=0,b4.3=0,b4.11=0,b4.12=0,b4.14=0,sigma.flow=100,
                         
                         b5.0=0,b5.2=0,b5.2q=0,b5.14=0,sigma.sal=100,
                   
                         b11.0=0,b11.12=0,b11.14=0,sigma.height=100)}

# MCMC settings to hand to JAGS
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain     
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

## RUN MODEL IN JAGS ###########################################################
added_paths_out <- jags(data=data, 
                             inits=inits, 
                             parameters.to.save=parameters, 
                             model.file="growth.added.paths.model.txt", 
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

sink("growth.final.model.txt")		
cat("
model {
### LIKELIHOODS  ###############################################################
	for(i in 1:N) {
	
## Growth ##	 
	  growth[i] ~ dnorm(growth.hat[i],tau.growth)                            
		growth.hat[i] <- b0.0 + 
	  	b0.2*temp[i] +
	  	b0.3*sub_prop[i] + 
	  	b0.4*flow_rate[i] + 
	  	b0.5*mean_sal[i] + 
	  	b0.9*mussel_size[i] 

## Experiment level varaibles ##	  	
#2. Temp
  temp[i] ~ dnorm(temp.hat[i], tau.temp)
  temp.hat[i] <- b2.0  

#5. Salinity
  mean_sal[i] ~dnorm(sal.hat[i],tau.sal)
  sal.hat[i] <-b5.0 +
    b5.2*temp[i] + b5.2q*temp[i]*temp[i] +
    b5.14*main[i]
	}

## Plot level variables  ##
for(j in 1:30) {
#3. SUbmergence Time
  sub_prop[j] ~ dnorm(sub.hat[j], tau.sub)
  sub.hat[j] <- b3.0 + 
   b3.12*ele[j] + 
   b3.13*head[j]

#4. Flow
  flow_rate[j] ~ dnorm(flow.hat[j], tau.flow)
    flow.hat[j] <- b4.0 + 
      b4.3*sub_prop[j] + # ADDED
      b4.13*head[j] + # ADDED
      b4.14*main[j] + # ADDED
      b4.12*ele[j]   # ADDED
	}
	
### PRECISION VARIABLES  #######################################################
	tau.growth <- 1/(sigma.growth*sigma.growth)
	tau.temp <- 1/(sigma.temp*sigma.temp)
	tau.sub <- 1/(sigma.sub*sigma.sub)
	tau.flow <- 1/(sigma.flow*sigma.flow)
	tau.sal <- 1/(sigma.sal*sigma.sal)

### PRIORS  ####################################################################
	b0.0 ~ dnorm(0,0.00001); b0.1 ~ dnorm(0,0.00001); b0.2 ~ dnorm(0,0.00001)
	b0.3 ~ dnorm(0,0.00001); b0.4 ~ dnorm(0,0.00001); b0.5 ~ dnorm(0,0.00001)
	b0.9 ~ dnorm(0,0.00001);sigma.growth ~ dunif(0,100)
	
	b2.0 ~ dnorm(0,0.00001); sigma.temp ~ dunif(0,100)
	
  b3.0 ~ dnorm(0,0.00001); b3.12 ~ dnorm(0,0.00001); b3.13 ~ dnorm(0,0.00001); sigma.sub ~ dunif(0,100)
  
  b4.0 ~ dnorm(0,0.00001); b4.3 ~ dnorm(0,0.00001); b4.12 ~ dnorm(0,0.00001)
  b4.13 ~ dnorm(0,0.00001);b4.14 ~ dnorm(0,0.00001); sigma.flow ~ dunif(0,100)
  
  b5.0 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);b5.2q ~ dnorm(0,0.00001)
  b5.14 ~ dnorm(0,0.00001); sigma.sal ~ dunif(0,100)

}
    ",fill=TRUE)
sink()
# end of JAGS code creation

### CREATING OBJECTS TO HAND TO JAGS  ###########################################

N=90  # number of seasons x plots
data = list(N = N,
            growth = response$growth_rate,
            temp = response$mean_water_temp, 
            sub_prop = response$sub_prop, 
            flow_rate = response$flow_rate, 
            mean_sal = response$mean_sal,
            mussel_size = response$mussel_size, 
            ele = response$ele, 
            head = response$head, 
            main = response$main)   

parameters <- c("b0.0","b0.2","b0.3","b0.4","b0.5","b0.9","sigma.growth",
                
                "b2.0","sigma.temp",
                
                "b3.0", "b3.12","b3.13","sigma.sub",
                
                "b4.0","b4.3","b4.12","b4.13","b4.14","sigma.flow",
                
                "b5.0","b5.2","b5.2q","b5.14","sigma.sal")

inits <- function(){list(b0.0=0,b0.2=0,b0.3=0,b0.4=0,b0.5=0,b0.9=0,sigma.growth=100, 
                         
                         b2.0=0, sigma.temp=100,
                         
                         b3.0=0, b3.12=0, b3.13=0, sigma.sub=100,
                         
                         b4.0=0,b4.3=0,b4.12=0,b4.13=0,b4.14=0,sigma.flow=100,
                         
                         b5.0=0,b5.2=0,b5.2q=0,b5.14=0,sigma.sal=100)}

# MCMC settings to hand to jags
nc <- 3      # number of chains
ni <- 50000  # number of samples for each chain   
nb <- 5000   # number of samples to discard for burn in
nt <- 5      # thinning rate

## RUN MODEL IN JAGS ###########################################################
final.out <- jags(data=data, 
                  inits=inits, 
                  parameters.to.save=parameters, 
                  model.file="growth.final.model.txt", 
                  n.chains=nc, 
                  n.iter=ni, 
                  n.burnin=nb, 
                  n.thin=nt, 
                  DIC=TRUE, 
                  parallel = T)

# Print some basic results
print(final.out,digits=4) 

## Export results for prediction maps ##
# saveRDS(final.out,"model_out/growth_SEM_final_output.rds")

# ## Export results for table ##
# hyp_out_final <- as.data.frame(final.out$summary)  # result summary
# hyp_out_final <- round(hyp_out_final,digits = 3)  # round results
# hyp_out_final$median_95 <- paste(hyp_out_final$`50%`," (",hyp_out_final$`2.5%`,",",hyp_out_final$`97.5%`,")",sep ="")  # create column for median and 95CI
# hyp_out_final <- hyp_out_final %>% select(median_95,Rhat,overlap0)
# write.csv(hyp_out_final,"~/Documents/School/Master's/Thesis documents/ch2 pub/results/model outputs/growth_hyp_final_out_water_temp.csv")

### PREDICTION EQUATIONS #######################################################

# Growth
growth.hat <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*response$mean_water_temp +
  final.out$mean$b0.3*response$sub_prop + 
  final.out$mean$b0.4*response$flow_rate +
  final.out$mean$b0.5*response$mean_sal + 
  final.out$mean$b0.9*response$mussel_size 

# Submergance
sub.hat <- final.out$mean$b3.0 + 
  final.out$mean$b3.12*response$ele + 
  final.out$mean$b3.13*response$head

# Flow
flow.hat <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*response$sub_prop + 
  final.out$mean$b4.12*response$ele + 
  final.out$mean$b4.13*response$head+
  final.out$mean$b4.14*response$main

# Salinity
sal.hat <- final.out$mean$b5.0 + 
  final.out$mean$b5.2*response$mean_water_temp + 
  final.out$mean$b5.2q*response$mean_water_temp^2 +
  final.out$mean$b5.14*response$main

### RESIDUALS ###################################################################
temp.res <- response$mean_water_temp		 
size.res <- response$mussel_size
ele.res <- response$ele
head.res <- response$head
main.res <- response$main
growth.res <- response$growth_rate-growth.hat 	 
sub.res <- response$sub_prop-sub.hat
flow.res <- response$flow_rate -flow.hat
sal.res <- response$mean_sal-sal.hat

# Residual matrix
final.res <- data.frame(temp.res,flow.res,sal.res,size.res,growth.res,sub.res,ele.res,head.res,main.res)

# Export growth residuals
resid_out <- data.frame(
  plot = response$plot,
  season = response$season,
  resid = growth.res
)
saveRDS(resid_out,"model_out/growth_resid.rds")

###  FINAL MODEL R2 ############################################################
summary(lm(response$growth_rate ~ growth.hat))
summary(lm(response$sub_prop~sub.hat))
summary(lm(response$mean_sal~sal.hat))
summary(lm(response$flow_rate~flow.hat))

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

### FOR GROWTH #################################################################

### TEMP ###

# min
growth.predict.2 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*min(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)	       

# max
growth.predict.3 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*max(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)      	  

# QUERY
(growth.predict.3-growth.predict.2)/(max(response$growth_rate)-min(response$growth_rate)) 

### SUB ###

# min
growth.predict.4 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*min(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)	       	

# max
growth.predict.5 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*max(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)   	

# QUERY
(growth.predict.5-growth.predict.4)/(max(response$growth_rate)-min(response$growth_rate))

### FLOW ###

# min
growth.predict.6 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*min(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)	       

# max
growth.predict.7 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*max(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)    	  

# QUERY
(growth.predict.7-growth.predict.6)/(max(response$growth_rate)-min(response$growth_rate))

### SAL ###

# min
growth.predict.8 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*min(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)	       	

# max
growth.predict.9 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*max(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)     	  

# QUERY
(growth.predict.9-growth.predict.8)/(max(response$growth_rate)-min(response$growth_rate))

### SIZE ###

# min
growth.predict.10 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*min(response$mussel_size)   	

# max
growth.predict.11 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*max(response$mussel_size)     	  

# QUERY
(growth.predict.11-growth.predict.10)/(max(response$growth_rate)-min(response$growth_rate))

### For SUB ####################################################################

### ELE ###

# min
sub.predict.2 <- final.out$mean$b3.0 + 
  final.out$mean$b3.12*min(response$ele) + 
  final.out$mean$b3.13*mean(response$head) 	       	

# max
sub.predict.3 <- final.out$mean$b3.0 + 
  final.out$mean$b3.12*max(response$ele) + 
  final.out$mean$b3.13*mean(response$head)    	 

# QUERY
(sub.predict.3-sub.predict.2)/(max(response$sub_prop)-min(response$sub_prop))

### HEAD ###

# min
sub.predict.4 <- final.out$mean$b3.0 + 
  final.out$mean$b3.12*mean(response$ele) + 
  final.out$mean$b3.13*min(response$head) 		       

# max
sub.predict.5 <- final.out$mean$b3.0 + 
  final.out$mean$b3.12*mean(response$ele) + 
  final.out$mean$b3.13*max(response$head)       	 

# QUERY
(sub.predict.5-sub.predict.4)/(max(response$sub_prop)-min(response$sub_prop))

### FOR FLOW ###################################################################
  
### SUB ###

# min
flow.predict.2 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*min(response$sub_prop) + 
  final.out$mean$b4.12*mean(response$ele) + 
  final.out$mean$b4.13*mean(response$head)+
  final.out$mean$b4.14*mean(response$main)	       	

# max
flow.predict.3 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*max(response$sub_prop) + 
  final.out$mean$b4.12*mean(response$ele) + 
  final.out$mean$b4.13*mean(response$head)+
  final.out$mean$b4.14*mean(response$main)    	 

# QUERY
(flow.predict.3-flow.predict.2)/(max(response$flow_rate)-min(response$flow_rate))

### ELE ###

# min
flow.predict.4 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*mean(response$sub_prop) + 
  final.out$mean$b4.12*min(response$ele) + 
  final.out$mean$b4.13*mean(response$head)+
  final.out$mean$b4.14*mean(response$main)      

# max
flow.predict.5 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*mean(response$sub_prop) + 
  final.out$mean$b4.12*max(response$ele) + 
  final.out$mean$b4.13*mean(response$head)+
  final.out$mean$b4.14*mean(response$main)     	  

# QUERY
(flow.predict.5-flow.predict.4)/(max(response$flow_rate)-min(response$flow_rate))

### HEAD ###

# min
flow.predict.6 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*mean(response$sub_prop) + 
  final.out$mean$b4.12*mean(response$ele) + 
  final.out$mean$b4.13*min(response$head)+
  final.out$mean$b4.14*mean(response$main)     	

# max
flow.predict.7 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*mean(response$sub_prop) + 
  final.out$mean$b4.12*mean(response$ele) + 
  final.out$mean$b4.13*max(response$head)+
  final.out$mean$b4.14*mean(response$main)     	  

# QUERY
(flow.predict.7-flow.predict.6)/(max(response$flow_rate)-min(response$flow_rate))

### MAIN ###

# min
flow.predict.8 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*mean(response$sub_prop) + 
  final.out$mean$b4.12*mean(response$ele) + 
  final.out$mean$b4.13*mean(response$head)+
  final.out$mean$b4.14*min(response$main)       	

# max
flow.predict.9 <- final.out$mean$b4.0 + 
  final.out$mean$b4.3*mean(response$sub_prop) + 
  final.out$mean$b4.12*mean(response$ele) + 
  final.out$mean$b4.13*mean(response$head)+
  final.out$mean$b4.14*max(response$main)      	 

# QUERY
(flow.predict.9-flow.predict.8)/(max(response$flow_rate)-min(response$flow_rate))

### FOR SAL ####################################################################

# Prediction Equations
sal.predict.1 <- sal.hat <- final.out$mean$b5.0 + 
  final.out$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
  final.out$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
  final.out$mean$b5.14*mean(response$main)    

### TEMP ###

# min
sal.predict.2 <- final.out$mean$b5.0 + 
  final.out$mean$b5.2*min(response$mean_water_temp,na.rm = T) + 
  final.out$mean$b5.2q*min(response$mean_water_temp,na.rm = T)^2+
  final.out$mean$b5.14*mean(response$main)	       

# max
sal.predict.3 <- final.out$mean$b5.0 + 
  final.out$mean$b5.2*max(response$mean_water_temp,na.rm = T) + 
  final.out$mean$b5.2q*max(response$mean_water_temp,na.rm = T)^2+
  final.out$mean$b5.14*mean(response$main)     	 

# QUERY
(sal.predict.3-sal.predict.2)/(max(response$mean_sal)-min(response$mean_sal))

### MAIN ###

# min
sal.predict.6 <- final.out$mean$b5.0 + 
  final.out$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
  final.out$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
  final.out$mean$b5.14*min(response$main)      	

# max
sal.predict.7 <- final.out$mean$b5.0 + 
  final.out$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
  final.out$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
  final.out$mean$b5.14*max(response$main) 

# QUERY
(sal.predict.7-sal.predict.6)/(max(response$mean_sal)-min(response$mean_sal))

################################################################################
## TOTAL EFFECTS QUERIES ##
################################################################################

### TEMP ###

# min
growth.predict.tot2 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*min(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*(final.out$mean$b5.0 + 
                         final.out$mean$b5.2*min(response$mean_water_temp,na.rm = T) + 
                         final.out$mean$b5.2q*min(response$mean_water_temp,na.rm = T)^2+
                         final.out$mean$b5.14*mean(response$main)) + 
  final.out$mean$b0.9*mean(response$mussel_size) 

# max
growth.predict.tot3 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*max(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*mean(response$flow_rate) +
  final.out$mean$b0.5*(final.out$mean$b5.0 + 
                         final.out$mean$b5.2*max(response$mean_water_temp,na.rm = T) + 
                         final.out$mean$b5.2q*max(response$mean_water_temp,na.rm = T)^2+
                         final.out$mean$b5.14*mean(response$main)) + 
  final.out$mean$b0.9*mean(response$mussel_size) 

# QUERY
(growth.predict.tot3-growth.predict.tot2)/(max(response$growth_rate)-min(response$growth_rate))

## ELE ##

# min
growth.predict.tot2 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*(final.out$mean$b3.0 + 
                         final.out$mean$b3.12*min(response$ele) + 
                         final.out$mean$b3.13*mean(response$head)) + 
  final.out$mean$b0.4*(final.out$mean$b4.0 + 
                         final.out$mean$b4.3*(final.out$mean$b3.0 + 
                                                final.out$mean$b3.12*min(response$ele) + 
                                                final.out$mean$b3.13*mean(response$head)) + 
                         final.out$mean$b4.12*min(response$ele) + 
                         final.out$mean$b4.13*mean(response$head)+
                         final.out$mean$b4.14*mean(response$main)) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size) 

# max
growth.predict.tot3 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*(final.out$mean$b3.0 + 
                         final.out$mean$b3.12*max(response$ele) + 
                         final.out$mean$b3.13*mean(response$head)) + 
  final.out$mean$b0.4*(final.out$mean$b4.0 + 
                         final.out$mean$b4.3*(final.out$mean$b3.0 + 
                                                final.out$mean$b3.12*max(response$ele) + 
                                                final.out$mean$b3.13*mean(response$head)) + 
                         final.out$mean$b4.12*max(response$ele) + 
                         final.out$mean$b4.13*mean(response$head)+
                         final.out$mean$b4.14*mean(response$main)) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size) 

# QUERY
(growth.predict.tot3-growth.predict.tot2)/(max(response$growth_rate)-min(response$growth_rate))

### HEAD ###

# min
growth.predict.tot2 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*(final.out$mean$b3.0 + 
                         final.out$mean$b3.12*mean(response$ele) + 
                         final.out$mean$b3.13*min(response$head)) + 
  final.out$mean$b0.4*(final.out$mean$b4.0 + 
                         final.out$mean$b4.3*(final.out$mean$b3.0 + 
                                                final.out$mean$b3.12*mean(response$ele) + 
                                                final.out$mean$b3.13*min(response$head)) + 
                         final.out$mean$b4.12*mean(response$ele) + 
                         final.out$mean$b4.13*min(response$head)+
                         final.out$mean$b4.14*mean(response$main)) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size)

# max
growth.predict.tot3 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*(final.out$mean$b3.0 + 
                         final.out$mean$b3.12*mean(response$ele) + 
                         final.out$mean$b3.13*max(response$head)) + 
  final.out$mean$b0.4*(final.out$mean$b4.0 + 
                         final.out$mean$b4.3*(final.out$mean$b3.0 + 
                                                final.out$mean$b3.12*mean(response$ele) + 
                                                final.out$mean$b3.13*max(response$head)) + 
                         final.out$mean$b4.12*mean(response$ele) + 
                         final.out$mean$b4.13*max(response$head)+
                         final.out$mean$b4.14*mean(response$main)) +
  final.out$mean$b0.5*mean(response$mean_sal) + 
  final.out$mean$b0.9*mean(response$mussel_size) 

# QUERY
(growth.predict.tot3-growth.predict.tot2)/(max(response$growth_rate)-min(response$growth_rate))

### MAIN ###

# min
growth.predict.tot2 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*min(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*(final.out$mean$b4.0 + 
                         final.out$mean$b4.3*mean(response$sub_prop) + 
                         final.out$mean$b4.12*mean(response$ele) + 
                         final.out$mean$b4.13*mean(response$head)+
                         final.out$mean$b4.14*min(response$main)) +
  final.out$mean$b0.5*(final.out$mean$b5.0 + 
                         final.out$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
                         final.out$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
                         final.out$mean$b5.14*min(response$main)) + 
  final.out$mean$b0.9*mean(response$mussel_size) 

# max
growth.predict.tot3 <- final.out$mean$b0.0 + 
  final.out$mean$b0.2*mean(response$mean_water_temp,na.rm = T) +
  final.out$mean$b0.3*mean(response$sub_prop) + 
  final.out$mean$b0.4*(final.out$mean$b4.0 + 
                         final.out$mean$b4.3*mean(response$sub_prop) + 
                         final.out$mean$b4.12*mean(response$ele) + 
                         final.out$mean$b4.13*mean(response$head)+
                         final.out$mean$b4.14*max(response$main)) +
  final.out$mean$b0.5*(final.out$mean$b5.0 + 
                         final.out$mean$b5.2*mean(response$mean_water_temp,na.rm = T) + 
                         final.out$mean$b5.2q*mean(response$mean_water_temp,na.rm = T)^2+
                         final.out$mean$b5.14*max(response$main)) + 
  final.out$mean$b0.9*mean(response$mussel_size) 

# QUERY
(growth.predict.tot3-growth.predict.tot2)/(max(response$growth_rate)-min(response$growth_rate))

