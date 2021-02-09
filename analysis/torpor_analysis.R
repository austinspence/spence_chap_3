##################################
## hummingbird organ analysis
## created on: Feb. 8th, 2021
## created by: Austin S.

# libraries
library(corrplot)
library(R2jags)
library(MCMCvis)
library(ggplot2)
library(dplyr)


# clean things up
rm(list = ls())

# good functions
stdize <- function(x){
  (x - mean(x))/(sd(x))
}

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


## read in the data
torpor <- read.csv("./data/working/working_torpor_data.csv")


# standardize elevation
torpor$std.elev <- stdize(torpor$elev)


## torpor use analysis
torpor.use.func.1 <- function(){
  
  #### Priors
  b.0 ~ dnorm(0, 0.001)
  b.temp ~ dnorm(0, 0.001)    # lowest ambient temperature
  b.juli ~ dnorm(0, 0.001)    # julian date of test
  b.elev ~ dnorm(0, 0.001)    # elevation of origin
  b.site ~ dnorm(0, 0.001)    # barcroft vs. bigpine
  b.capt ~ dnorm(0, 0.001)    # days in captivity
  
  tau ~ dgamma(0.001, 0.001)
  tau.ind ~ dgamma(0.001, 0.001)
  
  
  ## Process model
  for(i in 1:n.obs){
    torpor[i] ~ dbern(p[i])
    
    logit(p[i]) <- mu[i]
    
    mu[i] <- b.0 +
      b.temp * temp[i] +
      b.juli * juli[i] +
      b.elev * elev[i] +
      b.site * site[i] +
      b.capt * capt[i] +
      b.ind[ind[i]]
    
  }
  
  ## Random intercept to account for year 
  for(j in 1:n.ind){
    b.ind[j] ~ dnorm(0, tau.ind)
  }
  
  
}



## the data
torpor.use.list.1 <- list(
  torpor = torpor$torpor_use,
  temp = torpor$lowest_temp,
  juli = torpor$julian_testing_date,
  elev = torpor$std.elev,
  site = torpor$location,
  capt = torpor$days_in_captivity,
  ind = torpor$individual,
  n.obs = length(torpor$capture_id),
  n.ind = length(unique(torpor$individual))
)


## run the model
ptm <- proc.time()

torpor.use.mod.1 <- jags(data = torpor.use.list.1,
                     parameters.to.save = c("b.0",
                                            "b.temp",
                                            "b.elev",
                                            "b.juli",
                                            "b.site",
                                            "b.capt"),
                     model.file = torpor.use.func.1,
                     n.chains = 3, 
                     n.iter = 10000,
                     n.burnin = 1000,
                     n.thin = 10)

proc.time() -ptm

MCMCsummary(torpor.use.mod.1)
MCMCplot(torpor.use.mod.1, excl = "deviance", ref_ovl = TRUE)














