#Convert jags code from BORAL to NIMBLE code

#TODO: Compilation Issue with current model. Its probably the use of inprod function, need to debug that

#guide: 
#https://r-nimble.org/quick-guide-for-converting-from-jags-or-bugs-to-nimble
# CJ Brown
#2024-10-12

# Load the nimble package
library(nimble)
library(dplyr)
library(ggplot2)
library(tidyr)

url <- "https://raw.githubusercontent.com/cbrown5/BenthicLatent/refs/heads/master/data-raw/JuvUVCSites_with_ReefTypes_16Jun2016.csv"
dat <- read.csv(url)

#Boral functions
#have sourced them as the package won't open without access to jags
source("scripts/jags-to-nimble/boral-functions.R")

#
# Define the data 
#
#help file here: https://raw.githubusercontent.com/emitanaka/boral/refs/heads/master/man/make.jagsboralmodel.Rd
#response matrix
Y <- select(dat, pres.topa, pres.habili) %>%
  mutate(y3 = rpois(nrow(dat), 10), y4 = rpois(nrow(dat), 2)) %>%
  as.matrix()

n <- nrow(Y)
p <- ncol(Y)
X <- dat %>%
  mutate(logmindist = log(mindist)) %>%
  select(logmindist) %>%
  as.matrix()
num.X <- ncol(X)
num.lv <- 2

# Create the boral txt file
make.jagsboralmodel("poisson", 
                    num.X = num.X, #number covariates (not counting intercept)
                    X.ind = NULL, #indicator for which covariates for which responses (speciess)
                    num.traits = 0, 
                    which.traits = NULL,
                    lv.control = list(num.lv = num.lv, type = "independent"), 
                    row.eff = "none", #ie site random or fixed
                    row.ids = NULL, # rows for multiple row effects
                    offset = NULL, 
                    trial.size = 1, 
                    n, #number of rows in response matrix (ie num samples)
                    p, #number of columsn in response matrix (ie num spp)
                    model.name = "scripts/jags-to-nimble/jags-boral-model.txt",
                    prior.control = list(type = c("normal", "normal", "normal", "uniform"), 
                                         hypparams = c(10, 10, 10, 30), 
                                         ssvs.index = -1, #no svss if -1
                                         ssvs.g = 1e-6, 
                                         ssvs.traitsindex = -1))

#
# Convert to NIMBLE
#

## JAGS model written for boral version 2.0.2 on 2024-10-12 14:39:22.729514 ##

#1. Wrap your model code in nimbleCode({}), directly in R.
#2. Provide information about missing or empty indices
#3. Use nimbleMCMC() as the just-do-it way to run an MCMC, or specify MCMC
# SD instead of precision in dnorm
# I() function, how to implement in NIMBLE? 

mod1 <- nimbleCode({
  ## Data Level ## 
  for(i in 1:n) {
    for(j in 1:p) { 
      #something wrong here: inprod(X.coefs[j,1:nX],X[i,1:nX]) 
      #heck how inprod works
      #exp(lv.coefs[j,1] + )
      # intercept (lv.coefs[j,1] is the intercept for species j) + eta
      y[i,j] ~ dpois(exp(lv.coefs[j,1])) 
      
      y[i] ~ dnorm(beta0 + inprod(beta[1:p], x[i, 1:p]), sd = sigma)
      
    }
  }
  
  
  ## Process level and priors ##
  
  ## Separate species intercepts
  for(j in 1:p) { 
    lv.coefs[j,1] ~ dnorm(0, sd = 10) 
  } 
  
  for(j in 1:p) { 
    X.coefs[j,1] ~ dnorm(0,0.1) 
  } 
}
)


## constants, data, and initial values
constants <- list(n = n, p = 4, nX = num.X)
# inits <- list(beta0 = 0, beta1 = 0)
colnames(Y) <- NULL
dat2 <- list(y = Y, X = X)
## create the model object
inits <- list(X.coefs = matrix(0.1, nrow = p, ncol = num.X)
              )

bmod <- nimbleModel(code = mod1, constants = constants, data = dat2,
                        # inits = inits,
                        check = FALSE)


bspec <- configureMCMC(bmod)
# bspec$addSampler(type = 'slice', target = 'X.coefs')
# bspec$addSampler(type = 'slice', target = 'lv.coefs')
build_bmod <- buildMCMC(bspec)
cMod <- compileNimble(bmod)
cMCMC <- compileNimble(build_bmod, project = bmod)



