#Test run  with inprod in a model 

#mod1 works so inprod working now
#mod2 - couple issues to fix
# indexing of lv.coefs, the lv.coef matrix is the other way around to 
# boral because of hte way inprod works, so need to change the indexing in the loops
# get compile error, 
#Doesn't seem to be due to truncated normal or zero
# Is it beause of using lv.coefs twice, once fo rintercept and again for coefs? 

#guide: 
#https://r-nimble.org/quick-guide-for-converting-from-jags-or-bugs-to-nimble
# CJ Brown
#2024-10-12

# Load the nimble package
library(nimble)
library(dplyr)
library(ggplot2)
library(tidyr)

#make up some data
n <- 100
p <- 3
dat <- data.frame(x1 = rnorm(n), 
                  x2 = rnorm(n))
alpha <- c(5,-1, 0)
beta <- matrix(c(-1, 1, 2, -2, 0, -3), 
               nrow = 2)

nx <- 2
X <- model.matrix(~0+ x1 + x2, data = dat)
Y <- lambda <- matrix(NA, nrow = n, ncol = p)

for (ispp in 1:p){
  lambda[,ispp] <- exp(alpha[ispp] + beta[1:nx,ispp] %*% t(X))
  Y[,ispp] <- rpois(n, lambda[,ispp])
}

plot(X[,1], Y[,1])
plot(Y[,1], Y[,3])

mod1 <- nimbleCode({
  ## Data Level ## 
  for(i in 1:n) {
    for(j in 1:p) { 
      log(eta[i,j]) <- beta0[j] + inprod(beta[1:nx, j], x[i, 1:nx])
      y[i,j] ~ dpois(eta[i,j])
    }
  }
  
  for(j in 1:p) { 
    beta0[j] ~ dnorm(0, sd = 10) 
    
    for (k in 1:nx) {
      beta[k,j] ~ dnorm(0, sd = 10)
    }
  } 
  
}
)


## constants, data, and initial values
constants <- list(n = n, p = p, nx = 2)
# inits <- list(beta0 = 0, beta1 = 0)
colnames(Y) <- NULL
dat2 <- list(y = Y, x = X)
inits <- list(beta0 = alpha, beta1 = beta)

bmod <- nimbleModel(code = mod1, constants = constants, data = dat2,
                    inits = inits,
                    check = FALSE)


bspec <- configureMCMC(bmod)
# bspec$addSampler(type = 'slice', target = 'X.coefs')
# bspec$addSampler(type = 'slice', target = 'lv.coefs')
build_bmod <- buildMCMC(bspec)
cMod <- compileNimble(bmod)
cMCMC <- compileNimble(build_bmod, project = bmod)

samples <- runMCMC(cMCMC, niter = 5000, nburnin = 1000,
                   thin = 5, nchains = 1)
samples <- as.data.frame(samples)

probs <- c(0.025, 0.5, 0.975)
names(samples)
quantile(samples[,1], probs)
quantile(samples[,2], probs)
quantile(samples[,3], probs)
quantile(samples[,4], probs)
quantile(samples[,5], probs)
quantile(samples[,6], probs)
head(samples)

#
# Testrun LV model 
#

mod2 <- nimbleCode({
  ## Data Level ## 
  for(i in 1:n) {
    for(j in 1:p) { 
      log(eta[i,j]) <- lv.coef[1,j] + inprod(lv.coef[2:num.lv, j], lvs[i, 1:num.lv])
      y[i,j] ~ dpois(eta[i,j])
    }
  }
  
  ## Separate species intercepts
  for(j in 1:p) { 
    lv.coef[1, j] ~ dnorm(0, sd = 10) 
  }
  
  ## Constraints to 0 on upper diagonal
  # for(i in 1:(num.lv-1)) {
    # for(j in (i+2):(num.lv+1)) {
      lv.coefs[num.lv+1,1] <- 0
    # } 
  # }
  
  for(i in 1:num.lv) {
    lv.coefs[i+1, i] ~ dnorm(0, sd = 10) 
  } 
  
  ## Sign constraints on diagonal elements
  # for(i in 2:num.lv) { 
    # for(j in 2:i) { 
      # lv.coefs[i,j] ~ dnorm(0,10) #this was: dnorm(0,0.1)I(0,) in JAGS
      lv.coefs[2,2] ~ T(dnorm(0, 0.1), 0, a) #should truncate the normal
    # } 
  # } 
  
  ## Free lower diagonals
  # for(i in (num.lv+1):p) {
    for(j in 2:(num.lv+1)) {
      lv.coefs[j,3] ~ dnorm(0,0.1)
    } 
  #} 
  
}
)


## constants, data, and initial values
constants <- list(n = n, p = 3, num.lv = 2)
# inits <- list(beta0 = 0, beta1 = 0)
colnames(Y) <- NULL
dat2 <- list(y = Y)
inits <- list(beta0 = alpha, beta1 = beta)

bmod <- nimbleModel(code = mod2, constants = constants, data = dat2,
                    # inits = inits,
                    check = FALSE)


bspec <- configureMCMC(bmod)
# bspec$addSampler(type = 'slice', target = 'X.coefs')
# bspec$addSampler(type = 'slice', target = 'lv.coefs')
build_bmod <- buildMCMC(bspec)
cMod <- compileNimble(bmod)
cMCMC <- compileNimble(build_bmod, project = bmod)

samples <- runMCMC(cMCMC, niter = 5000, nburnin = 1000,
                   thin = 5, nchains = 1)
samples <- as.data.frame(samples)

probs <- c(0.025, 0.5, 0.975)
names(samples)
quantile(samples[,1], probs)
quantile(samples[,2], probs)
quantile(samples[,3], probs)
quantile(samples[,4], probs)
quantile(samples[,5], probs)
quantile(samples[,6], probs)
head(samples)
