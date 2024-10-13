#Convert jags code from BORAL to NIMBLE code

# CJ Brown
#2024-10-12

# Load the nimble package
library(nimble)
library(dplyr)
library(ggplot2)
library(tidyr)

url <- "https://raw.githubusercontent.com/cbrown5/BenthicLatent/refs/heads/master/data-raw/JuvUVCSites_with_ReefTypes_16Jun2016.csv"
dat <- read.csv(url)

#
# Define the model 
#

mod1 <- nimbleCode({
  beta0 ~ dnorm(0, sd = 10000)
  beta1 ~ dnorm(0, sd = 10000)
  for (i in 1:N) {
    log(lambda[i]) <- beta0 + beta1 * x[i]
    topa[i] ~ dpois(lambda[i])
    }
  }
)

## constants, data, and initial values
constants <- list(N = nrow(dat))
dat$x <- log(dat$mindist)

dat2 <- select(dat, topa = pres.topa, x)

inits <- list(beta0 = 0, beta1 = 0)

## create the model object
glmModel <- nimbleModel(code = mod1, constants = constants, data = dat2, 
                         inits = inits, check = FALSE)
spec <- configureMCMC(glmModel, nodes = NULL)
spec$addSampler(type = 'slice', target = 'beta0')
spec$addSampler(type = 'slice', target = 'beta1')
customGlmMCMC <- buildMCMC(spec)
CglmModel <- compileNimble(glmModel)
CglmMCMC <- compileNimble(customGlmMCMC, project = glmModel)

# Sample from the model
samples <- runMCMC(CglmMCMC, niter = 5000, nburnin = 1000,
                   thin = 5, nchains = 3)
samples <- as.data.frame(samples)
# Put chains into one column
samples2 <- samples %>%
  mutate(iteration = 1:nrow(samples)) %>%
  pivot_longer(cols = -iteration, 
               names_to = "parameter", 
               values_to = "value" 
               ) %>%
  separate(parameter, into = c("chain", "parameter"), sep = "\\.") %>%
  pivot_wider(names_from = parameter, values_from = value) 


## plot the results
samples2 %>%
  ggplot() + 
  aes(x = beta0) + 
  geom_density()

samples2 %>%
  ggplot() + 
  aes(x = beta1) + 
  geom_density()
