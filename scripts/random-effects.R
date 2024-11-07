# Explanation of random effects structure for GLMs
# CJ Brown
# 2024-11-07

library(tidyverse)

nsites <- 30
nquads <- 10
n <- nsites * nquads

#intercept
alpha <- 1 

#simulate site z effects
sigma2_site <- 1 #var
site_z <- rnorm(nsites, 0, sd = sqrt(sigma2_site))

beta <- -2
site_x <- rnorm(nsites, 0, 1)

# beta <- 1
# x <- rnorm(nsites, 0, 1)

# Simulate data
#create data.frame with site and quad IDs
dat <- expand.grid(quad = 1:nquads, site = 1:nsites)
dat$quad2 <- paste(dat$site, dat$quad, sep = "_")
# View(dat)
#add site z effects
dat$site_z <- rep(site_z, each = nquads)
dat$site_x <- rep(site_x, each = nquads)
# View(dat)

#simulate quad effects
sigma2_quad <- 2^2
dat$quad_errors <- rnorm(n, 0, sd = sqrt(sigma2_quad))
# View(dat)

#simulate response
dat$y <- alpha + beta*dat$site_x + dat$site_z + dat$quad_errors

mean(dat$y)
dat %>% group_by(site) %>% summarise(mean_y = mean(y)-alpha)
site_z

# fit a model with lme4
library(lme4)
m1 <- lmer(y ~ 1 +site_x+ (1|site), data = dat)
summary(m1)

#wrong way that doesn't account for psuedo-replication
m2 <- lm(y ~ 1 + site_x, data = dat)
summary(m2)

#variance components
#take site variance from model and residual variance, divide them
total_var <- 0.0105 + 1.0284
total_var
site_var <- 1.0284
site_var / total_var

# Plot the data
dat %>% ggplot(aes(x = site, y = y)) +
  geom_point(alpha = 0.5) +
  theme_minimal()

#plot response vs site_x
dat %>% ggplot(aes(x = site_x, y = y, color = factor(site))) +
  geom_point(alpha = 0.5) +
  theme_minimal()

# get BLUPs out of the model
plot(ranef(m1)$site[,1], site_z)
abline(0,1)
