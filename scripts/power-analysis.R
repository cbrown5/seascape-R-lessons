# Set seed for reproducibility
set.seed(123)
n <- 20
x <- rnorm(n)
# Simulate y as a function of x with some noise
y <- 0 * x + rnorm(n)
plot(x,y)  
  # Fit model and extract p-value for x coefficient
model <- lm(y ~ x)
summary(model)


# Function to simulate one dataset and return p-value
simulate_regression <- function(n = 20, b = 0) {
  x <- rnorm(n)
  # Simulate y as a function of x with some noise
    y <- b * x + rnorm(n)
  
  # Fit model and extract p-value for x coefficient
  model <- lm(y ~ x)
  summary(model)$coefficients[2, 4]
}

# Run 1000 simulations
n_sims <- 10000
p_values <- replicate(n_sims, simulate_regression())

# Calculate proportion of false positives (alpha = 0.05)
type_I_error_rate <- mean(p_values < 0.05)

# Print results
cat("Type I error rate:", round(type_I_error_rate, 3), "\n")
cat("Expected rate: 0.05\n")

# Optional: create histogram of p-values
hist(p_values, breaks = 30, main = "Distribution of p-values",
     xlab = "p-value", ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2)

#
# Function to simulate student t data and return p-value
#
simulate_t <- function(n = 200, df = 2) {
    x <- rnorm(n)
    # Generate t-distributed residuals with low degrees of freedom
    y <- rt(n, df = df) + rnorm(n)
    
    # Fit linear model (incorrectly assuming normality)
    model <- lm(y ~ x)
    
    summary(model)$coefficients[2, 4]
    #plot(model, 2)
}

# Run 1000 simulations with Poisson data
p_values_t <- replicate(n_sims, simulate_t(df = 1))

# Calculate Type I error rate for Poisson case
type_I_error_rate_simulate_t <- mean(p_values_t < 0.05)

# Print results for both cases
cat("\nNormal residuals Type I error rate:", round(type_I_error_rate, 3), "\n")
cat("student t data Type I error rate:", round(type_I_error_rate_simulate_t, 3), "\n")
cat("Expected rate: 0.05\n")
n = 200
 x <- rnorm(n)
# Generate t-distributed residuals with low degrees of freedom
y <- rt(n, df = 2) + rnorm(n)
    
    # Fit linear model (incorrectly assuming normality)
     model <- lm(y ~ x)
hist(resid(model))
plot(model, 2)

# Plot histograms side by side
par(mfrow = c(1, 2))
hist(p_values, breaks = 30, main = "Normal residuals",
     xlab = "p-value", ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2)

hist(p_values_t, breaks = 30, main = "Student t data",
     xlab = "p-value", ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2)

#
# Function to simulate Poisson data and return p-value
#
simulate_poisson <- function(n = 20, lambda = 1) {
    x <- rnorm(n)
    # Generate Poisson-distributed residuals
    y <- rpois(n, lambda = lambda)
    
    # Fit linear model (incorrectly assuming normality)
    model <- lm(y ~ x)
    
    summary(model)$coefficients[2, 4]
}

# Run 1000 simulations with Poisson data
p_values_poisson <- replicate(n_sims, simulate_poisson(lambda = 1))

# Calculate Type I error rate for Poisson case
type_I_error_rate_simulate_poisson <- mean(p_values_poisson < 0.05)

# Print results for both cases
cat("\nNormal residuals Type I error rate:", round(type_I_error_rate, 3), "\n")
cat("Poisson data Type I error rate:", round(type_I_error_rate_simulate_poisson, 3), "\n")
cat("Expected rate: 0.05\n")

# Plot histograms side by side
par(mfrow = c(1, 2))
hist(p_values, breaks = 30, main = "Normal residuals",
     xlab = "p-value", ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2)

hist(p_values_poisson, breaks = 30, main = "Poisson data",
     xlab = "p-value", ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2)

# ------------------------------- #
# POWER ANALYSIS 
# ------------------------------- #

#Example plot of x vs y
n <- 20
x <- rnorm(n)
# True Poisson regression: log(lambda) = beta0 + beta*x
lambda <- exp(0.5 + 0.2 * x)
y <- rpois(n, lambda)
plot(x, y)


# Power analysis comparing linear vs Poisson regression
simulate_with_effect <- function(n = 20, beta = 0.5) {
    x <- rnorm(n)
    # True Poisson regression: log(lambda) = beta0 + beta*x
    lambda <- exp(0.5 + beta * x)
    y <- rpois(n, lambda)

    # Fit both models
    lm_model <- lm(y ~ x)
    pois_model <- glm(y ~ x, family = poisson)
    
    # Return p-values for both models
    c(
        lm = summary(lm_model)$coefficients[2, 4],
        poisson = summary(pois_model)$coefficients[2, 4]
    )
}

# Run simulations
n_sims <- 5000
results <- replicate(n_sims, simulate_with_effect(n = 20, beta = 1))

# Calculate power for both methods
power_lm <- mean(results["lm",] < 0.05)
power_poisson <- mean(results["poisson",] < 0.05)

# Print results
cat("\nPower Analysis Results:\n")
cat("Linear regression power:", round(power_lm, 3), "\n")
cat("Poisson regression power:", round(power_poisson, 3), "\n")

# Visualize p-value distributions
par(mfrow = c(1, 2))
hist(results["lm",], breaks = 30, main = "Linear Regression p-values",
     xlab = "p-value", ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2)

hist(results["poisson",], breaks = 30, main = "Poisson Regression p-values",
     xlab = "p-value", ylab = "Frequency")
abline(v = 0.05, col = "red", lwd = 2)
