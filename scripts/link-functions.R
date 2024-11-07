# based on this blog: 
# https://www.seascapemodels.org/rstats/2018/10/16/understanding-the-glm-link.html


#CJ Brown
#2024-10-10


n <- 1000
set.seed(42)

x1 <- rpois(n, lambda = 1)
x10 <- rpois(n, lambda = 10)

mean(x1)
var(x1)
mean(x10)
var(x10)

par(mfrow = c(1,2))
hist(x1, xlim = c(0, 25), seq(0, 25, by = 1))
hist(x10, xlim = c(0, 25), seq(0, 25, by = 1))


# Linear models
# ymean[i] = a + b * x[i]

n <- 100
beta <- -8 
alpha <- 4
x <- seq(0, 1, length.out = n)

ymean <- alpha + beta*x

plot(x, ymean, type = 'l', 
     xlab = "pollution",
     ylab = "Mean Number fish")
abline(h = 0, lty = 2, lwd = 2)

set.seed(55)
yobs_normal <- ymean + rnorm(n)
plot(x, ymean, type = "l")
points(x, yobs_normal)
abline(h = 0, lty = 2, lwd = 2)
yobs_normal

## Simulate from a Poisson
gamma <- -3.2 #effect of pollution
1/exp(gamma)
alpha <- 4 #intercept - mean at zero pollution

yexp <- alpha*exp(gamma*x)
#Bit of algebra
# ymean = a*exp(gamma*x)
#log(ymean) = log(a) + gamma*x

plot(x, (yexp), type = 'l',
     ylim = c(-2, 4))
abline(h = 0, lty = 2)

yobs_pois <- rpois(n, yexp)
plot(x, yexp, type = "l", 
     xlab = "pllution",
     ylab = "number fish",
     ylim = c(0, 8))
points(x, yobs_pois)

# Doing it the righ tway around 
# fit the model 

m1 <- glm(yobs_pois ~ x,
          family = poisson(link = "log"))
coef(m1)
exp(coef(m1)[1])

1/exp(coef(m1)[2])

ypredict <- predict(m1, type = "response", se = TRUE)
plot(x, yexp, type = 'l', xlab = "Pollution level",
     ylab = "Number of fish counted",
     ylim = c(0, 8))
lines(x, ypredict$fit, lwd = 2, col = "red", lty = 2)
#Add lines for standard errors
lines(x, ypredict$fit + ypredict$se.fit, lty = 3, col = "red")
lines(x, ypredict$fit - ypredict$se.fit, lty = 3, col = "red")
#plot observations
# points(x, yobs_pois)
legend('topright', legend = c("True mean", "Estimated mean"),
       lwd = 2, lty = c(1,2), col = c("black", "red"))



y <- rpois(1000, 3)
hist(y)
hist(log(y+0.1))
