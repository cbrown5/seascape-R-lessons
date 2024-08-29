# Conditional programming

library(ggplot2)
# if statement

#Set this value to TRUE to save teh plots
do_plots <- TRUE


#bunch of data wrangling
dat <- data.frame(x = 1:10,
  y = rnorm(10, x^2, 5))

# Only make the plots if do_plots = TRUE
if (do_plots) {
  g1 <- ggplot(dat) +
    geom_point(aes(x = x, y = y))
  ggsave("Outputs/plot1.png", g1)
}


#
# Another use of ifs is to run some code 
# if some condition in the
# data is met
#

dat2 <- data.frame(x = runif(30, 18, 50))
dat2$x
max_age <- 28
min_age <- 18
#logical statement to ask about values of x

if (any(dat2$x > max_age)){
  stop("Chris is not that old!")
}
#Another way to make errors is 
# stopifnot(all(dat2$x > min_age), "Chris is not that young!")

# If-else

if (any(dat2$x > max_age)){
  print("Chris is not that old!")
} else {
  hist(dat2$x)
}


#application in for loops
nrow_dat2 <- nrow(dat2)
for (i in 1:nrow_dat2){
   print(i)
  if(dat2$x[i] > max_age){
   print(dat2$x[i])
  }
 
}


