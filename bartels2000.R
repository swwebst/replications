##################
# STEVEN WEBSTER #
#   MLE P.S. 4   #
#     R CODE     #
##################

############
# PART ONE #
############

rm(list=ls(all=TRUE))

library(ggplot2)
library(plyr)
library(reshape2)

set.seed(54231)

# CREATING CDF'S
# cdf for the scobit
# taken from Achen (2002), pg. 428
scobitcdf <- function(z, a){
  1- (1 / ((1 + exp(z)) ^ a))
}

# cdf for the logit
logitcdf <- function(z){
  1 / (1 + exp(-z))
}

# general Bernoulli distribution
bernoulli <- function(p, y){
  (p ^ y) * ((1 - p) ^ (1 - y))
}


# QUESTION ONE
# dgp
dgp_scobit <- function(n, a, b){
  x <- rnorm(n = n, mean = 0, sd = 1) # standard normal
  xb <- cbind(1, x) %*% b # creates vector of coefficients
  p <- scobitcdf(z = xb, a = a) # gives probability of "success" in dichotomous "y"
  y <- rbinom(n = n, size = 1, prob = p) # dummy variable "y"
  d <- data.frame(y = y, x = x)
}

d <- dgp_scobit(2000, 0.4, c(3, 7)) # data frame using scobit dgp

# QUESTION TWO
# scobit
scobit_loglike <- function(param, d){
  y <- d[1]
  x <- cbind(1, as.matrix(d[2]))
  b <- param[1:2]
  a <- param[3]
  xb <- x %*% b
  p <- scobitcdf(z = xb, a = a)
  sum(log(bernoulli(p = p, y = y)))
}

# logit
logit_loglike <- function(param, d) {
  y <- d[1]
  x <- cbind(1, as.matrix(d[2]))
  b <- param[1:2]
  xb <- x %*% b
  p <- logitcdf(z = xb)
  sum(log(bernoulli(p = p, y = y)))
}


# QUESTION THREE
# maximize the scobit
max_scobit <- function(d){
                optim(
                par = c(1,1,1),
                fn = scobit_loglike,
                hessian = TRUE,
                control = list(fnscale = -1),
                d = d
                )
}
max_scobit(d)
# this gives an alpha range of 2.7-6.37
# gives a beta of .42
# it does a *decent* job of recovering the "true" values

scobitd2 <- dgp_scobit(2000, 0.2, c(2,4))
scobitd3 <- dgp_scobit(2000, 0.3, c(5,9))

# maximize the scobit again
max_scobit <- function(d){
  optim(
    par = c(1,1,1),
    fn = scobit_loglike,
    hessian = TRUE,
    control = list(fnscale = -1),
    d = d
  )
}
max_scobit(scobitd2)

# maximize the scobit one more time
max_scobit <- function(d){
  optim(
    par = c(1,1,1),
    fn = scobit_loglike,
    hessian = TRUE,
    control = list(fnscale = -1),
    d = d
  )
}
max_scobit(scobitd3)


# QUESTION FOUR
# additional data sets
# must use the dgb_scobit function made earlier
d1 <- dgp_scobit(2000, 0.4, c(1,4))
d2 <- dgp_scobit(2000, 0.5, c(1,3))
d3 <- dgp_scobit(2000, 0.8, c(1,7))

# logit function
max_logit <- function(d){
                optim(
                  par = c(1,1),
                  fn = logit_loglike,
                  hessian = TRUE,
                  control = list(fnscale = -1),
                  d = d
  )
}

# logit and scobit models with d1, d2, and d3
s1 <- max_scobit(d1) # alpha = 0.43; beta ranges = .8 - 3.95
s2 <- max_scobit(d2) # alpha = 0.43; beta ranges = 1.38 - 3.49
s3 <- max_scobit(d3) # alpha = 0.43; beta ranges = 2.59 - 9.59
l1 <- max_logit(d1) # alpha = -0.53; beta = 2.57
l2 <- max_logit(d2) # alpha = -0.10; beta = 2.26
l3 <- max_logit(d3) # alpha = 0.33; beta = 6.54

# PREDICTED PROBABILITIES
# scobit
scobitpp <- function(param, d) {
  y <- d[1]
  x <- cbind(1, as.matrix(d[2]))
  b <- param[1:2]
  a <- param[3]
  xb <- x %*% b
  scobitcdf(z = xb, a = a)
}

# logit 
logitpp <- function(param, d) {
  y <- d[1]
  x <- cbind(1, as.matrix(d[2]))
  b <- param[1:2]
  xb <- x %*% b
  logitcdf(z = xb)
}

# store estimates
s1a <- s1$par
s2a <- s2$par
s3a <- s3$par
l1a <- l1$par
l2a <- l2$par
l3a <- l3$par

# lstplot = "logit, scobit, true value plot"
lstplot <- function(d, sparameters, lparameters, tparameters) {
  truepr <- scobitpp(param = tparameters, d = d)
  scobitpr <- scobitpp(param = sparameters, d = d)
  logitpr <- logitpp(param = lparameters, d = d)
  plotdataframe <- data.frame(x = d[2], truepr, scobitpr, logitpr)
  plotdataframe <- melt(plotdataframe, id = "x")
  names(plotdataframe) <- c("x", "model", "p.hat")
  levels(plotdataframe$model) <- c("true", "scobit", "logit")
  return(plotdataframe)
}

# data sets
firstdata <- lstplot(d1, s1a, l1a, c(1,4,0.4))
seconddata <- lstplot(d2, s2a, l2a, c(1,3,0.5))
thirddata <- lstplot(d3, s3a, l3a, c(1,7,0.8))

par(mfrow=c(2,2))
ggplot(firstdata, aes(x = x, y = p.hat, colour = model)) + geom_line() + ggtitle("Predicted Probabilities")
ggplot(seconddata, aes(x = x, y = p.hat, colour = model)) + geom_line() + ggtitle("Predicted Probabilities")
ggplot(thirddata, aes(x = x, y = p.hat, colour = model)) + geom_line() + ggtitle("Predicted Probabilities")

