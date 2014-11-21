######################
#   Steven Webster   # 
# MLE Referee Report #
#   Bartels (2000)   #
######################

library(foreign)
library(ggplot2)
library(plyr)

d <- read.dta(file.choose())

#########################
# REPLICATING TABLE ONE #
#########################

# Note: The analysis used in Batels (2000)
# excludes those who did note vote and 
# those who voted for a third-party candidate.

# creating variables for the model
# independent variables
d$strong[d$VCF0301=="7. Strong Republican"] <- 1
d$strong[d$VCF0301=="1. Strong Democrat"] <- -1
d$strong[d$VCF0301=="2. Weak Democrat"|d$VCF0301=="3. Independent - Democrat"|d$VCF0301=="4. Independent - Independent"|d$VCF0301=="5. Independent - Republican"|d$VCF0301=="6. Weak Republican"] <- 0

d$weak[d$VCF0301=="6. Weak Republican"] <- 1
d$weak[d$VCF0301=="2. Weak Democrat"] <- -1
d$weak[d$VCF0301=="7. Strong Republican"|d$VCF0301=="1. Strong Democrat"|d$VCF0301=="3. Independent - Democrat"|d$VCF0301=="4. Independent - Independent"|d$VCF0301=="5. Independent - Republican"] <- 0

d$leaners[d$VCF0301=="5. Independent - Republican"] <- 1
d$leaners[d$VCF0301=="3. Independent - Democrat"] <- -1
d$leaners[d$VCF0301=="1. Strong Democrat"|d$VCF0301=="7. Strong Republican"|d$VCF0301=="2. Weak Democrat"|d$VCF0301=="6. Weak Republican"|d$VCF0301=="4. Independent - Independent"] <- 0

# dependent variable
d$VCF0706 <- as.numeric(d$VCF0706)
d$VCF0706[d$VCF0706==1|d$VCF0706==4|d$VCF0706==5|d$VCF0706==6] <- NA
na.omit(d)

d$votegop <- as.numeric(d$VCF0706==3)

# bartels data
# must subset because Bartels (2000) only
# uses 1952-1996, while the dataset contains
# 1948 - 2012.
bartels.dat <- subset(d, VCF0004 > 1948 & VCF0004 < 1998)

# model = probit(votegop ~ strong + weak + leaners)

# bernoulli likelihood
bern.lik <- function(p,y) {
    lik <- p^y * (1-p)^(1-y)
    return(prod(lik))
}

# bernoulli log likelihood
bern.llik <- function(p,y) {
    llik <- (y * log(p)) + ((1-y) * log(1-p))
    return(sum(llik))
}

# probit cdf
probit.cdf <- function(param, d) {
  y <- d[1]
  x <- cbind(1, as.matrix(d[2]))
  b <- param[1:2]
  xb <- x %*% b
  p <- probit.cdf(z = xb)
  sum(log(bernoulli(p = p, y = y)))
}

# probit log-likelihood
probit.loglike <- function(param, d){
	y <- d[1]
	x <- cbind(1, as.matrix(d[2]))
	b <- param[1:2]
	xb <- x %*% b 
	p <- probit.cdf(z = xb)
	sum(log(bernoulli(p = p, y = y)))
}


# probit
bartels.probit <- optim(par=0.5,
             data=x,
             fn=probit.loglike,
             method="L-BFGS-B",
             hessian=T,
             control=list(fnscale=-1),
             lower=0.000001,
             upper=0.999999)