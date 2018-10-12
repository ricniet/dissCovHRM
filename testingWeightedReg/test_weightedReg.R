library(tidyverse)

set.seed(666)

N <- 100 # number of observations

# Aggregated data
aggregated <- 
  data.frame(x=1:5) %>%
  mutate( y = round(2 * x + 2 + rnorm(length(x)) ),
          freq = as.numeric(table(sample(1:5, N, 
                                         replace=TRUE, prob=c(.3, .4, .5, .4, .3))))
  )

individuals <- aggregated[ rep(1:5, aggregated$freq) , c("x", "y") ]

wei_lm = lm(y ~ x, data=aggregated, weight=freq)
no_weigh_lm = lm(y ~ x, data=individuals)

wei_lm$coefficients
no_weigh_lm$coefficients

## in jags

x <- aggregated$x
y <- aggregated$y
weight <- aggregated$freq
N <- length(y)

x2 <- individuals$x
y2 <- individuals$y
N2 <- length(individuals$y)

lin_wt_mod <- function() {
  
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau*weight[i])
    mu[i] <- beta[1] + beta[2] * x[i]
  }
  
  # Prior for beta
  for(j in 1:2){
    beta[j] ~ dnorm(0,0.0001)
  }
  
  # Prior for the inverse variance
  tau   ~ dgamma(0.001, 0.001)
  sigma     <- 1/sqrt(tau)
}

lin_mod <- function() {
  
  for (i in 1:N2) {
    y2[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta[1] + beta[2] * x2[i]
  }
  
  # Prior for beta
  for(j in 1:2){
    beta[j] ~ dnorm(0,0.0001)
  }
  
  # Prior for the inverse variance
  tau   ~ dgamma(0.001, 0.001)
  sigma     <- 1/sqrt(tau)
}

dat <- list("N","x","y","weight")
dat2 <- list("N2","x2","y2")
params <- c("beta","tau","sigma")

library(R2jags)
fit_wt_lm1 <- jags.parallel(data = dat, parameters.to.save = params, model.file = lin_wt_mod,
                n.chains = 3, n.iter = 3000, n.burnin = 1000, n.thin = 1, DIC = F)
fit_wt_lm1$BUGSoutput$summary

fit_lm1 <- jags.parallel(data = dat2, parameters.to.save = params, model.file = lin_mod,
                         n.chains = 3, n.iter = 3000, n.burnin = 1000, n.thin = 1, DIC = F)
fit_lm1$BUGSoutput$summary
