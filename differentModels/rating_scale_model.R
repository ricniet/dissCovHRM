library(R2jags); library(tictoc)

source('./dataGeneration/genModels.R')

rsmModel <- function () {
  
  for (i in 1:N) {
    for (j in 1:J) {
      Y[i,j] ~ dcat(p[i,j,1:K])
    }
  }
  
  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:K) {
        eta[i,j,k] <- alpha * (theta[i] - (beta[j] + kappa[k]))
        psum[i,j,k] <- sum(eta[i,j,1:k])
        exp.psum[i,j,k] <- exp(psum[i,j,k])
        p[i,j,k] <- exp.psum[i,j,k] / sum(exp.psum[i,j,1:K])
      }
    }
  }
  
  # prior: theta
  
  for (i in 1:N) {
    theta[i] ~ dnorm(0, 1)
  }
  
  kappa[1] <- 0
  
  for (k in 2:K) {
    kappa[k] ~ dnorm(0, prec.kappa)
  }
  
  beta[J] <- -sum(beta[1:(J-1)])
  
  for (j in 1:(J-1)) {
    beta[j] ~ dnorm(0, prec.beta)
  }
  
  alpha <- 1
  
  # transformations
  
  #prec.alpha <- pow(sd.alpha, -2)
  prec.beta <- pow(sd.beta, -2)
  prec.kappa <- pow(sd.kappa, -2)
}

simDat <- rsm.sim(N=300,J=18,K=5,twoPL=F)

Y <- simDat$Y
N <- nrow(Y)
J <- ncol(Y)
K <- length(unique(Y[,1]))
sd.alpha <- 1
sd.beta <- 1
sd.kappa <- 1

rsmDat <- list("Y","N","J","K","sd.beta","sd.kappa")

rsmParams <- c('theta','alpha','beta','kappa')

tic('JAGS Runtime')
jagsFit <- jags.parallel(
    data = rsmDat,
    parameters.to.save = rsmParams,
    model.file = rsmModel,
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 2,
    jags.module = 'glm'
  )
toc()
