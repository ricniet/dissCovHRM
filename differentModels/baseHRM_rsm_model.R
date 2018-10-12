library(R2jags); library(tictoc)

source('./dataGeneration/genModels.R')

baseHRM <- function () {
  
  ## Signal Detection Model
  
  for (i in 1:NN) { 
    
    x[i] ~ dcat(prob.sdt[i,])
    
    for (k in 1:K) {
      
      a[i,k] <- k - xi[subject[i],item[i]] - phi[rater[i]] 
      b[i,k] <- -pow(a[i,k],2)
      c[i,k] <- 1 / (2 * psi[rater[i]]^2)
      d[i,k] <- exp(b[i,k] * c[i,k])
      prob.sdt[i,k] <- d[i,k]/sum(d[i,])
    }}
  
  # RSM
  
  for (i in 1:N) {
    for (j in 1:J) {
      
      xi[i,j] ~ dcat(prob.irt[i,j,1:K])
      
      for (k in 1:K) {
        eta[i,j,k] <- alpha * (theta[i] - (beta[j] + kappa[k]))
        psum[i,j,k] <- sum(eta[i,j,1:k])
        exp.psum[i,j,k] <- exp(psum[i,j,k])
        prob.irt[i,j,k] <- exp.psum[i,j,k] / sum(exp.psum[i,j,1:K])
      }
    }
  }
  
  # Priors, Constraints, and Transformations
  
  # prior: rater stuff
  
  for (r in 1:R) {
    phi[r] ~ dnorm(0, prec.phi)
    psi[r] ~ dt(0,prec.psi,1);T(0,)
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
  prec.phi <- pow(sd.phi, -2)
  prec.psi <- pow(sd.psi, -2)
}

simDat <- baseHRM.sim(N=400,J=8,K=5,R=20,twoPL=F)

data <- simDat$less_data
N <- simDat$N
R <- simDat$R
K <- simDat$K
J <- simDat$J
NN <- nrow(data)
subject <- data$subject
x <- data$x
rater <- data$rater
item <- data$item
sd.beta <- 1
sd.kappa <- 1
sd.phi <- 1
sd.psi <- sqrt(15)

baseHRM_data <- list("subject","rater","x","item",
                    "NN","N","J","R","K","sd.beta",
                    "sd.kappa","sd.phi","sd.psi")

baseHRM_params <- c("beta","kappa","theta","phi","psi")

tic('JAGS Runtime')

baseHRM_jags <- jags.parallel(
    model.file=baseHRM,
    data=baseHRM_data,
    #inits=covHRM_inits,
    n.chains=2,
    n.burnin=3000,
    n.iter=6000,
    n.thin=6,
    parameters.to.save=baseHRM_params,
    jags.module = 'glm', 
    jags.seed = 11071987)
toc()

baseHRM_results <- data.frame(round(baseHRM_jags$BUGSoutput$summary,4))
