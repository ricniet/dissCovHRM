library(R2jags); library(tictoc); library(tidyverse)

source('./dataGeneration/genModels.R')

simDat <- covHRMf.onlyRho.sim(N=500,J=8,K=5,R=40,twoPL=F,eta=c(.88,1.1))
#data <- simDat$data
data <- simDat$less_data
subject <- data$subject
item <- data$item
rater <- data$rater
cov1 <- data$cov1
x <- data$x
rMatrix <- simDat$raterMatrix %>% arrange(cov1,rater)
#for constraints:
paste('phi[',rMatrix[which(rMatrix$cov1==1),]$rater[1:(length(rMatrix[which(rMatrix$cov1==1),]$rater)-1)],']',sep='',collapse='+')
paste('phi[',rMatrix[which(rMatrix$cov1==2),]$rater[1:(length(rMatrix[which(rMatrix$cov1==2),]$rater)-1)],']',sep='',collapse='+')
# number of subjects per rater
data %>% group_by(rater) %>% summarize(n_distinct(subject))

N <- simDat$N     
S <- 2  #no. of covariate levels
R <- simDat$R 
J <- simDat$J 
K <- simDat$K                
NN<- nrow(data)
sd.phi <- sqrt(3)
sd.psi <- sqrt(15)
sd.eta <- sqrt(3)
sd.beta <- 1
sd.kappa <- 1

covHRM_mod = function(){
  
  # Specify Likelihood
  
  ## Signal Detection Model
  
  for (i in 1:NN) { 
    
    x[i] ~ dcat(prob.sdt[i,])
    
    for (k in 1:K) {
      
      a[i,k] <- k - xi[subject[i],item[i]] - (phi[rater[i]] + eta[cov1[i]])
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
        etaI[i,j,k] <- alpha * (theta[i] - (beta[j] + kappa[k]))
        psum[i,j,k] <- sum(etaI[i,j,1:k])
        exp.psum[i,j,k] <- exp(psum[i,j,k])
        prob.irt[i,j,k] <- exp.psum[i,j,k] / sum(exp.psum[i,j,1:K])
      }
    }
  }
  
  # Priors, Constraints, and Transformations
  
  # Rater and covariate effects
  
  phi[19] <- -(phi[2] + phi[5] + phi[6] + phi[8] + phi[10] + phi[11] + phi[16] +
                 phi[18])
  
  phi[20] <- -(phi[1] + phi[3] + phi[4] + phi[7] + phi[9] + phi[12] + phi[13] +
                 phi[14] + phi[15] + phi[17])
  
  for(nu in 1:(R-2)){      
    phi[nu] ~ dnorm(0, prec.phi)
  }
  
  for (r in 1:R) {
    psi[r] ~ dt(0, prec.psi, 1);T(0,)
  }
  
  for(s in 1:S){                
    eta[s] ~ dnorm(0, prec.eta)
  }
  
  prec.psi <- pow(sd.psi,-2)
  prec.phi <- pow(sd.phi,-2)
  prec.eta <- pow(sd.eta,-2)
  
  #IRT parameters
  
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

covHRM_data <- list("subject","rater","x","item","cov1",
                    "NN","N","J","R","K","S",
                    "sd.phi","sd.psi","sd.eta","sd.beta","sd.kappa")

# Specify Initial Values

covHRM_inits <- function()
  list(
    rhocov =  rnorm(R,0,.5),
    zeta =  rep(1, each = R),
    theta.prec = rgamma(1,100,100),
    alpha = c(runif(J-1,.5,1.5),NA),
    b = rnorm(J,0,.5),
    g = matrix(rnorm(J*(K-1),0,0.25),nrow=J,ncol=K-1),
    theta = rnorm(N, 0, 1))

# Parameters to Track and Save

covHRM_params <- c("alpha","kappa","theta",
                   "phi","eta","psi")

# Specify MCMC run

tic('JAGS Runtime')

covHRM_jags <- jags.parallel(
    model.file=covHRM_mod,
    data=covHRM_data,
    #inits=covHRM_inits,
    n.chains=3,
    n.burnin=2500,
    n.iter=5000,
    n.thin=5,
    parameters.to.save=covHRM_params,
    jags.module = 'glm', 
    jags.seed = 11071987)
toc()

covHRM_results <- data.frame(round(covHRM_jags$BUGSoutput$summary,4))

save(simDat, file='./testResults/fcovHRM_onlyRHO/simDat_fcovHRM_onlyRHO_N300R20J12K5_5000i2500b1t.Rdat')
write.csv(covHRM_results, file='./testResults/fcovHRM_onlyRHO/N300R20J12K5_5000i2500b1t.csv')

save(covHRM_jags, file='./testResults/fcovHRM_onlyRHO/N300R20J12K5_5000i2500b1t.Rdat')
