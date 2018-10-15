library(R2jags); library(tictoc); library(tidyverse)

source('./dataGeneration/genModels.R')

#load('./dataGeneration/generatedDataSets/forTests/set2_simDat_fcovHRM_onlyRHO_N500R50J8K5.Rdat')
#simDat <- covHRMr.onlyRho.sim(N=200,J=8,K=5,R=20,twoPL=F)
#data <- simDat$data
data <- simDat$less_data
subject <- data$subject
item <- data$item
rater <- data$rater
rMatrix <- data %>% select(rater,cov1) %>% distinct()
rater2 <- rMatrix$rater
cov1 <- rMatrix$cov1
x <- data$x

N <- simDat$N     
S <- 2  #no. of covariate levels
R <- simDat$R 
J <- simDat$J 
K <- simDat$K                
NN<- nrow(data)
sd.phi <- 1
sd.psi <- sqrt(15)
sd.eta <- 1
sd.beta <- 1
sd.kappa <- 1

covHRM_mod = function(){
  
  # Specify Likelihood
  
  ## Signal Detection Model
  
  for (i in 1:NN) { 
    
    x[i] ~ dcat(prob.sdt[i,])
    
    for (k in 1:K) {
      
      a[i,k] <- k - xi[subject[i],item[i]] - rho[rater[i]] 
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
  
  for(nu in 1:R){      
    rho[nu] ~ dnorm(rho.term[nu], prec.rho)
    rho.term[nu] <- phi[rater2[nu]] + eta[cov1[nu]]
  }
  
  for (r in 1:R) {
    psi[r] ~ dt(0, prec.psi, 1);T(0,)
  }
  
  # CHANGE CHANGE:
  phi[48] <- -(phi[1]+phi[3]+phi[5]+phi[7]+phi[8]+phi[9]+phi[12]+phi[16]+
                 phi[18]+phi[21]+phi[23]+phi[25]+phi[28]+phi[29]+phi[33]+
                 phi[34]+phi[35]+phi[36]+phi[37]+phi[38]+phi[40]+phi[41]+
                 phi[43]+phi[44]+phi[45]+phi[46]+phi[47])
  
  phi[50] <- -(phi[2]+phi[4]+phi[6]+phi[10]+phi[11]+phi[13]+phi[14]+phi[15]+
                 phi[17]+phi[19]+phi[20]+phi[22]+phi[24]+phi[26]+phi[27]+
                 phi[30]+phi[31]+phi[32]+phi[39]+phi[42]+phi[49])
  
  # CHANGE CHANGE:
  for (nu in 1:(47)){      
    phi[nu] ~ dnorm(0, prec.phi)
  }
  
  phi[49] ~ dnorm(0, prec.phi)
  
  for(s in 1:S){                
    eta[s] ~ dnorm(0, prec.eta)
  }
  
  prec.rho <- pow(sd.rho,-2)
  prec.psi <- pow(sd.psi,-2)
  prec.phi <- pow(sd.phi,-2)
  prec.eta <- pow(sd.eta,-2)
  
  sd.rho ~ dt(0,1/15,1);T(0,)
  
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
  
  # save extra parameters
  
  etaDelta <- eta[1] - eta[2]
  
  # transformations
  
  #prec.alpha <- pow(sd.alpha, -2)
  prec.beta <- pow(sd.beta, -2)
  prec.kappa <- pow(sd.kappa, -2)
}

covHRM_data <- list("subject","rater","rater2","x","item","cov1",
                    "NN","N","J","R","K","S",
                    "sd.phi","sd.psi","sd.eta","sd.beta","sd.kappa")

# Specify Initial Values

covHRM_inits <- function()
  list(
    beta = c(rnorm(J-1,0,1.2),NA),
    kappa = c(NA,sort(rnorm(K-1,0,1.5))),
    theta = rnorm(N, 0, 2),
    # CHANGE CHANGE:
    phi = c(rnorm(47,0,1),NA,rnorm(1,0,1),NA),
    eta = rnorm(2,0,1.2),
    sd.rho = runif(1,.01,1),
    psi = runif(R,.5,1.5))

# Parameters to Track and Save

covHRM_params <- c("kappa","theta","beta",
                   "phi","eta","psi","rho","etaDelta",
                   'sd.rho')
# Specify MCMC run

tic('JAGS Runtime')

covHRM_jags <- jags.parallel(
  model.file=covHRM_mod,
  data=covHRM_data,
  inits=covHRM_inits,
  n.chains=3,
  n.burnin=10000,
  n.iter=20000,
  n.thin=20,
  parameters.to.save=covHRM_params,
  jags.module = 'glm', 
  jags.seed = 11071987)
toc()

covHRM_results <- data.frame(round(covHRM_jags$BUGSoutput$summary,4))

write.csv(covHRM_results, file='./justTests/test8/covHRMr_RHOonly_test8_N500R50J8K5_3ch.csv')
save(covHRM_jags, file='./justTests/test8/covHRMr_RHOonly_test8_N500R50J8K5_3ch.Rdat')

source('./justTests/assessTests.R')

assess(jagsObj = covHRM_jags, 
       R = R,
       filePath = './justTests/test8/',
       eval = c(1:47,49),
       random = T)
