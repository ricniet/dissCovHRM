library(R2jags); library(tictoc); library(tidyverse)

source('./dataGeneration/genModels.R')

load('./dataGeneration/generatedDataSets/forTests/set2_simDat_fcovHRM_onlyRHO_N500R50J8K5.Rdat')
#simDat <- covHRMf.onlyRho.sim(N=500,J=8,K=5,R=40,twoPL=F,eta=c(.88,1.1))
data <- simDat$data
#data <- simDat$less_data
subject <- data$subject
item <- data$item
rater <- data$rater
cov1 <- data$cov1
x <- data$x
rMatrix <- simDat$raterMatrix %>% arrange(cov1,rater)
#for constraints:
paste('phi[',rMatrix[which(rMatrix$cov1==1),]$rater[1:(length(rMatrix[which(rMatrix$cov1==1),]$rater)-1)],']',sep='',collapse='+')
paste('phi[',rMatrix[which(rMatrix$cov1==2),]$rater[1:(length(rMatrix[which(rMatrix$cov1==2),]$rater)-1)],']',sep='',collapse='+')
#for saving rho
paste('rho[',rMatrix[which(rMatrix$cov1==1),1],'] <- ','phi[',rMatrix[which(rMatrix$cov1==1),1],'] + eta[1]',collapse=' ',sep='')
paste('rho[',rMatrix[which(rMatrix$cov1==2),1],'] <- ','phi[',rMatrix[which(rMatrix$cov1==2),1],'] + eta[2]',collapse=' ',sep='')
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
  
  # save extra parameters
  
  etaDelta <- eta[1] - eta[2]
  
  # CHANGE CHANGE:
  rho[1] <- phi[1] + eta[1] 
  rho[3] <- phi[3] + eta[1] 
  rho[5] <- phi[5] + eta[1] 
  rho[7] <- phi[7] + eta[1] 
  rho[8] <- phi[8] + eta[1] 
  rho[9] <- phi[9] + eta[1] 
  rho[12] <- phi[12] + eta[1] 
  rho[16] <- phi[16] + eta[1] 
  rho[18] <- phi[18] + eta[1] 
  rho[21] <- phi[21] + eta[1] 
  rho[23] <- phi[23] + eta[1] 
  rho[25] <- phi[25] + eta[1] 
  rho[28] <- phi[28] + eta[1] 
  rho[29] <- phi[29] + eta[1] 
  rho[33] <- phi[33] + eta[1] 
  rho[34] <- phi[34] + eta[1] 
  rho[35] <- phi[35] + eta[1] 
  rho[36] <- phi[36] + eta[1] 
  rho[37] <- phi[37] + eta[1] 
  rho[38] <- phi[38] + eta[1] 
  rho[40] <- phi[40] + eta[1] 
  rho[41] <- phi[41] + eta[1] 
  rho[43] <- phi[43] + eta[1] 
  rho[44] <- phi[44] + eta[1] 
  rho[45] <- phi[45] + eta[1] 
  rho[46] <- phi[46] + eta[1] 
  rho[47] <- phi[47] + eta[1] 
  rho[48] <- phi[48] + eta[1]
  
  rho[2] <- phi[2] + eta[2] 
  rho[4] <- phi[4] + eta[2] 
  rho[6] <- phi[6] + eta[2] 
  rho[10] <- phi[10] + eta[2] 
  rho[11] <- phi[11] + eta[2] 
  rho[13] <- phi[13] + eta[2] 
  rho[14] <- phi[14] + eta[2] 
  rho[15] <- phi[15] + eta[2] 
  rho[17] <- phi[17] + eta[2] 
  rho[19] <- phi[19] + eta[2] 
  rho[20] <- phi[20] + eta[2] 
  rho[22] <- phi[22] + eta[2] 
  rho[24] <- phi[24] + eta[2] 
  rho[26] <- phi[26] + eta[2] 
  rho[27] <- phi[27] + eta[2] 
  rho[30] <- phi[30] + eta[2] 
  rho[31] <- phi[31] + eta[2] 
  rho[32] <- phi[32] + eta[2] 
  rho[39] <- phi[39] + eta[2] 
  rho[42] <- phi[42] + eta[2] 
  rho[49] <- phi[49] + eta[2] 
  rho[50] <- phi[50] + eta[2]
  
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
    beta = c(rnorm(J-1,0,1.2),NA),
    kappa = c(NA,sort(rnorm(K-1,0,1.5))),
    theta = rnorm(N, 0, 2),
    # CHANGE CHANGE:
    phi = c(rnorm(47,0,1),NA,rnorm(1,0,1),NA),
    eta = rnorm(2,0,1.2),
    psi = runif(R,.5,1.5))

# Parameters to Track and Save

covHRM_params <- c("kappa","theta","beta",
                   "phi","eta","psi","rho","etaDelta")

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

write.csv(covHRM_results, file='./justTests/test7/covHRMf_RHOonly_test7_N500R50J8K5_3ch20000i10000b20t.csv')
save(covHRM_jags, file='./justTests/test7/covHRMf_RHOonly_test7_N500R50J8K5_3ch20000i10000b20t.Rdat')

source('./justTests/assessTests.R')

assess(jagsObj = covHRM_jags, 
       R = R,
       filePath = './justTests/test7/',
       eval = c(1:47,49),
       random = F)
