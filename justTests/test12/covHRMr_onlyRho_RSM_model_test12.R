library(R2jags); library(tictoc); library(tidyverse)

#source('./dataGeneration/genModels.R')

load('./dataGeneration/generatedDataSets/forTests/set4_simDat_rcovHRM_onlyRHO_N600R75J8K5.Rdat')
#simDat <- covHRMr.onlyRho.sim(N=200,J=8,K=5,R=20,twoPL=F)
#data <- simDat$data
data <- simDat$less_data
subject <- data$subject
item <- data$item
rater <- data$rater
rMatrix <- simDat$raterMatrix %>% arrange(cov1,rater)
rater2 <- rMatrix$rater
cov1 <- rMatrix$cov1
x <- data$x

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
  phi[73] <- -(phi[4]+phi[7]+phi[8]+phi[9]+phi[10]+phi[14]+phi[15]+phi[16]+phi[18]+
                 phi[19]+phi[20]+phi[21]+phi[22]+phi[24]+phi[25]+phi[30]+phi[31]+
                 phi[33]+phi[34]+phi[35]+phi[36]+phi[37]+phi[38]+phi[39]+phi[42]+
                 phi[44]+phi[46]+phi[47]+phi[54]+phi[57]+phi[58]+phi[62]+phi[64]+
                 phi[66]+phi[69]+phi[71])
  
  phi[75] <- -(phi[1]+phi[2]+phi[3]+phi[5]+phi[6]+phi[11]+phi[12]+phi[13]+phi[17]+
                 phi[23]+phi[26]+phi[27]+phi[28]+phi[29]+phi[32]+phi[40]+phi[41]+
                 phi[43]+phi[45]+phi[48]+phi[49]+phi[50]+phi[51]+phi[52]+phi[53]+
                 phi[55]+phi[56]+phi[59]+phi[60]+phi[61]+phi[63]+phi[65]+phi[67]+
                 phi[68]+phi[70]+phi[72]+phi[74])
  
  # CHANGE CHANGE:
  for (nu in 1:(72)){      
    phi[nu] ~ dnorm(0, prec.phi)
  }
  
  phi[74] ~ dnorm(0, prec.phi)
  
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
  
  # rho2
  
  rho2[4] <- phi[4] + eta[1] 
  rho2[7] <- phi[7] + eta[1] 
  rho2[8] <- phi[8] + eta[1] 
  rho2[9] <- phi[9] + eta[1] 
  rho2[10] <- phi[10] + eta[1] 
  rho2[14] <- phi[14] + eta[1] 
  rho2[15] <- phi[15] + eta[1] 
  rho2[16] <- phi[16] + eta[1] 
  rho2[18] <- phi[18] + eta[1] 
  rho2[19] <- phi[19] + eta[1] 
  rho2[20] <- phi[20] + eta[1] 
  rho2[21] <- phi[21] + eta[1] 
  rho2[22] <- phi[22] + eta[1] 
  rho2[24] <- phi[24] + eta[1] 
  rho2[25] <- phi[25] + eta[1] 
  rho2[30] <- phi[30] + eta[1] 
  rho2[31] <- phi[31] + eta[1] 
  rho2[33] <- phi[33] + eta[1] 
  rho2[34] <- phi[34] + eta[1] 
  rho2[35] <- phi[35] + eta[1] 
  rho2[36] <- phi[36] + eta[1] 
  rho2[37] <- phi[37] + eta[1] 
  rho2[38] <- phi[38] + eta[1] 
  rho2[39] <- phi[39] + eta[1] 
  rho2[42] <- phi[42] + eta[1] 
  rho2[44] <- phi[44] + eta[1] 
  rho2[46] <- phi[46] + eta[1] 
  rho2[47] <- phi[47] + eta[1] 
  rho2[54] <- phi[54] + eta[1] 
  rho2[57] <- phi[57] + eta[1] 
  rho2[58] <- phi[58] + eta[1] 
  rho2[62] <- phi[62] + eta[1] 
  rho2[64] <- phi[64] + eta[1] 
  rho2[66] <- phi[66] + eta[1] 
  rho2[69] <- phi[69] + eta[1] 
  rho2[71] <- phi[71] + eta[1] 
  rho2[73] <- phi[73] + eta[1]
  
  rho2[1] <- phi[1] + eta[2] 
  rho2[2] <- phi[2] + eta[2] 
  rho2[3] <- phi[3] + eta[2] 
  rho2[5] <- phi[5] + eta[2] 
  rho2[6] <- phi[6] + eta[2] 
  rho2[11] <- phi[11] + eta[2] 
  rho2[12] <- phi[12] + eta[2] 
  rho2[13] <- phi[13] + eta[2] 
  rho2[17] <- phi[17] + eta[2] 
  rho2[23] <- phi[23] + eta[2] 
  rho2[26] <- phi[26] + eta[2] 
  rho2[27] <- phi[27] + eta[2] 
  rho2[28] <- phi[28] + eta[2] 
  rho2[29] <- phi[29] + eta[2] 
  rho2[32] <- phi[32] + eta[2] 
  rho2[40] <- phi[40] + eta[2] 
  rho2[41] <- phi[41] + eta[2] 
  rho2[43] <- phi[43] + eta[2] 
  rho2[45] <- phi[45] + eta[2] 
  rho2[48] <- phi[48] + eta[2] 
  rho2[49] <- phi[49] + eta[2] 
  rho2[50] <- phi[50] + eta[2] 
  rho2[51] <- phi[51] + eta[2] 
  rho2[52] <- phi[52] + eta[2] 
  rho2[53] <- phi[53] + eta[2] 
  rho2[55] <- phi[55] + eta[2] 
  rho2[56] <- phi[56] + eta[2] 
  rho2[59] <- phi[59] + eta[2] 
  rho2[60] <- phi[60] + eta[2] 
  rho2[61] <- phi[61] + eta[2] 
  rho2[63] <- phi[63] + eta[2] 
  rho2[65] <- phi[65] + eta[2] 
  rho2[67] <- phi[67] + eta[2] 
  rho2[68] <- phi[68] + eta[2] 
  rho2[70] <- phi[70] + eta[2] 
  rho2[72] <- phi[72] + eta[2] 
  rho2[74] <- phi[74] + eta[2] 
  rho2[75] <- phi[75] + eta[2]
  
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
    phi = c(rnorm(72,0,1),NA,rnorm(1,0,1),NA),
    eta = rnorm(2,0,1.2),
    sd.rho = runif(1,.01,1),
    psi = runif(R,.5,1.5))

# Parameters to Track and Save

covHRM_params <- c("kappa","theta","beta",
                   "phi","eta","psi","rho","rho2","etaDelta",
                   "sd.rho")
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

write.csv(covHRM_results, file='./justTests/test12/covHRMr_RHOonly_test12_N600R75J8K5_3ch.csv')
save(covHRM_jags, file='./justTests/test12/covHRMr_RHOonly_test12_N600R75J8K5_3ch.Rdat')

source('./justTests/assessTests.R')

assess(jagsObj = covHRM_jags, 
       R = R,
       filePath = './justTests/test9/',
       eval = c(1:47,49),
       random = T)
