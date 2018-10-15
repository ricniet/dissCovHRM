args <- commandArgs(trailingOnly = TRUE)
name <- ifelse(length(args)>0, args[1], "nodelist.txt")
testNum <- args[2]
modelType <- args[3]
dataDesc <- args[4]
simDesc <- args[5]

assign("hostnames", scan(name, what="", sep="\n"), envir = .GlobalEnv)
library(R2jags)
library(plyr)
library(reshape2)

# Function D: Bayesian/MCMC model -----------------------------------------

est_covHRM <- function(data) {

  subject <- data$subject
  item <- data$item
  rater <- data$rater
  cov1 <- data$cov1
  x <- data$x
  
  N <- length(unique(subject))     
  S <- 2  #no. of covariate levels
  R <- length(unique(rater))
  J <- length(unique(item))
  K <- length(unique(x))                
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
    
    phi[30] <- -(phi[2]+phi[5]+phi[7]+phi[8]+phi[9]+phi[10]+phi[11]+
                   phi[15]+phi[16]+phi[17]+phi[18]+phi[22]+phi[25]+phi[28]+phi[29])
    
    phi[27] <- -(phi[1]+phi[3]+phi[4]+phi[6]+phi[12]+phi[13]+phi[14]+
                   phi[19]+phi[20]+phi[21]+phi[23]+phi[24]+phi[26])
    
    for (nu in 1:(26)){      
      phi[nu] ~ dnorm(0, prec.phi)
    }
    
    for (nu in 28:(29)){      
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
    
    # save extra parameters
    
    etaDelta <- eta[1] - eta[2]
    
    rho[2] <- phi[2] + eta[1] 
    rho[5] <- phi[5] + eta[1] 
    rho[7] <- phi[7] + eta[1] 
    rho[8] <- phi[8] + eta[1] 
    rho[9] <- phi[9] + eta[1] 
    rho[10] <- phi[10] + eta[1] 
    rho[11] <- phi[11] + eta[1] 
    rho[15] <- phi[15] + eta[1] 
    rho[16] <- phi[16] + eta[1] 
    rho[17] <- phi[17] + eta[1] 
    rho[18] <- phi[18] + eta[1] 
    rho[22] <- phi[22] + eta[1] 
    rho[25] <- phi[25] + eta[1] 
    rho[28] <- phi[28] + eta[1] 
    rho[29] <- phi[29] + eta[1] 
    rho[30] <- phi[30] + eta[1]
    rho[1] <- phi[1] + eta[2] 
    rho[3] <- phi[3] + eta[2] 
    rho[4] <- phi[4] + eta[2] 
    rho[6] <- phi[6] + eta[2] 
    rho[12] <- phi[12] + eta[2] 
    rho[13] <- phi[13] + eta[2] 
    rho[14] <- phi[14] + eta[2] 
    rho[19] <- phi[19] + eta[2] 
    rho[20] <- phi[20] + eta[2] 
    rho[21] <- phi[21] + eta[2] 
    rho[23] <- phi[23] + eta[2] 
    rho[24] <- phi[24] + eta[2] 
    rho[26] <- phi[26] + eta[2] 
    rho[27] <- phi[27] + eta[2]
    
    # transformations
    
    #prec.alpha <- pow(sd.alpha, -2)
    prec.beta <- pow(sd.beta, -2)
    prec.kappa <- pow(sd.kappa, -2)
  }
  
  
  # Specify Initial Values
  
  covHRM_inits <- function()
    list(
      beta = c(rnorm(J-1,0,1.2),NA),
      kappa = c(NA,sort(rnorm(K-1,0,1.5))),
      theta = rnorm(N, 0, 2),
      phi = c(rnorm(26,0,1),NA,rnorm(2,0,1),NA),
      eta = rnorm(2,0,1.2),
      psi = runif(R,.5,1.5))
  
  covHRM_data <- list(subject=subject,rater=rater,x=x,item=item,cov1=cov1,
                      NN=NN,N=N,J=J,R=R,K=K,S=S,
                      sd.phi=sd.phi,sd.psi=sd.psi,sd.eta=sd.eta,sd.beta=sd.beta,sd.kappa=sd.kappa)
  
  assign("n.chains.global", 3, envir = .GlobalEnv)
  
  # Parameters to Track and Save
  
  covHRM_params <- c("kappa","theta","beta",
                     "phi","eta","psi","rho","etaDelta")
  
   library(dclone)
   
   cl <- makePSOCKcluster(hostnames)
   library(jags.parallel.tacc)
  
   covHRM_tacc <- jags.par(cl, #jags.par
                         data = covHRM_data,
                         inits = covHRM_inits,
                         parameters.to.save = covHRM_params,
                         model.file = covHRM_mod,
                         n.iter=20000,n.chains=n.chains.global,n.thin=20,n.burnin=10000
   )
  return(covHRM_tacc)
}

# Function F: Simulation driver, calls C/D --------------------------------

run_sim <- function(){
  data <- read.csv('set1_simDat_fcovHRM_onlyRHO_N500R40J8K5.csv')
  
  jagsObj <- est_covHRM(data=data) #estimation function
  newName <- paste0('jagsObj_test',testNum,'_',modelType,'_',dataDesc,'_',simDesc,sep = '')
  assign(newName, jagsObj)
  
  ## this code was added by TACC folks
  fileName=paste(newName,".Rdat",sep="")
  save(newName,data,file=fileName)
  
  #save(jags.hrm.manytimes,data_all,file=paste("mirtHRM_",rindex,".Rdat",sep=""))
  stat = data.frame(round(jagsObj$BUGSoutput$summary,4))
  names(stat)[5] <- "median"
  
  newName=paste0('jagsSummary_test',testNum,'_',modelType,'_',dataDesc,'_',simDesc,sep = '')
  fileName=paste(newName,".csv",sep="")
  write.csv(stat,fileName)
}

###############
## Main Call ##
###############

pt=proc.time()
v=run_sim()
print(proc.time()-pt)