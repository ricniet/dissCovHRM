library(R2jags); library(dplyr)

load("../Project5_MIRT/ncme_2017_workshop/code/covHRM/covHRM_gendata1.Rdat")

data <- covHRM_data$less_data

SACD <-  read.csv("../sacd/data/SACD144A_5tp.csv")

data <- SACD %>% filter(timepoint == 0, score > 0)

type <- as.numeric(as.factor(data$type))

data$item <- as.numeric(as.factor(data$item_ID))
#head(SACD)
item<- as.numeric(as.factor(data$item))
id<- as.numeric(as.factor(data$subject))
rater<- as.numeric(as.factor(data$rater))
type <- as.numeric(as.factor(data$type))
x <- data$x

#data$type <- ifelse(data$type == "p", 0, 1)

####
###### Define Variables
####

N <- length(unique(id))     #no. of students
R <- length(unique(rater))  #no. of raters
J <- length(unique(item))   #no. of [altruism] items
K <- max(x)                 #no. of categories
NN  <- dim(data)[1]  # N*R*J*M, extent of y

####
###### Estimation
####

lhrm_cov <- function(){  
  # likelihood - sdm part
  
  for (i in 1:NN) { # NN = N*J*R*M
    
    x[i] ~ dcat(p[i,])
    
    for (k in 1:K) {
      d[i,k] <- k - xi[id[i],item[i]] - rhocov[rater[i]] #rhocov = contribution to bias by pseudorater v
      z[i,k] <- exp(-d[i,k]*d[i,k]/2*exp(zeta[rater[i]])) #zeta = contribution to variance by pseudorater v
      p[i,k] <- z[i,k]/sum(z[i,])
    }}
  
  for(nu in 1:R){      
    rhocov[nu] ~ dnorm(alpha[type[nu]],prec.rhocov)               ###NOte that ratertype is a variable that should appear in the input file. It should indicate the type of rater, 0 (parents) or 1 (teacher)
    zeta[nu] ~ dnorm(taustar[type[nu]],prec.zeta)
  }
  
  for(g in 1:2){                ####This is 1:2 because there are TWO types of raters
    alpha[g]~dnorm(0,.1)
    taustar[g]~dnorm(0,.1)
  }        
  
  
  # likelihood - pcm part
  
    for (n in 1:N) {
      for(j in 1:J) {
        xi[n,j] ~ dcat(pcm[n,j,])
        psi[n,j,1] <- 0
        for (k in 1:(K-1)) {
          psi[n,j,k+1] <- a[j]*(th[n] - g[j,k])
        }
        for (k in 1:K) {
          term[n,j,k] <- exp(sum(psi[n,j,1:k]))
        }
        for (k in 1:K) {
          pcm[n,j,k] <- term[n,j,k]/sum(term[n,j,])
        }
      }
    }
  
  # priors - pcm part
  
  for (j in 1:J) {
    a[j] ~ dgamma(A.ALPHA,A.BETA)
    for (k in 1:(K-2)) {
      g[j,k] ~ dnorm(G.MEAN,1/(G.SD*G.SD))
    }
    g[j,3] <- -sum(g[j,1:2])
  }
  
  for (n in 1:N) {
    th[n] ~ dnorm(0,th.prec) 
  }
    
  prec.rhocov ~ dgamma(1,1)
  prec.zeta ~ dgamma(1,1)
  
  # for(nu in 1:R){   
  #    rhocov[nu] ~ dnorm(alpha[type[nu]],prec.rhocov)
  #    zeta[nu] ~ dnorm(taustar[type[nu]],prec.zeta)
  #  }
  
  #  for(g in 1:2){              ### THIS IS 2 AGAIN, SO DON"T CHANGE ANYTHING  HERE
  #    alpha[g]~dnorm(0,.1)
  #    taustar[g]~dnorm(0,.1)
  #  }
  
  
  th.prec ~ dgamma(TH.PREC.ALPHA,TH.PREC.BETA)
  
  # transformations -- for theta sd and for pcm item params...
  
  sd.theta <- 1/sqrt(th.prec)
  
  #   for (m in 1:M) {
  #     mean.th.m[m] <- mean(th[,m])
  #     sd.th.m[m] <- sd(th[,m])
  #   }
}

###
##### Initial Values
###

lhrm_cov_inits <- function()
  list(
    rhocov =  rep(0, each = R),
    zeta =  rep(1, each = R),
    th.prec = rgamma(1,100,100),
    a = runif(J,.5,1.5),
    b = rnorm(J,0,1),
    g = matrix(c(rnorm(J*(K-2),0,0.25),rep(NA,J)),nrow=J,ncol=K-1)
  )

###
##### Data Read-In
###



A.ALPHA=1
A.BETA =1
B.MEAN = 0
B.SD = 1
G.MEAN = 0
G.SD = 1
TH.PREC.ALPHA = 10
TH.PREC.BETA  = 10


lhrm_cov_data <- list("id", "rater", "x", "item", "type", 
                      "NN", "N", "J", "K", "R", "A.ALPHA",
                      "A.BETA","B.MEAN","B.SD","G.MEAN","G.SD",
                      "TH.PREC.ALPHA","TH.PREC.BETA")

lhrm_cov_params <- c("a","beta","g","sd.theta","th","rhocov","zeta","taustar")

###
##### For JAGS Run
###

covHRM_jags <- jags.parallel(inits = lhrm_cov_inits, data = lhrm_cov_data, 
                                parameters.to.save = lhrm_cov_params, n.iter = 3000,
                                model.file = lhrm_cov, n.burnin = 0, n.thin = 1,
                                jags.module = 'glm',
                                n.chains = 2)

test <- data.frame(round(covHRM_jags$BUGSoutput$summary,4))
library(ggmcmc)
s <- ggs(as.mcmc(lhrm_cov_run16))
ggmcmc(s, file="~/Desktop/trace_test.pdf", plot = "traceplot")
