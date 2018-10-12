
covHRMr.sim <- function(N,J,K,R,twoPL=F) {
  theta <- rnorm(N)
  alpha <- ifelse(twoPL==T,rlnorm(J,0,.5),1)
  beta <- rnorm(J)
  kappa <- sort(rnorm(K-1,0,1.1))
  
  gb <- array(NA,c(J,K-1))
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      gb[j,k] <- beta[j] + kappa[k]
    }
  }
  
  term <- array(NA, dim=c(N,J,K))
  term[,,1] <- 0
  
  for (n in 1:N) {
    term[n,,2:K] <- sweep(sweep(gb,1,theta[n],"-"),1,-alpha,"*")
  }
  
  Y <- matrix(NA,nrow=N,ncol=J)
  
  for (n in 1:N) {
    for (j in 1:J) {
      term[n,j,] <- exp(cumsum(term[n,j,]))
    }
    term[n,,] <- sweep(term[n,,],1,apply(term[n,,],1,sum),"/")
    
    for (j in 1:J) {
      Y[n,j] <- sample(1:K,1,prob=term[n,j,])
    }
  }
  
  # library(eRm)
  # pcm.fit <- PCM(Y-1)
  # thresholds(pcm.fit)
  
  # SDM
  
  covLevels <- 2
  t <- data.frame(
    a = factor(1:R),
    b = factor(sample(letters[1:covLevels],size=R,replace=T,prob=rep(1/covLevels,covLevels)))
    #factor(sample(1:covLevels, size=R,replace=T,prob=rep(1/covLevels,covLevels)))
  )
  
  D <- model.matrix(~ .+0, data=t, 
                    contrasts.arg = lapply(t, contrasts, contrasts=FALSE))
  V <- nrow(D)
  
  
  library(Matrix)
  rank <- rankMatrix(D)
  
  phi <- runif(R,0.4,0.9)
  psi <- runif(R,0.1,0.7)
  eta <- c(.67,.34)
  tau <- c(.5,.2)
  etaB <- as.matrix(c(phi,eta))
  tauB <- as.matrix(c(psi,tau))
  sd.eta <- .08
  sd.tau <- .15
  
  
  rho <- rnorm(V, D %*% etaB, sd.eta)
  omega <- rnorm(V, D %*% tauB, sd.tau)
  
  xi <- c(Y) # students then items, so length is N*J
  
  tau2 <- 1/(omega^2)
  
  X <- matrix(NA,nrow=N*J,ncol=R)
  
  p <- z <- d <- array(NA,dim=c(N*J,R,K))
  
  for (i in 1:(N*J)) {
    for (r in 1:R) {
      for (k in 1:K) {
        d[i,r,k] <- k - xi[i] - rho[r]
        z[i,r,k] <- exp(-d[i,r,k]*d[i,r,k]*tau2[r]/2)
      }
      for (k in 1:K) {
        p[i,r,k] <- z[i,r,k]/sum(z[i,r,])
      }
      X[i,r] <- sample(1:K,1,prob=p[i,r,])
    }
  }
  
  x <- c(X) # subjects in order, then 
  
  subject <- rep(1:N,J*R)
  rater <- rep(1:R,rep(N*J,R))
  item <- rep(rep(1:J,rep(N,J)),R)
  cov1 <- rep(1:2, each=N*J*(R/2))
  
  data <- data.frame(x=x, cov1,subject=subject,item=item,rater=rater)
  
  ## less_data (double-scoring)
  
  assign <- as.data.frame.matrix(t(replicate(N, sample(1:R,2))))
  assign$subject <- 1:N
  names(assign)[1:2] <- paste0('rater',1:2)
  
  less_data <- merge(data, assign, by='subject')
  less_data$marker1 <- with(less_data, ifelse(rater == rater1,1,0))
  less_data$marker2 <- with(less_data, ifelse(rater == rater2,1,0))
  #less_data$marker3 <- with(less_data, ifelse(rater == rater3,1,0))
  less_data$keeper <- apply(less_data[,7:8], 1, function(x) sum(x))
  less_data <- less_data[which(less_data$keeper != 0),] 
  less_data <- less_data[,1:4]
  less_data <- less_data[order(less_data$subject, less_data$rater, less_data$item),]
  
  output <- list()
  output$data <- data
  output$less_data <- less_data
  output$N <- N
  output$R <- R
  output$K <- K
  output$J <- J
  output$observedRatings <- X
  output$idealScores <- Y
  output$trueValues <- list(
    alpha = alpha,
    beta = beta,
    kappa = kappa,
    betaKappa = gb,
    theta = theta,
    phi = phi,
    psi = psi,
    eta = eta,
    tau = tau,
    rho = rho,
    omega = omega,
    sd.eta = sd.eta,
    sd.tau = sd.tau
  )
  output
}

covHRMf.sim <- function(N,J,K,R,twoPL=F) {
  theta <- rnorm(N)
  alpha <- ifelse(twoPL==T,rlnorm(J,0,.5),1)
  beta <- rnorm(J)
  kappa <- sort(rnorm(K-1,0,1.1))
  
  gb <- array(NA,c(J,K-1))
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      gb[j,k] <- beta[j] + kappa[k]
    }
  }
  
  term <- array(NA, dim=c(N,J,K))
  term[,,1] <- 0
  
  for (n in 1:N) {
    term[n,,2:K] <- sweep(sweep(gb,1,theta[n],"-"),1,-alpha,"*")
  }
  
  Y <- matrix(NA,nrow=N,ncol=J)
  
  for (n in 1:N) {
    for (j in 1:J) {
      term[n,j,] <- exp(cumsum(term[n,j,]))
    }
    term[n,,] <- sweep(term[n,,],1,apply(term[n,,],1,sum),"/")
    
    for (j in 1:J) {
      Y[n,j] <- sample(1:K,1,prob=term[n,j,])
    }
  }
  
  # library(eRm)
  # pcm.fit <- PCM(Y-1)
  # thresholds(pcm.fit)
  
  # SDM
  
  covLevels <- 2
  t <- data.frame(
    a = factor(1:R),
    b = factor(sample(letters[1:covLevels],size=R,replace=T,prob=rep(1/covLevels,covLevels)))
    #factor(sample(1:covLevels, size=R,replace=T,prob=rep(1/covLevels,covLevels)))
  )
  
  D <- model.matrix(~ .+0, data=t, 
                    contrasts.arg = lapply(t, contrasts, contrasts=FALSE))
  V <- nrow(D)
  
  
  library(Matrix)
  rank <- rankMatrix(D)
  
  phi <- runif(R,0.4,0.9)
  psi <- runif(R,0.1,0.7)
  eta <- c(.67,.34)
  tau <- c(.5,.2)
  etaB <- as.matrix(c(phi,eta))
  tauB <- as.matrix(c(psi,tau))
  sd.eta <- .08
  sd.tau <- .15
  
  
  rho <- D %*% etaB
  omega <- D %*% tauB
  
  xi <- c(Y) # students then items, so length is N*J
  
  tau2 <- 1/(omega^2)
  
  X <- matrix(NA,nrow=N*J,ncol=R)
  
  p <- z <- d <- array(NA,dim=c(N*J,R,K))
  
  for (i in 1:(N*J)) {
    for (r in 1:R) {
      for (k in 1:K) {
        d[i,r,k] <- k - xi[i] - rho[r]
        z[i,r,k] <- exp(-d[i,r,k]*d[i,r,k]*tau2[r]/2)
      }
      for (k in 1:K) {
        p[i,r,k] <- z[i,r,k]/sum(z[i,r,])
      }
      X[i,r] <- sample(1:K,1,prob=p[i,r,])
    }
  }
  
  x <- c(X) # subjects in order, then 
  
  subject <- rep(1:N,J*R)
  rater <- rep(1:R,rep(N*J,R))
  item <- rep(rep(1:J,rep(N,J)),R)
  cov1 <- rep(1:2, each=N*J*(R/2))
  
  data <- data.frame(x=x, cov1,subject=subject,item=item,rater=rater)
  
  ## less_data (double-scoring)
  
  assign <- as.data.frame.matrix(t(replicate(N, sample(1:R,2))))
  assign$subject <- 1:N
  names(assign)[1:2] <- paste0('rater',1:2)
  
  less_data <- merge(data, assign, by='subject')
  less_data$marker1 <- with(less_data, ifelse(rater == rater1,1,0))
  less_data$marker2 <- with(less_data, ifelse(rater == rater2,1,0))
  #less_data$marker3 <- with(less_data, ifelse(rater == rater3,1,0))
  less_data$keeper <- apply(less_data[,7:8], 1, function(x) sum(x))
  less_data <- less_data[which(less_data$keeper != 0),] 
  less_data <- less_data[,1:4]
  less_data <- less_data[order(less_data$subject, less_data$rater, less_data$item),]
  
  output <- list()
  output$data <- data
  output$less_data <- less_data
  output$N <- N
  output$R <- R
  output$K <- K
  output$J <- J
  output$observedRatings <- X
  output$idealScores <- Y
  output$trueValues <- list(
    alpha = alpha,
    beta = beta,
    kappa = kappa,
    betaKappa = gb,
    theta = theta,
    phi = phi,
    psi = psi,
    eta = eta,
    tau = tau,
    rho = rho,
    omega = omega
  )
  output
}

covHRMr.onlyRho.sim <- function(N,J,K,R,twoPL=F,eta,sd.eta) {
  theta <- rnorm(N)
  alpha <- ifelse(twoPL==T,rlnorm(J,0,.5),1)
  beta <- rnorm(J)
  kappa <- sort(rnorm(K-1,0,1.1))
  
  gb <- array(NA,c(J,K-1))
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      gb[j,k] <- beta[j] + kappa[k]
    }
  }
  
  term <- array(NA, dim=c(N,J,K))
  term[,,1] <- 0
  
  for (n in 1:N) {
    term[n,,2:K] <- sweep(sweep(gb,1,theta[n],"-"),1,-alpha,"*")
  }
  
  Y <- matrix(NA,nrow=N,ncol=J)
  
  for (n in 1:N) {
    for (j in 1:J) {
      term[n,j,] <- exp(cumsum(term[n,j,]))
    }
    term[n,,] <- sweep(term[n,,],1,apply(term[n,,],1,sum),"/")
    
    for (j in 1:J) {
      Y[n,j] <- sample(1:K,1,prob=term[n,j,])
    }
  }
  
  # library(eRm)
  # pcm.fit <- PCM(Y-1)
  # thresholds(pcm.fit)
  
  # SDM
  
  covLevels <- 2
  t <- data.frame(
    rater = factor(1:R),
    cov1 = factor(sample(letters[1:covLevels],size=R,replace=T,prob=rep(1/covLevels,covLevels)))
    #factor(sample(1:covLevels, size=R,replace=T,prob=rep(1/covLevels,covLevels)))
  )
  
  D <- model.matrix(~ .+0, data=t, 
                    contrasts.arg = lapply(t, contrasts, contrasts=FALSE))
  V <- nrow(D)
  
  t$cov1 <- ifelse(t$cov1=='a',1,2) 
  t$rater <- as.integer(as.character(t$rater))
  
  library(Matrix)
  rank <- rankMatrix(D)
  
  phi <- runif(R,0.4,0.9)
  psi <- runif(R,0.1,0.7)
  etaB <- as.matrix(c(phi,eta))
  
  rho <- rnorm(V, D %*% etaB, sd.eta)

  xi <- c(Y) # students then items, so length is N*J
  
  tau2 <- 1/(psi^2)
  
  X <- matrix(NA,nrow=N*J,ncol=R)
  
  p <- z <- d <- array(NA,dim=c(N*J,R,K))
  
  for (i in 1:(N*J)) {
    for (r in 1:R) {
      for (k in 1:K) {
        d[i,r,k] <- k - xi[i] - rho[r]
        z[i,r,k] <- exp(-d[i,r,k]*d[i,r,k]*tau2[r]/2)
      }
      for (k in 1:K) {
        p[i,r,k] <- z[i,r,k]/sum(z[i,r,])
      }
      X[i,r] <- sample(1:K,1,prob=p[i,r,])
    }
  }
  
  x <- c(X) # subjects in order, then 
  
  subject <- rep(1:N,J*R)
  rater <- rep(1:R,rep(N*J,R))
  item <- rep(rep(1:J,rep(N,J)),R)
  
  data <- data.frame(x=x,subject=subject,item=item,rater=rater)
  
  data <- data %>% left_join(t, by='rater')
  
  ## less_data (double-scoring)
  
  assign <- as.data.frame.matrix(t(replicate(N, sample(1:R,2))))
  assign$subject <- 1:N
  names(assign)[1:2] <- paste0('rater',1:2)
  
  less_data <- merge(data, assign, by='subject')
  less_data$marker1 <- with(less_data, ifelse(rater == rater1,1,0))
  less_data$marker2 <- with(less_data, ifelse(rater == rater2,1,0))
  #less_data$marker3 <- with(less_data, ifelse(rater == rater3,1,0))
  less_data$keeper <- apply(less_data[,8:9], 1, function(x) sum(x))
  less_data <- less_data[which(less_data$keeper != 0),] 
  less_data <- less_data[,1:5]
  less_data <- less_data[order(less_data$subject, less_data$rater, less_data$item),]
  
  output <- list()
  output$data <- data
  output$less_data <- less_data
  output$N <- N
  output$R <- R
  output$K <- K
  output$J <- J
  output$observedRatings <- X
  output$idealScores <- Y
  output$raterMatrix <- t
  output$rank <- rank
  output$trueValues <- list(
    alpha = alpha,
    beta = beta,
    kappa = kappa,
    betaKappa = gb,
    theta = theta,
    phi = phi,
    psi = psi,
    eta = eta,
    rho = rho,
    sd.eta = sd.eta
  )
  output
}

covHRMf.onlyRho.sim <- function(N,J,K,R,twoPL=F,eta) {
  theta <- rnorm(N)
  alpha <- ifelse(twoPL==T,rlnorm(J,0,.5),1)
  beta <- rnorm(J)
  kappa <- sort(rnorm(K-1,0,1.1))
  
  gb <- array(NA,c(J,K-1))
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      gb[j,k] <- beta[j] + kappa[k]
    }
  }
  
  term <- array(NA, dim=c(N,J,K))
  term[,,1] <- 0
  
  for (n in 1:N) {
    term[n,,2:K] <- sweep(sweep(gb,1,theta[n],"-"),1,-alpha,"*")
  }
  
  Y <- matrix(NA,nrow=N,ncol=J)
  
  for (n in 1:N) {
    for (j in 1:J) {
      term[n,j,] <- exp(cumsum(term[n,j,]))
    }
    term[n,,] <- sweep(term[n,,],1,apply(term[n,,],1,sum),"/")
    
    for (j in 1:J) {
      Y[n,j] <- sample(1:K,1,prob=term[n,j,])
    }
  }
  
  # library(eRm)
  # pcm.fit <- PCM(Y-1)
  # thresholds(pcm.fit)
  
  # SDM
  
  covLevels <- 2
  t <- data.frame(
    rater = factor(1:R),
    cov1 = factor(sample(letters[1:covLevels],size=R,replace=T,prob=rep(1/covLevels,covLevels)))
    #factor(sample(1:covLevels, size=R,replace=T,prob=rep(1/covLevels,covLevels)))
  )
  
  D <- model.matrix(~ .+0, data=t, 
                    contrasts.arg = lapply(t, contrasts, contrasts=FALSE))
  V <- nrow(D)
  
  t$cov1 <- ifelse(t$cov1=='a',1,2) 
  t$rater <- as.integer(as.character(t$rater))
  
  library(Matrix)
  rank <- rankMatrix(D)
  
  phi <- rnorm(R,0,.9)
  psi <- runif(R,0.1,0.7)
  etaB <- as.matrix(c(phi,eta))
  
  rho <- D %*% etaB
  
  xi <- c(Y) # students then items, so length is N*J
  
  tau2 <- 1/(psi^2)
  
  X <- matrix(NA,nrow=N*J,ncol=R)
  
  p <- z <- d <- array(NA,dim=c(N*J,R,K))
  
  for (i in 1:(N*J)) {
    for (r in 1:R) {
      for (k in 1:K) {
        d[i,r,k] <- k - xi[i] - rho[r]
        z[i,r,k] <- exp(-d[i,r,k]*d[i,r,k]*tau2[r]/2)
      }
      for (k in 1:K) {
        p[i,r,k] <- z[i,r,k]/sum(z[i,r,])
      }
      X[i,r] <- sample(1:K,1,prob=p[i,r,])
    }
  }
  
  x <- c(X) # subjects in order, then 
  
  subject <- rep(1:N,J*R)
  rater <- rep(1:R,rep(N*J,R))
  item <- rep(rep(1:J,rep(N,J)),R)
  
  data <- data.frame(x=x,subject=subject,item=item,rater=rater)
  
  data <- data %>% left_join(t, by='rater')
  
  ## less_data (double-scoring)
  
  assign <- as.data.frame.matrix(t(replicate(N, sample(1:R,2))))
  assign$subject <- 1:N
  names(assign)[1:2] <- paste0('rater',1:2)
  
  less_data <- merge(data, assign, by='subject')
  less_data$marker1 <- with(less_data, ifelse(rater == rater1,1,0))
  less_data$marker2 <- with(less_data, ifelse(rater == rater2,1,0))
  #less_data$marker3 <- with(less_data, ifelse(rater == rater3,1,0))
  less_data$keeper <- apply(less_data[,8:9], 1, function(x) sum(x))
  less_data <- less_data[which(less_data$keeper != 0),] 
  less_data <- less_data[,1:5]
  less_data <- less_data[order(less_data$subject, less_data$rater, less_data$item),]
  
  output <- list()
  output$data <- data
  output$less_data <- less_data
  output$N <- N
  output$R <- R
  output$K <- K
  output$J <- J
  output$observedRatings <- X
  output$idealScores <- Y
  output$raterMatrix <- t
  output$rank <- rank
  output$trueValues <- list(
    alpha = alpha,
    beta = beta,
    kappa = kappa,
    betaKappa = gb,
    theta = theta,
    phi = phi,
    psi = psi,
    eta = eta,
    rho = rho
  )
  output
}

rsm.sim <- function(N,J,K,twoPL=F) {
  theta <- rnorm(N)
  alpha <- ifelse(twoPL==T,rlnorm(J,0,.5),1)
  beta <- rnorm(J)
  kappa <- sort(rnorm(K-1,0,1.1))
  
  gb <- array(NA,c(J,K-1))
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      gb[j,k] <- beta[j] + kappa[k]
    }
  }
  
  term <- array(NA, dim=c(N,J,K))
  term[,,1] <- 0
  
  for (n in 1:N) {
    term[n,,2:K] <- sweep(sweep(gb,1,theta[n],"-"),1,-alpha,"*")
  }
  
  Y <- matrix(NA,nrow=N,ncol=J)
  
  for (n in 1:N) {
    for (j in 1:J) {
      term[n,j,] <- exp(cumsum(term[n,j,]))
    }
    term[n,,] <- sweep(term[n,,],1,apply(term[n,,],1,sum),"/")
    
    for (j in 1:J) {
      Y[n,j] <- sample(1:K,1,prob=term[n,j,])
    }
  }
  
  output <- list()
  output$Y <- Y
  output$N <- N
  output$K <- K
  output$J <- J
  output$trueValues <- list(
    alpha = alpha,
    beta = beta,
    kappa = kappa,
    betaKappa = gb,
    theta = theta
  )
  output
}

baseHRM.sim <- function(N,J,K,R,twoPL=F) {
  theta <- rnorm(N)
  alpha <- ifelse(twoPL==T,rlnorm(J,0,.5),1)
  beta <- rnorm(J)
  kappa <- sort(rnorm(K-1,0,1.1))
  
  gb <- array(NA,c(J,K-1))
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      gb[j,k] <- beta[j] + kappa[k]
    }
  }
  
  term <- array(NA, dim=c(N,J,K))
  term[,,1] <- 0
  
  for (n in 1:N) {
    term[n,,2:K] <- sweep(sweep(gb,1,theta[n],"-"),1,-alpha,"*")
  }
  
  Y <- matrix(NA,nrow=N,ncol=J)
  
  for (n in 1:N) {
    for (j in 1:J) {
      term[n,j,] <- exp(cumsum(term[n,j,]))
    }
    term[n,,] <- sweep(term[n,,],1,apply(term[n,,],1,sum),"/")
    
    for (j in 1:J) {
      Y[n,j] <- sample(1:K,1,prob=term[n,j,])
    }
  }
  
  # SDM
  
  phi <- runif(R,0.4,0.9)
  psi <- runif(R,0.1,0.7)
  
  xi <- c(Y) # students then items, so length is N*J
  
  tau2 <- 1/(psi^2)
  
  X <- matrix(NA,nrow=N*J,ncol=R)
  
  p <- z <- d <- array(NA,dim=c(N*J,R,K))
  
  for (i in 1:(N*J)) {
    for (r in 1:R) {
      for (k in 1:K) {
        d[i,r,k] <- k - xi[i] - phi[r]
        z[i,r,k] <- exp(-d[i,r,k]*d[i,r,k]*tau2[r]/2)
      }
      for (k in 1:K) {
        p[i,r,k] <- z[i,r,k]/sum(z[i,r,])
      }
      X[i,r] <- sample(1:K,1,prob=p[i,r,])
    }
  }
  
  x <- c(X) # subjects in order, then 
  
  subject <- rep(1:N,J*R)
  rater <- rep(1:R,rep(N*J,R))
  item <- rep(rep(1:J,rep(N,J)),R)
  
  data <- data.frame(x=x,subject=subject,item=item,rater=rater)
  
  ## less_data (double-scoring)
  
  assign <- as.data.frame.matrix(t(replicate(N, sample(1:R,2))))
  assign$subject <- 1:N
  names(assign)[1:2] <- paste0('rater',1:2)
  
  less_data <- merge(data, assign, by='subject')
  less_data$marker1 <- with(less_data, ifelse(rater == rater1,1,0))
  less_data$marker2 <- with(less_data, ifelse(rater == rater2,1,0))
  #less_data$marker3 <- with(less_data, ifelse(rater == rater3,1,0))
  less_data$keeper <- apply(less_data[,7:8], 1, function(x) sum(x))
  less_data <- less_data[which(less_data$keeper != 0),] 
  less_data <- less_data[,1:4]
  less_data <- less_data[order(less_data$subject, less_data$rater, less_data$item),]
  
  output <- list()
  output$data <- data
  output$less_data <- less_data
  output$N <- N
  output$R <- R
  output$K <- K
  output$J <- J
  output$observedRatings <- X
  output$idealScores <- Y
  output$trueValues <- list(
    alpha = alpha,
    beta = beta,
    kappa = kappa,
    betaKappa = gb,
    theta = theta,
    phi = phi,
    psi = psi
  )
  output
}
