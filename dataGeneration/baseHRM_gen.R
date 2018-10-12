gen_baseHRM_gpcm <- function(N, R, J, K, alpha, gamma, theta, phi, psi) {

  term <- array(NA,dim=c(N,J,K))
  term[,,1] <- 0
  
  for (n in 1:N) {
    term[n,,2:K] <- sweep(sweep(gamma,1,theta[n],"-"),1,-alpha,"*")
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
  
  xi <- c(Y) # students then items, so length is N*J
  
  tau.r <- 1/(psi^2)
  
  X <- matrix(NA,nrow=N*J,ncol=R)
  
  p <- z <- d <- array(NA,dim=c(N*J,R,K))
  
  for (i in 1:(N*J)) {
    for (r in 1:R) {
      for (k in 1:K) {
        d[i,r,k] <- k - xi[i] - phi[r]
        z[i,r,k] <- exp(-d[i,r,k]*d[i,r,k]*tau.r[r]/2)
      }
      for (k in 1:K) {
        p[i,r,k] <- z[i,r,k]/sum(z[i,r,])
      }
      X[i,r] <- sample(1:K,1,prob=p[i,r,])
    }
  }
  
  x <- c(X) # subjects in order, then 
  
  subject <- rep(1:N,J*R)
  rater <-   rep(1:R,rep(N*J,R))
  item <-    rep(rep(1:J,rep(N,J)),R)
  
  data <- data.frame(x=x, subject=subject, item=item, rater=rater)
  
  return(list(data=data,gamma=gamma,alpha=alpha,theta=theta,XI=Y,xi=xi,X=X,phi=phi,psi=psi,tau.r=tau.r))
}

load('~/Desktop/baseHRM_gendata.Rdat')
theta <- gen_baseHRM_data$theta
alpha <- gen_baseHRM_data$alpha
gamma <- gen_baseHRM_data$gamma
phi <- gen_baseHRM_data$phi
psi <- gen_baseHRM_data$psi

gen_baseHRM_data <- gen_baseHRM_gpcm(N = 500, R = 5, J = 8, K = 4, 
                                     alpha = alpha, gamma = gamma, theta = theta,
                                     psi = psi, phi = phi)

save(gen_baseHRM_data, file="../Project5_MIRT/ncme_2017_workshop/Data4Analysis/baseHRM_gendata.Rdat")

dat <- dat %>% spread(key=item, value=x)


res <- immer_HRM(dat = dat[,-c(1:2)], pid = dat$subject, rater = dat$rater,
          iter=3000, burnin=500, est.a = T, est.phi = 'r', est.sigma = F, est.mu = F)

summary(res)
plot(mod1,layout=2,ask=TRUE)

phi1 <- matrix(rep(phi,J),ncol=R,byrow=T)
psi1 <- matrix(rep(psi,J),ncol=R,byrow=T)

dat <- simulate_HRM(theta = theta, a = alpha, b = gamma, phi = phi, psi = psi)
