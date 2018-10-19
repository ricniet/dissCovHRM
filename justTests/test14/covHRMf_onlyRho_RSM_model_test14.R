library(R2jags); library(tictoc); library(tidyverse)

source('./dataGeneration/genModels.R')

load('./dataGeneration/generatedDataSets/forTests/set5_simDat_fcovHRM_onlyRHO_N5000R200J5K5.Rdat')
#simDat <- covHRMf.onlyRho.sim(N=500,J=8,K=5,R=40,twoPL=F,eta=c(.88,1.1))
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
#for saving rho
for (i in 1:length(rMatrix[which(rMatrix$cov1==1),1])) {
  cat('rho[',rMatrix[which(rMatrix$cov1==1),1][i],'] <- ',
      'phi[',rMatrix[which(rMatrix$cov1==1),1][i],'] + eta[1]',
      sep='',collapse='\n')
}
for (i in 1:length(rMatrix[which(rMatrix$cov1==2),1])) {
  cat('rho[',rMatrix[which(rMatrix$cov1==2),1][i],'] <- ',
      'phi[',rMatrix[which(rMatrix$cov1==2),1][i],'] + eta[2]',
      sep='',collapse='\n')
}
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
  phi[200] <- -(phi[1]+phi[2]+phi[3]+phi[4]+phi[9]+phi[10]+phi[11]+phi[16]+
                  phi[20]+phi[21]+phi[22]+phi[23]+phi[25]+phi[26]+phi[29]+
                  phi[33]+phi[35]+phi[36]+phi[39]+phi[40]+phi[43]+phi[44]+
                  phi[47]+phi[48]+phi[50]+phi[51]+phi[53]+phi[57]+phi[58]+
                  phi[59]+phi[60]+phi[62]+phi[63]+phi[66]+phi[69]+phi[70]+
                  phi[72]+phi[73]+phi[74]+phi[76]+phi[79]+phi[81]+phi[82]+
                  phi[84]+phi[86]+phi[88]+phi[89]+phi[90]+phi[91]+phi[92]+
                  phi[93]+phi[95]+phi[96]+phi[97]+phi[98]+phi[100]+phi[102]+
                  phi[103]+phi[104]+phi[108]+phi[110]+phi[111]+phi[112]+
                  phi[115]+phi[116]+phi[117]+phi[119]+phi[120]+phi[122]+
                  phi[123]+phi[124]+phi[125]+phi[126]+phi[127]+phi[128]+
                  phi[129]+phi[133]+phi[136]+phi[137]+phi[140]+phi[141]+
                  phi[144]+phi[146]+phi[148]+phi[151]+phi[152]+phi[154]+
                  phi[155]+phi[157]+phi[158]+phi[159]+phi[160]+phi[164]+
                  phi[165]+phi[170]+phi[171]+phi[172]+phi[173]+phi[176]+
                  phi[180]+phi[183]+phi[186]+phi[192]+phi[194]+phi[196]+
                  phi[197]+phi[199])
  
  phi[198] <- -(phi[5]+phi[6]+phi[7]+phi[8]+phi[12]+phi[13]+phi[14]+phi[15]+
                  phi[17]+phi[18]+phi[19]+phi[24]+phi[27]+phi[28]+phi[30]+
                  phi[31]+phi[32]+phi[34]+phi[37]+phi[38]+phi[41]+phi[42]+
                  phi[45]+phi[46]+phi[49]+phi[52]+phi[54]+phi[55]+phi[56]+
                  phi[61]+phi[64]+phi[65]+phi[67]+phi[68]+phi[71]+phi[75]+
                  phi[77]+phi[78]+phi[80]+phi[83]+phi[85]+phi[87]+phi[94]+
                  phi[99]+phi[101]+phi[105]+phi[106]+phi[107]+phi[109]+
                  phi[113]+phi[114]+phi[118]+phi[121]+phi[130]+phi[131]+
                  phi[132]+phi[134]+phi[135]+phi[138]+phi[139]+phi[142]+
                  phi[143]+phi[145]+phi[147]+phi[149]+phi[150]+phi[153]+
                  phi[156]+phi[161]+phi[162]+phi[163]+phi[166]+phi[167]+
                  phi[168]+phi[169]+phi[174]+phi[175]+phi[177]+phi[178]+
                  phi[179]+phi[181]+phi[182]+phi[184]+phi[185]+phi[187]+
                  phi[188]+phi[189]+phi[190]+phi[191]+phi[193]+phi[195])
  
  # CHANGE CHANGE:
  for (nu in 1:(197)){      
    phi[nu] ~ dnorm(0, prec.phi)
  }
  
  phi[199] ~ dnorm(0, prec.phi)
  
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
  rho[2] <- phi[2] + eta[1]
  rho[3] <- phi[3] + eta[1]
  rho[4] <- phi[4] + eta[1]
  rho[9] <- phi[9] + eta[1]
  rho[10] <- phi[10] + eta[1]
  rho[11] <- phi[11] + eta[1]
  rho[16] <- phi[16] + eta[1]
  rho[20] <- phi[20] + eta[1]
  rho[21] <- phi[21] + eta[1]
  rho[22] <- phi[22] + eta[1]
  rho[23] <- phi[23] + eta[1]
  rho[25] <- phi[25] + eta[1]
  rho[26] <- phi[26] + eta[1]
  rho[29] <- phi[29] + eta[1]
  rho[33] <- phi[33] + eta[1]
  rho[35] <- phi[35] + eta[1]
  rho[36] <- phi[36] + eta[1]
  rho[39] <- phi[39] + eta[1]
  rho[40] <- phi[40] + eta[1]
  rho[43] <- phi[43] + eta[1]
  rho[44] <- phi[44] + eta[1]
  rho[47] <- phi[47] + eta[1]
  rho[48] <- phi[48] + eta[1]
  rho[50] <- phi[50] + eta[1]
  rho[51] <- phi[51] + eta[1]
  rho[53] <- phi[53] + eta[1]
  rho[57] <- phi[57] + eta[1]
  rho[58] <- phi[58] + eta[1]
  rho[59] <- phi[59] + eta[1]
  rho[60] <- phi[60] + eta[1]
  rho[62] <- phi[62] + eta[1]
  rho[63] <- phi[63] + eta[1]
  rho[66] <- phi[66] + eta[1]
  rho[69] <- phi[69] + eta[1]
  rho[70] <- phi[70] + eta[1]
  rho[72] <- phi[72] + eta[1]
  rho[73] <- phi[73] + eta[1]
  rho[74] <- phi[74] + eta[1]
  rho[76] <- phi[76] + eta[1]
  rho[79] <- phi[79] + eta[1]
  rho[81] <- phi[81] + eta[1]
  rho[82] <- phi[82] + eta[1]
  rho[84] <- phi[84] + eta[1]
  rho[86] <- phi[86] + eta[1]
  rho[88] <- phi[88] + eta[1]
  rho[89] <- phi[89] + eta[1]
  rho[90] <- phi[90] + eta[1]
  rho[91] <- phi[91] + eta[1]
  rho[92] <- phi[92] + eta[1]
  rho[93] <- phi[93] + eta[1]
  rho[95] <- phi[95] + eta[1]
  rho[96] <- phi[96] + eta[1]
  rho[97] <- phi[97] + eta[1]
  rho[98] <- phi[98] + eta[1]
  rho[100] <- phi[100] + eta[1]
  rho[102] <- phi[102] + eta[1]
  rho[103] <- phi[103] + eta[1]
  rho[104] <- phi[104] + eta[1]
  rho[108] <- phi[108] + eta[1]
  rho[110] <- phi[110] + eta[1]
  rho[111] <- phi[111] + eta[1]
  rho[112] <- phi[112] + eta[1]
  rho[115] <- phi[115] + eta[1]
  rho[116] <- phi[116] + eta[1]
  rho[117] <- phi[117] + eta[1]
  rho[119] <- phi[119] + eta[1]
  rho[120] <- phi[120] + eta[1]
  rho[122] <- phi[122] + eta[1]
  rho[123] <- phi[123] + eta[1]
  rho[124] <- phi[124] + eta[1]
  rho[125] <- phi[125] + eta[1]
  rho[126] <- phi[126] + eta[1]
  rho[127] <- phi[127] + eta[1]
  rho[128] <- phi[128] + eta[1]
  rho[129] <- phi[129] + eta[1]
  rho[133] <- phi[133] + eta[1]
  rho[136] <- phi[136] + eta[1]
  rho[137] <- phi[137] + eta[1]
  rho[140] <- phi[140] + eta[1]
  rho[141] <- phi[141] + eta[1]
  rho[144] <- phi[144] + eta[1]
  rho[146] <- phi[146] + eta[1]
  rho[148] <- phi[148] + eta[1]
  rho[151] <- phi[151] + eta[1]
  rho[152] <- phi[152] + eta[1]
  rho[154] <- phi[154] + eta[1]
  rho[155] <- phi[155] + eta[1]
  rho[157] <- phi[157] + eta[1]
  rho[158] <- phi[158] + eta[1]
  rho[159] <- phi[159] + eta[1]
  rho[160] <- phi[160] + eta[1]
  rho[164] <- phi[164] + eta[1]
  rho[165] <- phi[165] + eta[1]
  rho[170] <- phi[170] + eta[1]
  rho[171] <- phi[171] + eta[1]
  rho[172] <- phi[172] + eta[1]
  rho[173] <- phi[173] + eta[1]
  rho[176] <- phi[176] + eta[1]
  rho[180] <- phi[180] + eta[1]
  rho[183] <- phi[183] + eta[1]
  rho[186] <- phi[186] + eta[1]
  rho[192] <- phi[192] + eta[1]
  rho[194] <- phi[194] + eta[1]
  rho[196] <- phi[196] + eta[1]
  rho[197] <- phi[197] + eta[1]
  rho[199] <- phi[199] + eta[1]
  rho[200] <- phi[200] + eta[1]
  
  rho[5] <- phi[5] + eta[2]
  rho[6] <- phi[6] + eta[2]
  rho[7] <- phi[7] + eta[2]
  rho[8] <- phi[8] + eta[2]
  rho[12] <- phi[12] + eta[2]
  rho[13] <- phi[13] + eta[2]
  rho[14] <- phi[14] + eta[2]
  rho[15] <- phi[15] + eta[2]
  rho[17] <- phi[17] + eta[2]
  rho[18] <- phi[18] + eta[2]
  rho[19] <- phi[19] + eta[2]
  rho[24] <- phi[24] + eta[2]
  rho[27] <- phi[27] + eta[2]
  rho[28] <- phi[28] + eta[2]
  rho[30] <- phi[30] + eta[2]
  rho[31] <- phi[31] + eta[2]
  rho[32] <- phi[32] + eta[2]
  rho[34] <- phi[34] + eta[2]
  rho[37] <- phi[37] + eta[2]
  rho[38] <- phi[38] + eta[2]
  rho[41] <- phi[41] + eta[2]
  rho[42] <- phi[42] + eta[2]
  rho[45] <- phi[45] + eta[2]
  rho[46] <- phi[46] + eta[2]
  rho[49] <- phi[49] + eta[2]
  rho[52] <- phi[52] + eta[2]
  rho[54] <- phi[54] + eta[2]
  rho[55] <- phi[55] + eta[2]
  rho[56] <- phi[56] + eta[2]
  rho[61] <- phi[61] + eta[2]
  rho[64] <- phi[64] + eta[2]
  rho[65] <- phi[65] + eta[2]
  rho[67] <- phi[67] + eta[2]
  rho[68] <- phi[68] + eta[2]
  rho[71] <- phi[71] + eta[2]
  rho[75] <- phi[75] + eta[2]
  rho[77] <- phi[77] + eta[2]
  rho[78] <- phi[78] + eta[2]
  rho[80] <- phi[80] + eta[2]
  rho[83] <- phi[83] + eta[2]
  rho[85] <- phi[85] + eta[2]
  rho[87] <- phi[87] + eta[2]
  rho[94] <- phi[94] + eta[2]
  rho[99] <- phi[99] + eta[2]
  rho[101] <- phi[101] + eta[2]
  rho[105] <- phi[105] + eta[2]
  rho[106] <- phi[106] + eta[2]
  rho[107] <- phi[107] + eta[2]
  rho[109] <- phi[109] + eta[2]
  rho[113] <- phi[113] + eta[2]
  rho[114] <- phi[114] + eta[2]
  rho[118] <- phi[118] + eta[2]
  rho[121] <- phi[121] + eta[2]
  rho[130] <- phi[130] + eta[2]
  rho[131] <- phi[131] + eta[2]
  rho[132] <- phi[132] + eta[2]
  rho[134] <- phi[134] + eta[2]
  rho[135] <- phi[135] + eta[2]
  rho[138] <- phi[138] + eta[2]
  rho[139] <- phi[139] + eta[2]
  rho[142] <- phi[142] + eta[2]
  rho[143] <- phi[143] + eta[2]
  rho[145] <- phi[145] + eta[2]
  rho[147] <- phi[147] + eta[2]
  rho[149] <- phi[149] + eta[2]
  rho[150] <- phi[150] + eta[2]
  rho[153] <- phi[153] + eta[2]
  rho[156] <- phi[156] + eta[2]
  rho[161] <- phi[161] + eta[2]
  rho[162] <- phi[162] + eta[2]
  rho[163] <- phi[163] + eta[2]
  rho[166] <- phi[166] + eta[2]
  rho[167] <- phi[167] + eta[2]
  rho[168] <- phi[168] + eta[2]
  rho[169] <- phi[169] + eta[2]
  rho[174] <- phi[174] + eta[2]
  rho[175] <- phi[175] + eta[2]
  rho[177] <- phi[177] + eta[2]
  rho[178] <- phi[178] + eta[2]
  rho[179] <- phi[179] + eta[2]
  rho[181] <- phi[181] + eta[2]
  rho[182] <- phi[182] + eta[2]
  rho[184] <- phi[184] + eta[2]
  rho[185] <- phi[185] + eta[2]
  rho[187] <- phi[187] + eta[2]
  rho[188] <- phi[188] + eta[2]
  rho[189] <- phi[189] + eta[2]
  rho[190] <- phi[190] + eta[2]
  rho[191] <- phi[191] + eta[2]
  rho[193] <- phi[193] + eta[2]
  rho[195] <- phi[195] + eta[2]
  rho[198] <- phi[198] + eta[2]
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
    phi = c(rnorm(197,0,1),NA,rnorm(1,0,1),NA),
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

write.csv(covHRM_results, file='./justTests/test14/covHRMf_RHOonly_test14_N5000R200J5K5_3ch20000i10000b20t.csv')
save(covHRM_jags, file='./justTests/test14/covHRMf_RHOonly_test14_N5000R200J5K5_3ch20000i10000b20t.Rdat')

source('./justTests/assessTests.R')

assess(jagsObj = covHRM_jags, 
       R = R,
       filePath = './justTests/test14/',
       eval = c(1:197,199),
       random = F)
