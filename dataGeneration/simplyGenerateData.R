library(tidyverse)

source('./dataGeneration/genModels.R')

# fixed model

#load('./testResults/fcovHRM_onlyRHO/simDat_fcovHRM_onlyRHO_N500R40J8K5.Rdat')
simDat <- covHRMf.onlyRho.sim(N=500,J=8,K=5,R=50,twoPL=F,eta=c(.88,1.1))
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

save(simDat, file='./dataGeneration/generatedDataSets/forTests/set2_simDat_fcovHRM_onlyRHO_N500R50J8K5.Rdat')

# random model

#load('./testResults/fcovHRM_onlyRHO/simDat_fcovHRM_onlyRHO_N500R40J8K5.Rdat')
simDat <- covHRMr.onlyRho.sim(N=500,J=8,K=5,R=50,twoPL=F,eta=c(.88,1.1),sd.eta=0.8)
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

save(simDat, file='./dataGeneration/generatedDataSets/forTests/set3_simDat_rcovHRM_onlyRHO_N500R50J8K5.Rdat')
