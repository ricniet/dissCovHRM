
assess <- function(jagsObj, R, filePath, eval, random) {
  library(mcmcplots)
  library(pryr)
  t <- as.mcmc(jagsObj)
  
  # traceplot
  
  trace.rho %<a-% { traplot(
    mcmc = t,
    parms = paste('rho[',1:R,']',sep=''))}
  
  trace.phi %<a-% { traplot(
    mcmc = t,
    parms = paste('phi[',1:R,']',sep=''))}
  
  trace.eta %<a-% { traplot(
    mcmc = t,
    parms = paste('eta[',1:2,']',sep=''))}
  
  trace.etaDelta %<a-% { traplot(
    mcmc = t,
    parms = 'etaDelta')}
  
  if (random==T) {
    trace.sdRho %<a-% { traplot(
      mcmc = t,
      parms = 'sd.rho')}
  }
  
  # density
  
  density.rho %<a-% { denplot(
    mcmc = t,
    parms = paste('rho[',1:R,']',sep=''))}
  
  density.phi %<a-% { denplot(
    mcmc = t,
    parms = paste('phi[',1:R,']',sep=''))}
  
  density.eta %<a-% { denplot(
    mcmc = t,
    parms = paste('eta[',1:2,']',sep=''))}
  
  density.etaDelta %<a-% { denplot(
    mcmc = t,
    parms = 'etaDelta')}
  
  if (random==T) {
    density.sdRho %<a-% { traplot(
      mcmc = t,
      parms = 'sd.rho')}
  }
  
  # save traceplots
  
  png(paste(filePath,"resultsPlots/trace_rho.png",sep=""),width = 20, height = 20, units = 'in', res = 90)
  trace.rho
  dev.off()
  
  png(paste(filePath,"resultsPlots/trace_phi.png",sep=""),width = 20, height = 20, units = 'in', res = 90)
  trace.phi
  dev.off()
  
  png(paste(filePath,"resultsPlots/trace_eta.png",sep=""))
  trace.eta
  dev.off()
  
  png(paste(filePath,"resultsPlots/trace_etaDelta.png",sep=""))
  trace.etaDelta
  dev.off()
  
  if (random==T) {
    png(paste(filePath,"resultsPlots/trace_sdRho.png",sep=""))
    trace.sdRho
    dev.off()
  }
  
  # save density plots
  
  png(paste(filePath,"resultsPlots/density_rho.png",sep=""),width = 20, height = 20, units = 'in', res = 90)
  density.rho
  dev.off()
  
  png(paste(filePath,"resultsPlots/density_phi.png",sep=""),width = 20, height = 20, units = 'in', res = 90)
  density.phi
  dev.off()
  
  png(paste(filePath,"resultsPlots/density_eta.png",sep=""))
  density.eta
  dev.off()
  
  png(paste(filePath,"resultsPlots/density_etaDelta.png",sep=""))
  density.etaDelta
  dev.off()
  
  if (random==T) {
    png(paste(filePath,"resultsPlots/density_sdRho.png",sep=""))
    density.sdRho
    dev.off()
  }
  
  # calculate bias for eval

  biasRho <- covHRM_results[which(substr(row.names(covHRM_results),1,3)=='rho'),1][eval] - simDat$trueValues$rho[eval]
  aveBiasRho <- round(mean(covHRM_results[which(substr(row.names(covHRM_results),1,3)=='rho'),1][eval] - simDat$trueValues$rho[eval]),3)
  sdBiasRho <- round(sd(covHRM_results[which(substr(row.names(covHRM_results),1,3)=='rho'),1][eval] - simDat$trueValues$rho[eval]),3)
  
  # if (random==T) {
  #   biasSdRho <- round(mean(covHRM_results[which(substr(row.names(covHRM_results),1,3)=='sd.'),1][eval] - simDat$trueValues$rho[eval]),3)
  # }
  
  biasPhi <- covHRM_results[which(substr(row.names(covHRM_results),1,3)=='phi'),1][eval] - simDat$trueValues$phi[eval]
  aveBiasPhi <- round(mean(covHRM_results[which(substr(row.names(covHRM_results),1,3)=='phi'),1][eval] - simDat$trueValues$phi[eval]),3)
  sdBiasPhi <- round(sd(covHRM_results[which(substr(row.names(covHRM_results),1,3)=='phi'),1][eval] - simDat$trueValues$phi[eval]),3)
  
  biasEta1 <- round(covHRM_results[which(row.names(covHRM_results)=='eta[1]'),1] - simDat$trueValues$eta[1],3)
  biasEta2 <- round(covHRM_results[which(row.names(covHRM_results)=='eta[2]'),1] - simDat$trueValues$eta[2],3)
  
  biasEtaDelta <- round(covHRM_results[which(row.names(covHRM_results)=='etaDelta'),1] - (simDat$trueValues$eta[1] - simDat$trueValues$eta[2]),3)
  
  cat("RHO, Average bias (SD):\n",sprintf("%s (%s)",aveBiasRho,sdBiasRho),"\n\n",
    "PHI, Average bias (SD):\n",sprintf("%s (%s)",aveBiasPhi,sdBiasPhi),"\n\n",
    "ETA1, Bias:\n",sprintf("%s",biasEta1),"\n\n",
    "ETA2, Bias:\n",sprintf("%s",biasEta2),"\n\n",
    "etaDELTA, Bias:\n",sprintf("%s",biasEtaDelta),"\n\n",
    file=paste(filePath,"biasResults.txt",sep="")
  )
}

assess_etaOnly <- function(jagsObj, R, filePath, eval, random) {
  library(mcmcplots)
  library(pryr)
  t <- as.mcmc(jagsObj)
  
  # traceplot
  
  trace.eta %<a-% { traplot(
    mcmc = t,
    parms = paste('eta[',1:2,']',sep=''))}
  
  trace.etaDelta %<a-% { traplot(
    mcmc = t,
    parms = 'etaDelta')}
  
  if (random==T) {
    trace.sdRho %<a-% { traplot(
      mcmc = t,
      parms = 'sd.rho')}
  }
  
  # density
  
  density.eta %<a-% { denplot(
    mcmc = t,
    parms = paste('eta[',1:2,']',sep=''))}
  
  density.etaDelta %<a-% { denplot(
    mcmc = t,
    parms = 'etaDelta')}
  
  if (random==T) {
    density.sdRho %<a-% { traplot(
      mcmc = t,
      parms = 'sd.rho')}
  }
  
  # save traceplots
  
  png(paste(filePath,"resultsPlots/trace_eta.png",sep=""))
  trace.eta
  dev.off()
  
  png(paste(filePath,"resultsPlots/trace_etaDelta.png",sep=""))
  trace.etaDelta
  dev.off()
  
  if (random==T) {
    png(paste(filePath,"resultsPlots/trace_sdRho.png",sep=""))
    trace.sdRho
    dev.off()
  }
  
  # save density plots
  
  png(paste(filePath,"resultsPlots/density_eta.png",sep=""))
  density.eta
  dev.off()
  
  png(paste(filePath,"resultsPlots/density_etaDelta.png",sep=""))
  density.etaDelta
  dev.off()
  
  if (random==T) {
    png(paste(filePath,"resultsPlots/density_sdRho.png",sep=""))
    density.sdRho
    dev.off()
  }
  
  # calculate bias for eval
  
  biasEta1 <- round(covHRM_results[which(row.names(covHRM_results)=='eta[1]'),1] - simDat$trueValues$eta[1],3)
  biasEta2 <- round(covHRM_results[which(row.names(covHRM_results)=='eta[2]'),1] - simDat$trueValues$eta[2],3)
  
  biasEtaDelta <- round(covHRM_results[which(row.names(covHRM_results)=='etaDelta'),1] - (simDat$trueValues$eta[1] - simDat$trueValues$eta[2]),3)
  
  cat("ETA1, Bias:\n",sprintf("%s",biasEta1),"\n\n",
      "ETA2, Bias:\n",sprintf("%s",biasEta2),"\n\n",
      "etaDELTA, Bias:\n",sprintf("%s",biasEtaDelta),"\n\n",
      file=paste(filePath,"biasResults.txt",sep="")
  )
}
