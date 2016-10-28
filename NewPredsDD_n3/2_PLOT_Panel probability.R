#####################################################################
###############################################################
##### Plot Popl'n average effect


########## SCATTER PLOT 


# Population-average parameters
ests <- out1$BUGSoutput$summary[c("mu.alpha","mu.beta"),1]

# Fake predictor
X3 <- seq(min(dat$x), max(dat$x), length=100)


# Linear predictor for population average effect
lpredPopAve <- ests[1] + ests[2]*X3 
plot(plogis(lpredPopAve)~X3)

# Number of simulations
N <- out1$BUGSoutput$n.sims

linPredPopAve <- matrix(NA, ncol=length(X3), nrow=N ) #container for predicted values
dim(linPredPopAve)

mu.A <- out1$BUGSoutput$sims.list$mu.alpha
mu.B <- out1$BUGSoutput$sims.list$mu.beta

for(i in 1:N){
  for(t in 1:length(X3)){
    linPredPopAve[i,t] <- mu.A[i] + mu.B[i]*X3[t] 
  }
}

dim(linPredPopAve)

probPredPopAve <- plogis(linPredPopAve)

meanProbPopAve <- apply(probPredPopAve, 2, mean)
upperCI.PopAve <- apply(probPredPopAve, 2, quantile, probs=c(0.95) )
lowerCIA.PopAve <- apply(probPredPopAve, 2, quantile, probs=c(0.05) )

