# rm(list=ls())
library(R2jags)
library(lme4)
library(MCMCpack)
library(arm)
library(lattice)
library(PerformanceAnalytics)
library(MuMIn)

# Read in data
dat <- read.csv('walleye_recruitment_new_predictors_correct_mean_GDD.csv')
dim(dat)
head(dat)

# Read in bass data
bassDat <- read.csv('fall_lmb8_cpue_lakemean.csv')

dim(bassDat)
length(unique(bassDat$WBIC))

dat <- dat[dat$WBIC %in% bassDat$WBIC,]
dim(dat)

# Remove lakes with missing data
dat <- dat[complete.cases(dat),]
length(unique(dat$WBIC))

length(unique(dat$WBIC))
head(dat)

# Rename and scale covarites
dat$y <- dat$recruitment
# Log-transform and then standardize level-1 predictors
dat$gdd <- log(dat$GDD_wtr_5c) 
dat$gddz <- as.numeric(scale(dat$gdd))
dat$x <- dat$gddz

dat$icez <- as.numeric(scale(dat$ice_off_julian_date))

dat$spring <-log(dat$coef_var_30.60)
dat$springz <- as.numeric(scale(dat$spring))

cor(dat$gddz, dat$icez)

dat$lakenum <- as.numeric(as.factor(as.numeric(dat$WBIC)))
summary(dat)

# renumber lakenum
dat$lakenum <- as.numeric(as.factor(as.numeric(dat$lakenum)))

samp.sizes <- as.numeric(table(dat$lakenum))
mean(samp.sizes)
median(samp.sizes)
range(samp.sizes)

mean(dat$GDD_wtr_5c)
range(dat$GDD_wtr_5c)
median(dat$GDD_wtr_5c)

#################################################################
########## BUGS CODE ############################################
#################################################################

# Define the model in the BUGS language and write a text file
sink("model.txt")
cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dbin(p[i],1)				# distributional assumption
    p[i] <- exp(lp[i])/(1+exp(lp[i])) # logit link function
    lp[i] <- alpha[group[i]] + beta[group[i]] * x[i]	# linear predictor    
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta[j] <- BB[j,2]
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivriate normal
    
    BB.hat[j,1] <- mu.alpha + gamma0b * z1[j] + gamma0b2 * z2[j] + gamma0b3 * z3[j] + gamma0b4 * z4[j]
    BB.hat[j,2] <- mu.beta + gamma1b * z1[j] + gamma1b2 * z2[j] + gamma1b3 * z3[j] + gamma1b4 * z4[j]
    }
    
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta ~ dnorm(0, 0.0001)
    gamma0b ~ dnorm(0, 0.0001)
    gamma0b2 ~ dnorm(0, 0.0001)
    gamma0b3 ~ dnorm(0, 0.0001)
    gamma0b4 ~ dnorm(0, 0.0001)
    gamma1b ~ dnorm(0, 0.0001)
    gamma1b2 ~ dnorm(0, 0.0001)
    gamma1b3 ~ dnorm(0, 0.0001)
    gamma1b4 ~ dnorm(0, 0.0001)
    
    
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()



# Number of parameters
K <- 2

# Create identity matrix for Wishart dist'n
#!!!!!!!Number of parameters to estimate (K)

W <- diag(K)

# Level-2 covariates
# lake area
area <- as.numeric(by(dat$area.hectares, dat$lakenum, mean)) 
area <- as.numeric(scale(log(area)))
hist(area)

area2 <- as.numeric(by(dat$area.hectares, dat$lakenum, mean)) 
range(area2)
mean(area2)
median(area2)

# Conductivity
cond <- as.numeric(by(dat$Conductance, dat$lakenum, mean)) 
cond <- as.numeric(scale(log(cond)))
hist(cond)

cond2 <- as.numeric(by(dat$Conductance, dat$lakenum, mean))
range(cond2)
mean(cond2)
median(cond2)
# sd(dat$Conductance)
# Bass CPUE
bassDat2 <- bassDat[bassDat$WBIC %in% dat$WBIC,]
dim(bassDat2)
hist(bassDat2$mean.cpue.fall)
bass <- as.numeric(scale(log(bassDat2$mean.cpue.fall+0.5))) # add a constant so we can log-transform
hist(bass)

range(bassDat2$mean.cpue.fall)
mean(bassDat2$mean.cpue.fall)
median(bassDat2$mean.cpue.fall)

# Latitude
lat <- as.numeric(by(dat$lat, dat$lakenum, mean)) 
hist(lat)

# GDD
gdd <- as.numeric(by(dat$GDD_wtr_5c, dat$lakenum, mean)) 
gdd2 <- as.numeric(scale(log(gdd)))
mat1 <- cbind(area, cond, bass, gdd2) # gdd and lat correlated: r = -0.75
cor(mat1)

gdd.raw <- as.numeric(by(dat$GDD_wtr_5c, dat$lakenum, mean)) 
range(gdd.raw)
mean(gdd.raw)
median(gdd.raw)


# Number of lakes
J <- length(unique(dat$lakenum))


# Load data
data <- list(y = dat$y, group = dat$lakenum, n = dim(dat)[1], J = J,
             x=dat$gddz, K=K, W=W, z1 = area, z2=cond, z3=bass, z4=gdd2 )


# Get data in form for lmer
dat$cond <- cond[dat$lakenum]
dat$area <- area[dat$lakenum]
dat$bass <- bass[dat$lakenum]
dat$gdd2 <- gdd2[dat$lakenum]


# Get marginal (fixed effects only) and conditional (both fixed and random effects) r2 for logistic model
# Nakagawa, S, Schielzeth, H. (2013). A general and simple method for obtaining R² from Generalized Linear Mixed-effects Models. Methods in Ecology and Evolution 4: 133–142
m1 <- glmer(y ~ 1 + x + area + x:area + cond + x:cond + bass + x:bass + (1+x|lakenum), 
            control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),
            family = binomial, data=dat)
summary(m1)
# R2 (use MuMln package - )
r2s <- r.squaredGLMM(m1)
r2s

# m2 <- glmer(y ~ 1 + x + area + x:area + cond + x:cond + bass + x:bass + gdd2 + x:gdd2 + (1+x|lakenum),
#             control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),
#             family = binomial, data=dat)
# summary(m2)
# # R2 (use MuMln package - )
# r22s <- r.squaredGLMM(m2)
# r22s

# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1), mu.beta=rnorm(1), 
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)),
        gamma0b = rnorm(1), gamma1b=rnorm(1),gamma0b2=rnorm(1),gamma0b3=rnorm(1),gamma0b4=rnorm(1),
        gamma1b2=rnorm(1),gamma1b3=rnorm(1),gamma1b4=rnorm(1) )
}


# Parameters monitored
parameters <- c("mu.alpha","mu.beta","BB", "Sigma.B","gamma0b","gamma1b","gamma0b2","gamma0b3","gamma0b4","gamma1b2","gamma1b3","gamma1b4")


# MCMC settings
ni <- 60000
nt <- 2
nb <- 30000
nc <- 3


start.time = Sys.time()         # Set timer 
# Call JAGS from R 

out1 <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize posteriors
print(out1, dig = 3)

# Find which parameters, if any, have Rhat > 1.1
which(out1$BUGSoutput$summary[, c("Rhat")] > 1.1)

### Calculate the probability that the effect is the direction of the estimated posterior mean
slope.means <- out1$BUGSoutput$mean$BB[,2]
probs1 <- numeric()
for(i in 1:length(slope.means)){
if( slope.means[i] > 0){
  probs1[i] <- mean(out1$BUGSoutput$sims.list$BB[,i,2] > 0)
} else {
  probs1[i] <- mean(out1$BUGSoutput$sims.list$BB[,i,2] < 0)
 }	
} # close loop


probsPos <- numeric()
probsNeg <- numeric()
for(i in 1:length(slope.means)){
  if( slope.means[i] > 0){
    probsPos[i] <- mean(out1$BUGSoutput$sims.list$BB[,i,2] > 0)
  } else {
    probsNeg[i] <- mean(out1$BUGSoutput$sims.list$BB[,i,2] < 0)
  }	
} # close loop

length(probsPos)
length(probsNeg)



probsPos <- probsPos[!is.na(probsPos)]
length(probsPos)
range(probsPos)
mean(probsPos)
median(probsPos)
sum(probsPos > 0.7)


probsNegInd <- probsNeg
# Indicator for which lakes have a negative effect of DD
probsNegInd2 <- as.numeric(!is.na(probsNegInd))
length(probsNegInd2)
sum(probsNegInd2)

probsNeg <- probsNeg[!is.na(probsNeg)]
length(probsNeg)
range(probsNeg)
mean(probsNeg)
median(probsNeg)
sum(probsNeg > 0.7)



# probs1

betas <- out1$BUGSoutput$mean$BB[,2]
betasCIsL <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.05)) 
betasCIsU <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.95)) 
sigBetas <- betasCIsL * betasCIsU > 0
sum(sigBetas) # was 25, now 6

betasCIsL80 <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.10)) 
betasCIsU80 <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.90)) 
sigBetas80 <- betasCIsL80 * betasCIsU80 > 0
sum(sigBetas80) # was 93, now 55

hist(betas)
abline(v=0)

# Number of negative and positve slopes
length(which(betas > 0)) # 135; 161
length(which(betas < 0)) # 229; 198

which(betas > 0)
which(betas > 0.95)
which(betas < -1.4)
betas[186]
sigBetas[186]

# covariate order: area, cond, bass, GDD
mean(out1$BUGSoutput$sims.list$mu.alpha) 
mean(out1$BUGSoutput$sims.list$gamma0b) # sig (+)
mean(out1$BUGSoutput$sims.list$gamma0b2) # NS (+)
mean(out1$BUGSoutput$sims.list$gamma0b3) # Sig (-)
mean(out1$BUGSoutput$sims.list$gamma0b4) # NS (-)

quantile(out1$BUGSoutput$sims.list$mu.alpha, c(0.05, 0.95)) 
quantile(out1$BUGSoutput$sims.list$gamma0b, c(0.05, 0.95)) # sig (+)
quantile(out1$BUGSoutput$sims.list$gamma0b2, c(0.05, 0.95)) # NS
quantile(out1$BUGSoutput$sims.list$gamma0b3, c(0.05, 0.95)) # Sig (-)
quantile(out1$BUGSoutput$sims.list$gamma0b4, c(0.05, 0.95)) # NS
# mean(out1$BUGSoutput$sims.list$gamma0b3) 

mean(out1$BUGSoutput$sims.list$mu.beta)
mean(out1$BUGSoutput$sims.list$gamma1b) # NS (-)
mean(out1$BUGSoutput$sims.list$gamma1b2) # Sig (-) - NS 
mean(out1$BUGSoutput$sims.list$gamma1b3) # Sig (-)
mean(out1$BUGSoutput$sims.list$gamma1b4) # NS (-)

quantile(out1$BUGSoutput$sims.list$mu.beta, c(0.05, 0.95)) 
# quantile(out1$BUGSoutput$sims.list$mu.beta, c(0.1, 0.9)) # 80% CI
quantile(out1$BUGSoutput$sims.list$gamma1b, c(0.05, 0.95)) # NS
quantile(out1$BUGSoutput$sims.list$gamma1b2, c(0.05, 0.975)) # Sig (-) - NS 
quantile(out1$BUGSoutput$sims.list$gamma1b3, c(0.05, 0.95)) # Sig (-)
# quantile(out1$BUGSoutput$sims.list$gamma1b3, c(0.1, 0.9)) # Sig (-)
# quantile(out1$BUGSoutput$sims.list$gamma1b3, c(0.05, 0.95)) # Sig at 90% level
quantile(out1$BUGSoutput$sims.list$gamma1b4, c(0.05, 0.95)) # NS
# mean(out1$BUGSoutput$sims.list$gamma1b3) 


# Probability effects are in direction of mean
mean(out1$BUGSoutput$sims.list$mu.alpha < 0) 
mean(out1$BUGSoutput$sims.list$gamma0b > 0) # sig (+)
mean(out1$BUGSoutput$sims.list$gamma0b2 > 0) # NS (+)
mean(out1$BUGSoutput$sims.list$gamma0b3 < 0) # Sig (-)
mean(out1$BUGSoutput$sims.list$gamma0b4 < 0) # NS (-)


mean(out1$BUGSoutput$sims.list$mu.beta < 0)
mean(out1$BUGSoutput$sims.list$gamma1b < 0)
mean(out1$BUGSoutput$sims.list$gamma1b2 < 0)
mean(out1$BUGSoutput$sims.list$gamma1b3 < 0)
mean(out1$BUGSoutput$sims.list$gamma1b4 < 0)

# Calculate correlation between varying slopes and intercepts
sd.a <- out1$BUGSoutput$mean$Sigma.B[1,1]
sd.b <- out1$BUGSoutput$mean$Sigma.B[2,2]
cov.ab <- out1$BUGSoutput$mean$Sigma.B[1,2]

corr.ab <- cov.ab/(sd.a * sd.b)
corr.ab

summary <- out1$BUGSoutput$summary
write.csv(summary,'final_model_summary.csv',row.names = T)


# Sort dat by lakenum (or WBIC)
dat2 <- dat[order(dat$lakenum),]
head(dat2)

wbics <- unique(dat2$WBIC)
write.csv(wbics,'final_wbics2.csv',row.names = F)

ProbsForG <- data.frame(wbics, probs1, probsNegInd2)
colnames(ProbsForG)<- c('WBIC','probability','Negative_Indicator')
dim(ProbsForG)
head(ProbsForG)

write.csv(ProbsForG, "ProbsForG.csv", row.names=F)

head(ProbsForG)
t1 <- ProbsForG

t1$ProbDecline <- ifelse(t1$Negative_Indicator==0, (1-t1$probability), t1$probability)




