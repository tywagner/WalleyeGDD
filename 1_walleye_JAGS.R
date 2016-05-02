# rm(list=ls())
library(R2jags)
library(lme4)
library(MCMCpack)
library(arm)
library(lattice)
library(PerformanceAnalytics)
library(MuMIn)

# Read in summary data 
summary <- read.csv('summary_walleye_recruitment.csv')

table(summary$recruitment.surveys)

# Subsample WBICs based on number of recruitment surveys
samp <- summary[ which(summary$recruitment.surveys > 0), ]
# table(samp$recruitment.surveys)


# Read in data
dat1 <- read.csv('walleye_recruitment_GDD.csv')
dim(dat1)
head(dat1)

# Recruitment data subsampled for WBICs with most recruitment surveys
dat <- dat1[dat1$WBIC %in% samp$WBIC,]
dim(dat)

length(unique(dat$WBIC))

# Read in bass data
bassDat <- read.csv('fall_lmb8_cpue_lakemean.csv')

dat <- dat[dat$WBIC %in% bassDat$WBIC,]
dim(dat)

length(unique(dat$WBIC))
head(dat)

# Rename and scale covarites
dat$y <- dat$recruitment.binary
# Log-transform and then standardize GDD
dat$x1 <- log(dat$GDD_wtr_5c) 
dat$x <- (dat$x1 - mean(dat$x1))/sd(dat$x1)
# summary(dat$x)
# dat$area <- as.numeric(scale(log(dat$area.hectares)))
# # hist(dat$area)
# dat$cond <- as.numeric(scale(log(dat$Conductance)))
# dat$bass <- dat$mean.lmb.standardized # Should probably log-transform and then standardize
dat$lake <- dat$Lake
# dat$lat <- as.numeric(scale(dat$latitude))
dat$lakenum <- as.numeric(as.factor(as.numeric(dat$WBIC)))
summary(dat)


# Remove lakes with missing data
dat <- dat[complete.cases(dat),]
length(unique(dat$lakenum))

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

range(dat$area.hectares)
mean(dat$area.hectares)
median(dat$area.hectares)

# Conductivity
cond <- as.numeric(by(dat$Conductance, dat$lakenum, mean)) 
cond <- as.numeric(scale(log(cond)))
hist(cond)

range(dat$Conductance)
mean(dat$Conductance)
median(dat$Conductance)

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
mat1 <- cbind(area, cond, bass, lat, gdd2) # gdd and lat correlated: r = -0.75
cor(mat1)

range(dat$GDD_wtr_5c)
mean(dat$GDD_wtr_5c)
median(dat$GDD_wtr_5c)


# Number of lakes
J <- length(unique(dat$lakenum))


# Load data
data <- list(y = dat$y, group = dat$lakenum, n = dim(dat)[1], J = J,
             x=dat$x, K=K, W=W, z1 = area, z2=cond, z3=bass, z4=gdd2 )


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


betas <- out1$BUGSoutput$mean$BB[,2]
betasCIsL <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.025)) 
betasCIsU <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.975)) 
sigBetas <- betasCIsL * betasCIsU > 0
sum(sigBetas)

betasCIsL80 <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.10)) 
betasCIsU80 <- apply(out1$BUGSoutput$sims.list$BB[,,2],2,quantile, c(0.90)) 
sigBetas80 <- betasCIsL80 * betasCIsU80 > 0
sum(sigBetas80)

hist(betas)
abline(v=0)

# Number of negative and positve slopes
length(which(betas > 0)) # 135
length(which(betas < 0)) # 229

which(betas > 0)
which(betas > 0.95)
which(betas < -1.4)
betas[186]
sigBetas[186]

# covariate order: area, cond, bass, GDD
mean(out1$BUGSoutput$sims.list$mu.alpha) 
mean(out1$BUGSoutput$sims.list$gamma0b) # sig (+)
mean(out1$BUGSoutput$sims.list$gamma0b2) # NS
mean(out1$BUGSoutput$sims.list$gamma0b3) # Sig (-)
mean(out1$BUGSoutput$sims.list$gamma0b4) # NS

quantile(out1$BUGSoutput$sims.list$mu.alpha, c(0.025, 0.975)) 
quantile(out1$BUGSoutput$sims.list$gamma0b, c(0.025, 0.975)) # sig (+)
quantile(out1$BUGSoutput$sims.list$gamma0b2, c(0.025, 0.975)) # NS
quantile(out1$BUGSoutput$sims.list$gamma0b3, c(0.025, 0.975)) # Sig (-)
quantile(out1$BUGSoutput$sims.list$gamma0b4, c(0.025, 0.975)) # NS
# mean(out1$BUGSoutput$sims.list$gamma0b3) 

mean(out1$BUGSoutput$sims.list$mu.beta)
mean(out1$BUGSoutput$sims.list$gamma1b) # NS
mean(out1$BUGSoutput$sims.list$gamma1b2) # Sig (-) - NS 
mean(out1$BUGSoutput$sims.list$gamma1b3) # Sig (-)
mean(out1$BUGSoutput$sims.list$gamma1b4) # NS

quantile(out1$BUGSoutput$sims.list$mu.beta, c(0.025, 0.975)) 
quantile(out1$BUGSoutput$sims.list$gamma1b, c(0.025, 0.975)) # NS
quantile(out1$BUGSoutput$sims.list$gamma1b2, c(0.025, 0.975)) # Sig (-) - NS 
quantile(out1$BUGSoutput$sims.list$gamma1b3, c(0.025, 0.975)) # Sig (-)
quantile(out1$BUGSoutput$sims.list$gamma1b4, c(0.025, 0.975)) # NS
# mean(out1$BUGSoutput$sims.list$gamma1b3) 


# Calculate correlation between varying slopes and intercepts
sd.a <- out1$BUGSoutput$mean$Sigma.B[1,1]
sd.b <- out1$BUGSoutput$mean$Sigma.B[2,2]
cov.ab <- out1$BUGSoutput$mean$Sigma.B[1,2]

corr.ab <- cov.ab/(sd.a * sd.b)
corr.ab

summary <- out1$BUGSoutput$summary
write.csv(summary,'final_model_summary.csv',row.names = T)

wbics <- unique(dat$WBIC)
write.csv(wbics,'final_wbics.csv',row.names = F)

