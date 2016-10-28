# rm(list=ls())
library(R2jags)
library(lme4)
library(MCMCpack)
library(arm)
library(lattice)


# Read in summary data 
datP <- read.csv('gdd_change_by_WBIC.csv')

dim(datP)
length(unique(datP$WBIC))
head(datP)

bassDat <- read.csv('fall_lmb8_cpue_lakemean.csv')
dim(bassDat)
datP <- datP[datP$WBIC %in% bassDat$WBIC,]
dim(datP)


datP$x1 <- log(datP$mean.gdd.early) 
datP$x <- (datP$x1 - mean(dat$x1))/sd(dat$x1)

#### Early GDD
gdd_early <- datP$x

pred_early <- array(NA,c(out$BUGSoutput$n.sim, dim(datP)[1] ) )
dim(pred_early)
for(i in 1:dim(datP)[1]){
  pred_early[,i] <- plogis( out$BUGSoutput$sims.list$BB[,i,1] + out$BUGSoutput$sims.list$BB[,i,2] * gdd_early[i] )
}
head(pred_early)


# Posterior mean
mean_early <- apply(pred_early, 2, mean)

# CIs
ci_early <- apply(pred_early, 2, quantile, c(0.975, 0.025))
dim(ci_early)

# Get all ests in one place
ci_earlyT <- t(ci_early)
dim(ci_earlyT)
ests_early <- cbind(datP$WBIC, mean_early, ci_earlyT,sigBetas,sigBetas80)
colnames(ests_early) <- c('WBIC','mean_early','97.5%','2.5%','sig95','sig80')

#### Late GDD
datP$x2 <- log(datP$mean.gdd.late) 
datP$x_late <- (datP$x2 - mean(dat$x1))/sd(dat$x1)

gdd_late <- datP$x_late


pred_late <- array(NA,c(out$BUGSoutput$n.sim, dim(datP)[1] ) )
dim(pred_late)
for(i in 1:dim(datP)[1]){
  pred_late[,i] <- plogis( out$BUGSoutput$sims.list$BB[,i,1] + out$BUGSoutput$sims.list$BB[,i,2] * gdd_late[i] )

}
head(pred_late)

# Posterior mean
mean_late <- apply(pred_late, 2, mean)

# CIs
ci_late <- apply(pred_late, 2, quantile, c(0.975, 0.025))

# Get all ests in one place
ci_lateT <- t(ci_late)
dim(ci_lateT)
ests_late <- cbind(datP$WBIC,mean_late, ci_lateT,sigBetas,sigBetas80)
colnames(ests_late) <- c('WBIC','mean_late','97.5%','2.5%','sig95','sig80')
# Save output
saveRDS(ests_early, file="early_estimates.rds")
saveRDS(ests_late, file="late_estimates.rds")


# testrun <- readRDS("late_estimates.rds")
# head(testrun)

num <- dim(datP)[1]

size.labels = 1
size.text = 1

probDiff <- mean_late - mean_early

# Look at estimates for those with negative slopes
probDiff[which(betas>0)]
betas[21]
betas[which(betas>0)]

probDiffColor <- rep('blue',num)
probDiffColor[probDiff < 0] <- 'red'

sigColor <- rep(0.3,num)
sigColor[sigBetas80==1] <- 1.5
sigColor[sigBetas == 1] <- 4



## Plot
res <- 6
name_figure <- "gdd_change_plot.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(1,1,0.5,1), oma=c(2,2.5,0.1,0.1) )

plot(1:num,mean_early,ylim=c(0,1),xlab='',ylab='', type='n',axes=F)
axis(side=1,cex.axis=size.text, mgp=c(0,0.1,0),tck= -0.008 ) 
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) 

points(1:num,mean_early, pch=1, cex=0.8)
points(1:num,mean_late, pch=6, cex=0.8) # triangle - late

arrows(1:num, mean_early, 1:num, mean_late, col= probDiffColor,code=2, angle=15,length=0.1,lwd=sigColor)

mtext('Lake number',side=1,outer=T,adj=0.5,cex=1,line=0.05)
mtext('Probability of successful walleye recruitment',side=2,outer=T,adj=0.5,cex=1.0,line=1.12)
box()
par(def.par)
dev.off()
