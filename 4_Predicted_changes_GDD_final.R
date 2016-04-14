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
ests_early <- cbind(datP$WBIC, mean_early, ci_earlyT)
colnames(ests_early) <- c('WBIC','mean_early','97.5%','2.5%')

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
ests_late <- cbind(datP$WBIC,mean_late, ci_lateT)
colnames(ests_late) <- c('WBIC','mean_late','97.5%','2.5%')
# Save output
saveRDS(ests_early, file="early_estimates.rds")
saveRDS(ests_late, file="late_estimates.rds")

# testrun <- readRDS("late_estimates.rds")
# head(testrun)

num <- dim(datP)[1]

size.labels = 1
size.text = 1

probDiff <- mean_late - mean_early

probDiffColor <- rep('blue',num)
probDiffColor[probDiff < 0] <- 'red'


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

arrows(1:num, mean_early, 1:num, mean_late, col= probDiffColor,code=2, angle=15,length=0.1,lwd=0.8)

mtext('Lake number',side=1,outer=T,adj=0.5,cex=1,line=0.05)
mtext('Probability of successful recruitment',side=2,outer=T,adj=0.5,cex=1.0,line=1.12)
box()
par(def.par)
dev.off()



#####################################################
########### PLOT ####################################
#####################################################
res <- 6
name_figure <- "gdd_change_plot2.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

nf <- layout(matrix(c(1:2), nrow = 1, ncol=2,byrow=TRUE), widths = c(0.7, 0.3))
layout.show(nf)
par(mar = c(1, 1, 1, 0) + 0.2,oma=c(2,2,1,0))
def.par <- par(no.readonly = TRUE)     #

size.labels = 1
size.text = 1

x.label <- 'Probability of successful recruitment'
y.label <- 'Lake number'

# # Data from early model
# rows <- 2060:2122
# Plot.data <- out$summary[rows,]
# plotting.region <- range(Plot.data[,c(3,7)])
# 
# # Data from ignore TP < 1 ug/L analysis
# rows2 <- 2055:2117
# Plot.data2 <- out2$summary[rows2,]
# plotting.region2 <- range(Plot.data2[,c(3,7)])
nlake <- 50 # num
### axis label options
spc <- 0.06
lab <- 1:nlake
cex <- 0.5
adj <- 0
const <- 0.5


###
plot(c(0, 1), c(1-const,nlake + const), 
     axes=F, xlab='',ylab='',type='n')
axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:nlake),labels=F,tck= -0.01)

text(par("usr")[1] - spc,1:nlake,srt = 0, adj =adj,
     labels = lab, xpd = TRUE,cex=0.8)

# 95% CIs for early analysis
segments(x0=ci_early[2,1:nlake], x1=ci_early[1,1:nlake],
         y0=1:nlake, y1=1:nlake, col='black',lwd=1)


## Estiamtes from early model
points(mean_early[1:nlake], 1:nlake, col='black',cex=0.5, pch=16)

## Estimates from late model
segments(x0=ci_late[2,1:nlake], x1=ci_late[1,1:nlake],
         y0=1:nlake+const, y1=1:nlake+const, col='black',lwd=1,lty=2)
points(mean_late[1:nlake], 1:nlake+const, col='firebrick4',cex=0.5, pch=15)
# 95% CIs


# Add x- and y-axis lables
mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T, adj=0.35)
mtext(y.label, line = 0.8, side = 2, cex = size.text, outer=T)

# abline(h=0)
box()

### ADD legend
par(mar = c(0, 0, 1, 2))
plot(1:3, rnorm(3), pch = 1, lty = 1, ylim=c(-2,2), type = "n", axes = FALSE, ann = FALSE)
legend(1, 1, c("Early GDD", "Late GDD"), pch = c(16,15), lty = c(1,2), cex=0.8, col=c('black','firebrick4'),
       x.intersp = 0.5)
###

par(def.par)
dev.off()





