########## SCATTER PLOT FOR LEVEL 2 MODEL FIT#########
out <- out1


# Select random slopes or ints 
mean.beta <- out$BUGSoutput$mean$BB[,1]

######## AREA 
# Fake data to predict
fake1 <- seq(min(area), max(area), length=50)


# Obtain 90% CIs for fitted line
est.lineA <- matrix(NA, ncol=length(fake1), nrow=out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:out$BUGSoutput$n.sims ){
  for(t in 1:length(fake1) ){
    est.lineA[i,t] <- out$BUGSoutput$sims.list$mu.alpha[i] + out$BUGSoutput$sims.list$gamma0b[i] * fake1[t] 
      }
}

# CIs for fitted values
fit1 <- apply(est.lineA, 2, mean )
upper.CIA <- apply(est.lineA, 2, quantile, probs=c(0.975) )
lower.CIA <- apply(est.lineA, 2, quantile, probs=c(0.025) )

## Grab 90% CIs for beta's
u.alpha <- numeric(length(mean.beta) )
l.alpha <- numeric(length(mean.beta) )
for(i in 1:length(mean.beta)) { 
  u.alpha[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.975) )
  l.alpha[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.025) )
}
########


######## Conductivity 
# Fake data to predict
fake1a <- seq(min(cond), max(cond), length=50)

# Obtain 90% CIs for fitted line
est.lineAa <- matrix(NA, ncol=length(fake1), nrow=out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:out$BUGSoutput$n.sims ){
  for(t in 1:length(fake1) ){
    est.lineAa[i,t] <- out$BUGSoutput$sims.list$mu.alpha[i] + out$BUGSoutput$sims.list$gamma0b2[i] * fake1a[t] 
  }
}

# CIs for fitted values
fit1a <- apply(est.lineAa, 2, mean )
upper.CIAa <- apply(est.lineAa, 2, quantile, probs=c(0.975) )
lower.CIAa <- apply(est.lineAa, 2, quantile, probs=c(0.025) )

## Grab 90% CIs for beta's
u.alphaa <- numeric(length(mean.beta) )
l.alphaa <- numeric(length(mean.beta) )
for(i in 1:length(mean.beta)) { 
  u.alphaa[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.975) )
  l.alphaa[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.025) )
}
########

######## BASS 
# Fake data to predict
fake1ab <- seq(min(bass), max(bass), length=50)

# Obtain 90% CIs for fitted line
est.lineAab <- matrix(NA, ncol=length(fake1), nrow=out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:out$BUGSoutput$n.sims ){
  for(t in 1:length(fake1) ){
    est.lineAab[i,t] <- out$BUGSoutput$sims.list$mu.alpha[i] + out$BUGSoutput$sims.list$gamma0b3[i] * fake1ab[t] 
  }
}

# CIs for fitted values
fit1ab <- apply(est.lineAab, 2, mean )
upper.CIAab <- apply(est.lineAab, 2, quantile, probs=c(0.975) )
lower.CIAab <- apply(est.lineAab, 2, quantile, probs=c(0.025) )

## Grab 90% CIs for beta's
u.alphaab <- numeric(length(mean.beta) )
l.alphaab <- numeric(length(mean.beta) )
for(i in 1:length(mean.beta)) { 
  u.alphaab[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.975) )
  l.alphaab[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.025) )
}
########


######## GDD 
# Fake data to predict
fake1abc <- seq(min(gdd2), max(gdd2), length=50)

# Obtain 90% CIs for fitted line
est.lineAabc <- matrix(NA, ncol=length(fake1), nrow=out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:out$BUGSoutput$n.sims ){
  for(t in 1:length(fake1) ){
    est.lineAabc[i,t] <- out$BUGSoutput$sims.list$mu.alpha[i] + out$BUGSoutput$sims.list$gamma0b4[i] * fake1abc[t] 
  }
}

# CIs for fitted values
fit1abc <- apply(est.lineAabc, 2, mean )
upper.CIAabc <- apply(est.lineAabc, 2, quantile, probs=c(0.975) )
lower.CIAabc <- apply(est.lineAabc, 2, quantile, probs=c(0.025) )

## Grab 90% CIs for beta's
u.alphaabc <- numeric(length(mean.beta) )
l.alphaabc <- numeric(length(mean.beta) )
for(i in 1:length(mean.beta)) { 
  u.alphaabc[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.975) )
  l.alphaabc[i] <- quantile(out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.025) )
}
########


####### FIGURE WITH CRI's

res <- 6
name_figure <- "Intercepts_panel.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)

nf <- layout(matrix(c(1:4),nrow=2,ncol=2,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.2,0.1,0.1,0.1),oma=c(1.3,2.5,0,1),mai=c(0.30,0.1,0.05,0) )


size.labels = 1
size.text = 1

x.label = expression(paste('Standardized ', log[e], '- lake area'))
x.label2 = expression(paste('Standardized ', log[e], '- conductance'))
x.label3 = expression(paste('Standardized ', log[e], '- bass CPE'))
x.label4 = expression(paste('Standardized ', log[e], '- average DD'))
y.label <- 'Logit-mean probability of recruitment success'

# y.label <- expression(paste('Logit-mean probability of recruitment success [',alpha[j],']'))

### AREA
plot(mean.beta ~ area,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',
     ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.1,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c( lower.CIA[i.for] , upper.CIA[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(area, mean.beta,pch=16,cex=0.8)

segments(x0=area, x1=area,
         y0=l.alpha, y1=u.alpha, col='black',lwd=1)


lines(fake1,fit1, lwd = 3, col="black", lty = 1)

mtext(x.label, line = 1.12, side = 1, cex = size.text)
# mtext(y.label, line = 1.1, side = 2, cex = size.text)
text(-2.0,6.5,'(A)')

box()

### COND
plot(mean.beta ~ cond,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',
     ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.1,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 , labels=F) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1a)
i.back <- order(fake1a, decreasing = TRUE )
x.polygon <- c( fake1a[i.for] , fake1a[i.back] )
y.polygon <- c( lower.CIAa[i.for] , upper.CIAa[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(cond, mean.beta,pch=16,cex=0.8)

segments(x0=cond, x1=cond,
         y0=l.alphaa, y1=u.alphaa, col='black',lwd=1)


lines(fake1a,fit1a, lwd = 3, col="black", lty = 1)

mtext(x.label2, line = 1.12, side = 1, cex = size.text)
#mtext(y.label, line = 1.1, side = 2, cex = size.text)
text(-2.4,6.5,'(B)')

box()


### BASS
plot(mean.beta ~ bass,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',
     ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.1,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1ab)
i.back <- order(fake1ab, decreasing = TRUE )
x.polygon <- c( fake1ab[i.for] , fake1ab[i.back] )
y.polygon <- c( lower.CIAab[i.for] , upper.CIAab[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(bass, mean.beta,pch=16,cex=0.8)

segments(x0=bass, x1=bass,
         y0=l.alphaab, y1=u.alphaab, col='black',lwd=1)


lines(fake1ab,fit1ab, lwd = 3, col="black", lty = 1)

mtext(x.label3, line = 1.12, side = 1, cex = size.text)
mtext(y.label, line = 1.1, side = 2, cex = size.text, outer=T, at=0.5)
text(-1.3,6.5,'(C)')
box()


### GDD
plot(mean.beta ~ gdd2,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',
     ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.1,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1,labels=F ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1abc)
i.back <- order(fake1abc, decreasing = TRUE )
x.polygon <- c( fake1abc[i.for] , fake1abc[i.back] )
y.polygon <- c( lower.CIAabc[i.for] , upper.CIAabc[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(gdd2, mean.beta,pch=16,cex=0.8)

segments(x0=gdd2, x1=gdd2,
         y0=l.alphaabc, y1=u.alphaabc, col='black',lwd=1)


lines(fake1abc,fit1abc, lwd = 3, col="black", lty = 1)

mtext(x.label4, line = 1.12, side = 1, cex = size.text)
#mtext(y.label, line = 1.1, side = 2, cex = size.text)
text(-2.0,6.5,'(D)')
box()

par(def.par)
dev.off()








