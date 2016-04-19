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



#################################
#@@@@@@@ EDU-Specific
#################################

# Extract EDU-specific parameters
GroupCoef <- matrix(out1$BUGSoutput$summary[1:(2*(length(unique(dat$lakenum)) )),1], c(length(unique(dat$lakenum)),2), byrow=F)


# Container for predicted values
linPredGroup  <- array(NA, c(out1$BUGSoutput$n.sims,length(X3),J)) 
dim(linPredGroup)

# Put each groups MCMC draws for all 3 parameters in its own list
group.params <- list()
for(m in 1:J){
  group.params[[m]] <- out1$BUGSoutput$sims.list$BB[,m,1:2]
}
# str(group.params)


# group.params[[1]][,1]

for(p in 1:J){ # loop over groups (J)
  for(i in 1:out1$BUGSoutput$n.sims ){  
    for(t in 1:length(X3)){
      linPredGroup[i,t,p] <- group.params[[p]][i,1] + group.params[[p]][i,2]*X3[t] 

    }	  
  }
}

# Convert to probability scale
probPredGroup <- plogis(linPredGroup)

# Store posterior means
meanProbGroup <- array(NA, c(1,length(X3),J) )

upperCI.Group <- array(NA, c(1,length(X3),J) )
lowerCI.Group <- array(NA, c(1,length(X3),J) )

for(i in 1:J){
  # Means
  meanProbGroup[,,i] <- apply(probPredGroup[,,i], 2, mean )

  # 90% CIs for fitted values
  upperCI.Group[,,i] <- apply(probPredGroup[,,i], 2, quantile, probs=c(0.95) )
  lowerCI.Group[,,i] <- apply(probPredGroup[,,i], 2, quantile, probs=c(0.05) )
}

# Save objects for plotting
saveRDS(meanProbGroup, file="meanProbGroup.rds")
saveRDS(upperCI.Group, file="upperCI.Group.rds")
saveRDS(lowerCI.Group, file="lowerCI.Group")

## Color for panel plots
sigColorP <- rep('lightgray',num)
sigColorP[sigBetas80==1] <- 'indianred1'
sigColorP[sigBetas == 1] <- 'lightblue'


################# PLOT ############################
res <- 6
name_figure <- "lake_panel_1.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 0.8
x.label = expression(paste('Standardized ', log[e], '-transformed GDD'))
# x.label = 'Standardized temperature'
y.label = 'Probability of successful walleye recruitment'


nf <- layout(matrix(c(1:100),nrow=10,ncol=10,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Group-specific plots

Ymin <- 0
Ymax <- 1

x <- dat$x
y <- dat$y

# ylabnums <- seq(min(x), max(x), length=5)

for(i in 1:100){
  
  plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
  
  if( i <=90){
	axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
	} else {
	axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	

  if( i ==1 | i==11 | i==21 | i==31 | i==41| i==51 |i==61 | i==71 | i==81 | i==91){
	axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1)
	} else {
	axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  # Add credible region
  i.for <- order(X3)
  i.back <- order(X3, decreasing = TRUE )
  x.polygon <- c( X3[i.for] , X3[i.back] )
  y.polygon <- c( lowerCI.Group[,,i][i.for] , upperCI.Group[,,i][i.back] )
  polygon( x.polygon , y.polygon , col = sigColorP[i], border = NA )
  
  # Add posterior means
  lines(X3, meanProbGroup[,,i],lwd=1, col='black',lty=1)

  text(0.1,0.9,i,cex=0.8)
 #text(0.3,0.9,paste('(',sum2$bkt.n[i],')'),cex=0.8)
  box()
  
}


mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
# 
# box()

par(def.par)
dev.off()




################# PLOT ############################
res <- 6
name_figure <- "lake_panel_2.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 0.8
x.label = expression(paste('Standardized ', log[e], '-transformed GDD'))
# x.label = 'Standardized temperature'
y.label = 'Probability of successful walleye recruitment'


nf <- layout(matrix(c(1:100),nrow=10,ncol=10,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Group-specific plots

Ymin <- 0
Ymax <- 1

x <- dat$x
y <- dat$y

# ylabnums <- seq(min(x), max(x), length=5)

for(i in c(101:200) ){
  
  plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
  
  if( i <=191){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==101 | i==111 | i==121 | i==131 | i==141| i==151 |i==161 | i==171 | i==181 | i==191){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1)
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  # Add credible region
  i.for <- order(X3)
  i.back <- order(X3, decreasing = TRUE )
  x.polygon <- c( X3[i.for] , X3[i.back] )
  y.polygon <- c( lowerCI.Group[,,i][i.for] , upperCI.Group[,,i][i.back] )
  polygon( x.polygon , y.polygon , col = sigColorP[i] , border = NA )
  
  # Add posterior means
  lines(X3, meanProbGroup[,,i],lwd=1, col='black',lty=1)
  
  text(0.1,0.9,i,cex=0.8)
  #text(0.3,0.9,paste('(',sum2$bkt.n[i],')'),cex=0.8)
  box()
  
}


mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
# 
# box()

par(def.par)
dev.off()




################# PLOT ############################
res <- 6
name_figure <- "lake_panel_3.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 0.8
x.label = expression(paste('Standardized ', log[e], '-transformed GDD'))
# x.label = 'Standardized temperature'
y.label = 'Probability of successful walleye recruitment'


nf <- layout(matrix(c(1:100),nrow=10,ncol=10,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Group-specific plots

Ymin <- 0
Ymax <- 1

x <- dat$x
y <- dat$y

# ylabnums <- seq(min(x), max(x), length=5)

for(i in c(201:300) ){
  
  plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
  
  if( i <=291){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==201 | i==211 | i==221 | i==231 | i==241| i==251 |i==261 | i==271 | i==281 | i==291){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1)
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  # Add credible region
  i.for <- order(X3)
  i.back <- order(X3, decreasing = TRUE )
  x.polygon <- c( X3[i.for] , X3[i.back] )
  y.polygon <- c( lowerCI.Group[,,i][i.for] , upperCI.Group[,,i][i.back] )
  polygon( x.polygon , y.polygon , col = sigColorP[i] , border = NA )
  
  # Add posterior means
  lines(X3, meanProbGroup[,,i],lwd=1, col='black',lty=1)
  
  text(0.1,0.9,i,cex=0.8)
  #text(0.3,0.9,paste('(',sum2$bkt.n[i],')'),cex=0.8)
  box()
  
}

# Poplulation-average plot

# plot(x,Y,  axes=F,ylim=c(Ymin,Ymax),xlim=c(min(x), max(x)), ylab='', xlab='', type='n')
# 
# axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
# axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, las=1, labels=F)
# 
# 
# i.for <- order(X3)
# i.back <- order(X3, decreasing = TRUE )
# x.polygon <- c( X3[i.for] , X3[i.back] )
# y.polygon <- c( lowerCIA.PopAve[i.for] , upperCI.PopAve[i.back] )
# polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
# 
# lines(X3,meanProbPopAve, lwd = 1, col="black", lty = 1)
# 
# text(0.2, 0.9, "All data", cex=0.8)
# text(0.6,0.9,'(7021)',cex=0.8)
# 
mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
# 
# box()

par(def.par)
dev.off()



################# PLOT ############################
res <- 6
name_figure <- "lake_panel_4.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 0.8
x.label = expression(paste('Standardized ', log[e], '-transformed GDD'))
# x.label = 'Standardized temperature'
y.label = 'Probability of successful walleye recruitment'


nf <- layout(matrix(c(1:100),nrow=10,ncol=10,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Group-specific plots

Ymin <- 0
Ymax <- 1

x <- dat$x
y <- dat$y

# ylabnums <- seq(min(x), max(x), length=5)

for(i in c(301:364) ){
  
  plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
  
  if( i <=354){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==301 | i==311 | i==321 | i==331 | i==341| i==351 |i==361){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1)
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  # Add credible region
  i.for <- order(X3)
  i.back <- order(X3, decreasing = TRUE )
  x.polygon <- c( X3[i.for] , X3[i.back] )
  y.polygon <- c( lowerCI.Group[,,i][i.for] , upperCI.Group[,,i][i.back] )
  polygon( x.polygon , y.polygon , col = sigColorP[i] , border = NA )
  
  # Add posterior means
  lines(X3, meanProbGroup[,,i],lwd=1, col='black',lty=1)
  
  text(0.1,0.9,i,cex=0.8)
  #text(0.3,0.9,paste('(',sum2$bkt.n[i],')'),cex=0.8)
  box()
  
}

 
mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T, adj=0.7)
mtext(x.label, line = -12.5, side = 1, cex = size.text, outer=T)
# 
# box()

par(def.par)
dev.off()

# rm(linPredPopAve, probPredPopAve, linPredGroup, probPredGroup)
