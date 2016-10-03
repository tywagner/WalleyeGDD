#########################
######Packages
#########################

library(maps)
library(maptools)
# options(device="quartz")


#########################
######Functions
#########################

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
  # checking arguments
  if(missing(loc))  stop("loc is missing")
  if(missing(size))  stop("size is missing")
  # default colors are white and black
  if(missing(cols)) cols <- rep(c("white","black"),8)
  # calculating coordinates of polygons
  radii <- rep(size/c(1,4,2,4),4)
  x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
  y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
  # drawing polygons
  for (i in 1:15) {
    x1 <- c(x[i],x[i+1],loc[1])
    y1 <- c(y[i],y[i+1],loc[2])
    polygon(x1,y1,col=cols[i])
  }
  # drawing the last polygon
  polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
  # drawing letters
  b <- c("E","N","W","S")
  for (i in 0:3) text((.3+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],(.3+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                      cex=cex,col="black")
}


# Grab slopes estimates
slope.ests <- out$BUGSoutput$mean$BB[,2]
# Color ramp for slopes
color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

Col <- color.gradient(slope.ests)

# 
# res <- 6
# name_figure <- "WI_map_MSB_BIO.png"
png(
  "WI_map_MSB_BIO.png",
  width     = 2,
  height    = 2,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
def.par <- par(no.readonly = TRUE)

#Setup the figure 173mm (W) x 120mm (H) 
# dev.new(width=173/25.4,height=120/25.4)
m=cbind(c(1,2,3,4),c(5,5,5,5))
# par(oma=c(.8,.3,.2,.2),family="Arial",ps=10)
# par(oma=c(.8,.3,.2,.2),ps=10)
par(oma=c(2.0,1,.2,.2),ps=10)
layout(m,widths=c(2,2),heights=c(1,1,1,1))
par(mar=c(.5,2,0,0))

##### 3 trend figures...pull out and plot the raw data for each one 
plot(x,y,ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
i.for <- order(X3)
i.back <- order(X3, decreasing = TRUE )
x.polygon <- c( X3[i.for] , X3[i.back] )
y.polygon <- c( lowerCI.Group[,,174][i.for] , upperCI.Group[,,174][i.back] )
polygon( x.polygon , y.polygon , col = sigColorP[i], border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,174],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
# text(-2.5, 0.95, '(A)')
axis(side=1,labels = FALSE,tck=-0.03, cex.axis=0.5)
axis(side=2,labels = FALSE,tck=-0.03, cex.axis=0.5)
axis(side=2,line = -.75,lwd=0, cex.axis=0.5)
box(lwd=1)


##### 3 trend figures...pull out and plot the raw data for each one 
plot(x,y,ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
i.for <- order(X3)
i.back <- order(X3, decreasing = TRUE )
x.polygon <- c( X3[i.for] , X3[i.back] )
y.polygon <- c( lowerCI.Group[,,186][i.for] , upperCI.Group[,,186][i.back] )
polygon( x.polygon , y.polygon , col = sigColorP[i], border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,186],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
# text(-2.5, 0.95, '(B)')
axis(side=1,labels = FALSE,tck=-0.03, cex.axis=0.5)
axis(side=2,labels = FALSE,tck=-0.03, cex.axis=0.5)
axis(side=2,line = -.75,lwd=0, cex.axis=0.5)
box(lwd=1)


plot(x,y,ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
i.for <- order(X3)
i.back <- order(X3, decreasing = TRUE )
x.polygon <- c( X3[i.for] , X3[i.back] )
y.polygon <- c( lowerCI.Group[,,27][i.for] , upperCI.Group[,,27][i.back] )
polygon( x.polygon , y.polygon , col = sigColorP[i], border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,27],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
# text(-2.5, 0.95, '(C)')
axis(side=1,labels = FALSE,tck=-0.03, cex.axis=0.5)
axis(side=2,labels = FALSE,tck=-0.03, cex.axis=0.5)
axis(side=2,line = -.75,lwd=0, cex.axis=0.5)
box(lwd=1)
mtext(side=2,"Probability of successful walleye recruitment",line=-0.4,at=0.5, outer=T, cex=0.4)


plot(x,y,ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
i.for <- order(X3)
i.back <- order(X3, decreasing = TRUE )
x.polygon <- c( X3[i.for] , X3[i.back] )
y.polygon <- c( lowerCI.Group[,,2][i.for] , upperCI.Group[,,2][i.back] )
polygon( x.polygon , y.polygon , col = sigColorP[i], border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,2],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
# axis(side=1,labels = TRUE,tck=-0.03)
# axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
# text(-2.5, 0.95, '(D)')
axis(side=2,labels = FALSE,tck=-0.03, cex.axis=0.5)
axis(side=2,line = -.75,lwd=0, cex.axis=0.5)
axis(side=1,mgp=c(1,0,0),tck= -0.01, cex.axis=0.5)
box(lwd=1)
mtext(side=1,expression(paste('Standardized ', log[e], '-transformed DD')),line=1.0,outer=T, at=0.28, cex=0.4)


#Generate left panel of plot with all lakes in WI plotted by cluster ID
map('county',region="Wisconsin",col=c('grey90'),fill=TRUE,resolution = 0,mar=c(0,0,0,0),border="grey50")
points(dat$longitude,dat$latitude,pch=18,col=Col,cex=1.2)

# legend("bottomleft",title="Slope",legend=c(1:10),col =rbPal(10),pch=20)

# legend('bottomleft',legend = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6","Cluster 7","Cluster 8"),cex=0.75,bty="n",pch=c(15,16,17,18,15,16,17,18),
#        col = c(rgb(1,0,0,1),rgb(0,1,0,1),rgb((255/255),(128/255),(0/255),1),rgb((102/255),(0/255),(204/255),1),rgb(0,0,(204/255),1),
#                rgb(0,0,0,1),rgb((0/255),(255/255),(255/255),1),rgb((255/255),(0/255),(127/255),1)),pt.cex = c(.75,.75,.75,1,.75,.75,.75,1))
map.scale(grconvertX(0.06,"npc"), grconvertY(.08, "npc"),col="black", metric = TRUE, ratio=FALSE, 
          relwidth=0.15,cex=0.4)
northarrow(loc = c(-87.6,46.3),size = .25,cex = 0.5)

#######Add Arrows to plot##########
arrows(-95,47.4,-89.2,46.11,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=1)
arrows(-95,45.4,-91.21,45.66,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=1)
arrows(-95,43.5,-88.9,45.54,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=1)
arrows(-95,41,-88.01,43.64,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=1)

par(def.par)
dev.off()




