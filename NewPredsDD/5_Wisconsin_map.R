#########################
######Packages
#########################

library(maps)
library(maptools)
library(rgeos)

coords <- read.csv('wisconsin_wbic_locations.csv')
head(coords)
dim(coords)
coords <- coords[,-1]

length(unique(coords$WBIC))
length(unique(dat$WBIC))

coordsPlot <- coords[coords$WBIC %in% dat$WBIC,]
dim(coordsPlot)
length(unique(coordsPlot$WBIC))

coordsPlot <- unique(coordsPlot[,1:3])
length(unique(coordsPlot$WBIC))

coordsPlot <- coordsPlot[order(coordsPlot$WBIC), ]
head(coordsPlot)

# Create dataframe of WBICs and probs
probsWBIC <- data.frame(probs1,wbics)
as.numeric(probsWBIC$wbics)

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


# Lakes to plot (order in dataframe, not WBIC)
lakePlot1 <- 185
lakePlot2 <- 135
lakePlot3 <- 52
lakePlot4 <- 86

res <- 6
name_figure <- "WI_map.png"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
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
y.polygon <- c( lowerCI.Group[,,lakePlot1][i.for] , upperCI.Group[,,lakePlot1][i.back] )
polygon( x.polygon , y.polygon , col = "lightgray", border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,lakePlot1],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
text(-2.5, 0.95, '(A)')
axis(side=1,labels = FALSE,tck=-0.03)
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=2,line = -.75,lwd=0)

text(-2.5, 0.88, round(probsWBIC[probsWBIC$wbics==wbics[lakePlot1],]$probs1,2))

box(lwd=1)


##### 3 trend figures...pull out and plot the raw data for each one 
plot(x,y,ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
i.for <- order(X3)
i.back <- order(X3, decreasing = TRUE )
x.polygon <- c( X3[i.for] , X3[i.back] )
y.polygon <- c( lowerCI.Group[,,lakePlot2][i.for] , upperCI.Group[,,lakePlot2][i.back] )
polygon( x.polygon , y.polygon , col = "lightgray", border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,lakePlot2],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
text(-2.5, 0.95, '(B)')
axis(side=1,labels = FALSE,tck=-0.03)
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=2,line = -.75,lwd=0)

text(-2.5, 0.88, round(probsWBIC[probsWBIC$wbics==wbics[lakePlot2],]$probs1,2))


box(lwd=1)


plot(x,y,ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
i.for <- order(X3)
i.back <- order(X3, decreasing = TRUE )
x.polygon <- c( X3[i.for] , X3[i.back] )
y.polygon <- c( lowerCI.Group[,,lakePlot3][i.for] , upperCI.Group[,,lakePlot3][i.back] )
polygon( x.polygon , y.polygon , col = "lightgray", border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,lakePlot3],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
text(-2.5, 0.95, '(C)')
axis(side=1,labels = FALSE,tck=-0.03)
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=2,line = -.75,lwd=0)

text(-2.5, 0.88, round(probsWBIC[probsWBIC$wbics==wbics[lakePlot3],]$probs1,2))


box(lwd=1)
mtext(side=2,"Probability of successful walleye recruitment",line=-0.4,at=0.5, outer=T)


plot(x,y,ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
i.for <- order(X3)
i.back <- order(X3, decreasing = TRUE )
x.polygon <- c( X3[i.for] , X3[i.back] )
y.polygon <- c( lowerCI.Group[,,lakePlot4][i.for] , upperCI.Group[,,lakePlot4][i.back] )
polygon( x.polygon , y.polygon , col = "lightgray", border = NA )
# Add posterior means
lines(X3, meanProbGroup[,,lakePlot4],lwd=1, col='black',lty=1)
# title(main="WBIC 8",line=-1,adj=.2)
# axis(side=1,labels = TRUE,tck=-0.03)
# axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
text(-2.5, 0.95, '(D)')
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=2,line = -.75,lwd=0)
axis(side=1,mgp=c(1,0,0),tck= -0.01)

text(-2.5, 0.88, round(probsWBIC[probsWBIC$wbics==wbics[lakePlot4],]$probs1,2))


box(lwd=1)
mtext(side=1,expression(paste('Standardized ', log[e], '-transformed DD')),line=1.0,outer=T, at=0.28)


#Generate left panel of plot with all lakes in WI plotted by cluster ID
map('county',region="Wisconsin",col=c('grey90'),fill=TRUE,resolution = 0,mar=c(0,0,0,0),border="grey50")
points(coordsPlot$lon,coordsPlot$lat,pch=18,col='black',cex=1.2)

# legend('bottomleft',legend = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6","Cluster 7","Cluster 8"),cex=0.75,bty="n",pch=c(15,16,17,18,15,16,17,18),
#        col = c(rgb(1,0,0,1),rgb(0,1,0,1),rgb((255/255),(128/255),(0/255),1),rgb((102/255),(0/255),(204/255),1),rgb(0,0,(204/255),1),
#                rgb(0,0,0,1),rgb((0/255),(255/255),(255/255),1),rgb((255/255),(0/255),(127/255),1)),pt.cex = c(.75,.75,.75,1,.75,.75,.75,1))
map.scale(grconvertX(0.06,"npc"), grconvertY(.08, "npc"),col="black", metric = TRUE, ratio=FALSE, 
          relwidth=0.15,cex=1)
northarrow(loc = c(-87.6,46.3),size = .25,cex = 0.75)

#######Add Arrows to plot##########
# arrows(-95,47.4,-89.2,46.11,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)
# arrows(-95,45.4,-91.21,45.66,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)
# arrows(-95,43.5,-88.9,45.54,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)
# arrows(-95,41,-88.01,43.64,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)

arrows(-95,48.4,coords[coords$WBIC==wbics[lakePlot1],]$lon, coords[coords$WBIC==wbics[lakePlot1],]$lat,
       xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)
arrows(-95,45.4,coords[coords$WBIC==wbics[lakePlot2],]$lon, coords[coords$WBIC==wbics[lakePlot2],]$lat,
       xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)
arrows(-95,42.5,coords[coords$WBIC==wbics[lakePlot3],]$lon, coords[coords$WBIC==wbics[lakePlot3],]$lat,
       xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)
arrows(-95,40, coords[coords$WBIC==wbics[lakePlot4],]$lon, coords[coords$WBIC==wbics[lakePlot4],]$lat,
       xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05, lty=1, lwd=2)


par(def.par)
dev.off()


#########################################
# Mess around to find lakes to plot etc.
which(betas > 0.85)
which(betas > 0.85)
# 173, 185
which(betas < -0.9)
# 52 53 58 59
dat[dat$WBIC==wbics[2275100],]

coords[coords$WBIC==wbics[lakePlot1],]$lon
coords[coords$WBIC==wbics[185],]
coords[coords$WBIC==wbics[52],]
coords[coords$WBIC==wbics[58],]

# probsWBIC[probsWBIC$wbics==wbics[2275100],]

probsWBIC[probsWBIC$wbics==wbics[334],]

probsWBIC[probsWBIC$wbics==59300,]

which(probsWBIC$wbics==2392000)
