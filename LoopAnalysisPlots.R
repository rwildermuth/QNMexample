#---------------------------------------------------------
# Author: Robert Wildermuth (rwildermuth@umassd.edu), Gavin Fay (gfay@umassd.edu)
# Created: 6/6/2016/2016
# Last Modified: 6/22/2020

# Description:
# Code from Gavin Fay is modified to create heat maps of the positive 
# and negative influences and weights for loop analysis. Plots 
# are also constructed to summarize positive, negative, and neutral 
# effects on system components.

#---------------------------------------------------------
library(shape)
library(RColorBrewer)

# image.scale() from http://menugget.blogspot.de/2013/12/new-version-of-imagescale-function.html
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}

#---------------------------------------------------------
# Generalized code from Gavin to plot weights and trends together
# RPW: still need to figure out appropriate legend code
PlotLoopResults <- function(adjMat, wMat, rNames){
  #image showing sign of adjoint and weights
  colfunc <- colorRampPalette(c("white", "steelblue"))
  layout(matrix(c(1,1,2,3), nrow=2, ncol=2), widths=c(4,lcm(3)), heights=c(1, 4))
  layout.show(3)
  par(mar=c(1.1,20,16,1.1))
  image(1:nrow(adjMat),1:nrow(adjMat),t(wMat[nrow(adjMat):1,]),col = colfunc(10), ylab='', xlab='', xaxt='n', yaxt='n')
  #colorlegend(col = colfunc(10), zlim = c(0,1), zlevels = 11, dz = 0.1,
  #            zval = seq(0,1,by=0.1), log = FALSE, posx = c(0.9, 0.93), 
  #            posy = c(0.05, 0.9), main = NULL, main.cex = 1.0, 
  #            main.col = "black", lab.col = "black", 
  #            digit = 1, left = FALSE)
  # add axes with component lables
  axis(3, at=1:nrow(adjMat), labels=FALSE, tck=0)
  text(x=1:nrow(adjMat), y=par("usr")[2]+0.5,
       labels=rNames, srt=-45, pos = 2, offset = -0.25, xpd=NA, cex = 1.5)
  #axis(side = 1,
  #   at = 1:nrow(adjMat),
  #   labels = rNames,
  #   tck=0, las=2, cex.axis=1.5)
  axis(side = 2,
     at = nrow(adjMat):1,
     labels = rNames,
     tck=0, las=2, cex.axis=1.5)
  #text(1,0,"+",col="black",cex=2)
  # index of entries in 'adjMat' that are positive
  pick <- which(adjMat>0)
  # add the positive effects
  x <- ceiling(pick/nrow(adjMat))
  y <- (nrow(adjMat)+1)-pick%%nrow(adjMat)#31
  y[y==(nrow(adjMat)+1)] <- 1
  text(x,y,"+",col="black",cex=1.5)
  # index of entries in 'adjMat' that are negative
  pick <- which(adjMat<0)
  # add the negative effects
  x <- ceiling(pick/nrow(adjMat))
  y <- (nrow(adjMat)+1)-pick%%nrow(adjMat)#31
  y[y==(nrow(adjMat)+1)] <- 1
  text(x,y,"-",col="black",cex=1.5)
  
  # index of entries in 'adjMat' that are neutral
  #pick <- which(adjMat==0)
  # add the negative effects
  #x <- ceiling(pick/nrow(adjMat))
  #y <- (nrow(adjMat)+1)-pick%%nrow(adjMat)#31
  #y[y==(nrow(adjMat)+1)] <- 1
  #text(x,y,".",col="black",cex=1.5)
  par(mar=c(0.5,0.5,0.5,0.5))
  plot.new()
  #text(x=0.5, y=0.75, labels="+  Positive\n   Response", cex=1.5)
  #text(x=0.5, y=0.5, labels="-  Negative\n    Response", cex = 1.5)
  
  par(mar=c(1.1,0.5,3,6))
  image.scale(z=t(wMat[nrow(adjMat):1,]), col = colfunc(10), 
              breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), axis.pos=4, add.axis=TRUE)
  mtext("Reliability Weight", side=4, line=3, cex=1.25)
 
}

#---------------------------------------------------------

# create plot with the Bear Island model
adjMat <- read.csv('adjointMat.csv', header=TRUE)
rNames <- adjMat[,1]
adjMat <- as.matrix(adjMat[,2:23])
row.names(adjMat) <- rNames
weightMat <- read.csv('weightedFB.csv', header=TRUE)
weightMat <- as.matrix(weightMat[,2:23])
row.names(weightMat) <- rNames

# Plot overall single-node positive press perturbations
PlotLoopResults(adjMat = adjMat, wMat = weightMat, rNames = rNames)



#---------------------------------------------------------
# Revised Fig 2 plot (alternative style)
# Modified 9/14/2017
simpAdjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/simpAdjMat2.0.csv',
                       header=TRUE)
rNames <- simpAdjMat[,1]
simpAdjMat <- as.matrix(simpAdjMat[,2:27])
row.names(simpAdjMat) <- rNames

simpWMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/simpWeight2.0.csv',
                     header=TRUE)
simpWMat <- as.matrix(simpWMat[,2:27])
row.names(simpWMat) <- rNames

# implement colors that are colorblind-friendly and print friendly from http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
colKey <- c("#beaed4", "#fdc086", "#7fc97f")
xMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)
colMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)
pchMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)

# create matrix indicating character types and colors
for(i in 1:ncol(simpWMat)){
  for(j in 1:nrow(simpWMat)){
    # if(simpWMat[j,i]<0.5){
    #   pchMat[j,i] <- 21
    # } else {
    #   pchMat[j,i] <- 16
    # }
    
    if(simpAdjMat[j,i] < 0){
      colMat[j,i] <- colKey[1]
      pchMat[j,i] <- 25
    } else if(simpAdjMat[j,i] > 0){
      colMat[j,i] <- colKey[3]
      pchMat[j,i] <- 24
    } else {
      colMat[j,i] <- colKey[2]
      pchMat[j,i] <- 21
    }
    
  }
  xMat[, i] <- i
}
png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/SimpleAdjPlot.png",
    width = 15, height = 11, units = 'in', res = 800)
layout(matrix(c(1,1,2,3), nrow=2, ncol=2), widths=c(5,1), heights=c(1, 4))
layout.show(3)
par(mar=c(1,20.5,16,0.1))
plot(as.vector(t(xMat)), rev(rep(c(1:nrow(simpWMat)), each = nrow(simpWMat))), pch = as.vector(t(pchMat)), 
     col = as.vector(t(colMat)), bg = as.vector(t(colMat)), cex = as.vector(t(simpWMat))+1,
     ylab='', xlab='', xaxt='n', yaxt='n')
axis(3, at=1:nrow(simpAdjMat), labels=FALSE, tck=0)
text(x=1:nrow(simpAdjMat), y=par("usr")[2]+0.5,
     labels=c("Tidal Forcing", "Winds", "Air Temperature", "Source Water Proportions", "Precipitation", 
              "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity",
              "Stratification", "Habitat: Pelagic", "Habitat: Nearshore", "Habitat: Seafloor & Demersal", 
              "Protected Species", "* Forage Fish", "Groundfish", "Fished Invertebrates", "Copepods & Micronekton",
              "Gelatinous Zooplankton", "Benthos", "Mid Atlantic Groundfish", "Primary Production", 
              "Detritus & Bacteria", "Recreational Groundfish Fishery", "Commercial Fishery", 
              "Cultural Practices & Attachments"), 
     srt=-45, pos = 2, offset = -0.25, xpd=NA, cex = 1.5)
axis(side = 2,
     at = nrow(simpAdjMat):1,
     labels = c("Tidal Forcing", "Winds", "Air Temperature", "Source Water Proportions", "Precipitation", 
                "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity",
                "Stratification", "* Habitat: Pelagic", "* Habitat: Nearshore", "* Habitat: Seafloor & Demersal", 
                "* Protected Species", "* Forage Fish", "* Groundfish", "* Fished Invertebrates", "Copepods & Micronekton",
                "Gelatinous Zooplankton", "Benthos", "Mid Atlantic Groundfish", "Primary Production", 
                "Detritus & Bacteria", "* Recreational Groundfish Fishery", "* Commercial Fishery", 
                "* Cultural Practices & Attachments"),
     tck=0, las=2, cex.axis=1.5)
polygon(x = c(11.5, 11.5, 13.5, 13.5), y=c(0.5,26.5,26.5,0.5), lty = 2)
polygon(x = c(24.5, 24.5, 25.5, 25.5), y=c(0.5,26.5,26.5,0.5), lty = 2)
#brackets(11.5, -2, 13.5, -2, lwd=2, type=1)
par(mar=c(0.5,0.5,0.5,0.5))
plot.new()
par(mar=c(20,0.1,20,6.5))
plot(x=rep(c(1,2,3),3), y=rep(c(1.5,2.5,3.5), each=3), pch=rep(c(25,21,24), each=3), 
     col=rep(colKey, each=3), bg=rep(colKey, each=3), cex=rep(c(2,1.5,1), 3), ylab='', xlab='', xaxt='n', yaxt='n', bty="n", 
     xlim = c(0,4), ylim=c(1,5))
axis(4, at=c(1.5, 2.5, 3.5), labels=c("Decrease", "No change", "Increase"), tck=0, las=1, cex.axis=1.25)
axis(1, at=1:3, labels=c(1, 0.5, 0), tck=-0.05, las=1, cex.axis=1.25)
mtext("Reliability\nWeight", side=1, line=4, cex=1.25)
dev.off()
#---------------------------------------------------------

