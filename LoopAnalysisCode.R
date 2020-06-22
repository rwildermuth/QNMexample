#---------------------------------------------------------
# Author: Robert Wildermuth, rwildermuth@umassd.edu
# Created: 2/24/2016
# Last Modified: 6/22/2020

# Description:
# Imports the interaction matrix for the Bear Island conceptual model (Summerhayes & Elton 1923)
# with hypothetical additions for reindeer to ensure stability
# and adapts the data to a (-1,0,1)-specified community matrix for loop 
# analysis. After the community matrix ('commMat') is created, loops and 
# paths are defined and enumerated, feedback and stability are evaluated, 
# and press perturbation is performed using the make.adjoint() function.
# The total feedback and weighted community matrices are also calculated 
# to determine effect determinacy (essentially number of + or - effects 
# devided by the total number of paths possible). All matrices are then 
# saved to file.

#---------------------------------------------------------

### LoopAnalyst QNMs

library(LoopAnalyst)
library(DiagrammeR)
library(data.table)

# Load the Mental Modeler matrix
mmMat <- read.csv('BearIsland.csv', header = TRUE)
# 22 x 23 data.frame

rNames <- mmMat[, 'X']
mmMat <- as.matrix(mmMat[, -1])
# 22 x 22 matrix
rownames(mmMat) <- rNames

# Translate into a (-1,0,1) community matrix
# intermediate function to assess values for positive or negative influence
influenceFxn <- function(x){
  if(is.na(x)){
    x <- 0
  } else if(x < 0){
    x <- -1
  } else if(x > 0){
    x <- 1
  } else if(x == 0){
    x <- 0
  } else {
    x <- NA # double check that all values are being captured
  }
  return(x)
}

test2 <- sapply(mmMat, influenceFxn)
class(test2)

# Must transpose to keep actors in columns and receivers in rows
commMat <- matrix(test2, nrow = nrow(mmMat), ncol = ncol(mmMat), byrow=TRUE)
rownames(commMat) <- colnames(mmMat)#rNames
colnames(commMat) <- colnames(mmMat)
# 22 x 22 matrix

# Add self-limitation to each component
for(i in 1:nrow(commMat)){
  commMat[i,i] <- -1
}
any(is.na(commMat))

# gives number of loops in the matrix (number of elements in the list)
# each list item contains the loop pathway described by the community 
# members (row or column number in 'commMat')
enumerate.loops(commMat)

# gives number of pathways between community members 'i' and 'j'
enumerate.paths(commMat, i = 3, j = 12) # from mineral salts to dung 

# determine feedback and stability
feedback(commMat)
det(-commMat) # added negative interactions from Reindeer (Hodkinson & Coulson 2004 OIKOS)
#(-1^(nrow(commMat)+1))*det(commMat)
# assess Lyapunov stability (Dambacher et al. 2003 Am Nat)
eigVals <- eigen(commMat, symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)

# save the community matrix
#write.csv(commMat, "BearIslandcommMat2.0.csv')

# negative inverse matrix
#solve(-commMat)

# Get adjoint matrix
adjMat <- make.adjoint(commMat, status = TRUE)

#write.csv(adjMat, 'adjointMat.csv')

# Get absolute feedback matrix
totMat <- make.T(commMat)

#write.csv(totMat, 'totMat.csv')

# get the weighted feedback matrix
weightedFB <- make.wfm(commMat, status = TRUE)

#write.csv(weightedFB, 'weightedFB.csv')


### QPress QNMs
########################################################################
# All code was originally obtained from  Melbourne-Thomas et al. 2012. #
########################################################################
# Adapted from Jon Reum's instructions for IEA QNM Workshop
# by Robert Wildermuth, 8/29/2018

# Download and install the package "QPress" if you haven't already
# Depending on your machine, you may have to go to www.rforge.net/QPress and download and install by hand
#install.packages("https://www.rforge.net/QPress/snapshot/QPress_0.21.tar.gz")
#install.packages('QPress',,'http://www.rforge.net/')

#Load the library
library(QPress)

###########################################################
# Take a look at the Blue King Crab QNM. 
#     A. Download and install the freeware program Dia (https://sourceforge.net/projects/dia-installer/)
#     B. Open up the Dia file
#     C. Links terminating in an arrow indicate a positive effect of a node (origin) on another (arrow terminal)
#        Links terminating in a dot indicate a negative effect.  
#        Links that are solid are "certain" to occur; dashed are "uncertain", but if they do occur their sign is known.
#     D. We can add new nodes, move links around here, etc.. The saved .dia file is loaded into R and analyzed.
###############################################################

#-------------------------------------------------------------------------------
# First create a Dia file to convert from mental modeler CSV to an "edge list" for QPress
# Taken from: https://github.com/NOAA-EDAB/QNM/blob/master/looping_qpress_Rpath.R
#function to create signed digraph from Mental Modeler 
MM2Qpress <- function(data){
  
  mental.sheet <- as.data.table(data)
  names(mental.sheet)[1] <- 'Box'
  n <- nrow(mental.sheet)
  box.names <- names(mental.sheet)[which(names(mental.sheet) != 'Box')]
  model <- c()
  for(i in 1:n){
    pos <- which(mental.sheet[i, .SD, .SDcol = box.names] > 0)
    if(length(pos) > 0){
      pos.interaction <- paste(mental.sheet[i, Box], '->', 
                               names(mental.sheet[i, pos + 1, with = F]), 
                               sep = '')
      model <- append(model, pos.interaction)
    }
    
    neg <- which(mental.sheet[i, .SD, .SDcol = box.names] < 0)
    if(length(neg) > 0){
      neg.interaction <- paste(mental.sheet[i, Box], '-*', 
                               names(mental.sheet[i, neg + 1, with = F]), 
                               sep = '')
      model <- append(model, neg.interaction)
    }
  }
  return(model)
}
#function to create qnm models for qpress from signed digraphs
make.qnm<-function(modmat){
  q<-MM2Qpress(modmat)
  qnm<-dput(q)
  qnm.mod<-parse.digraph(qnm)
  qnm.model<-enforce.limitation(qnm.mod)
}

# Use previously defined community matrix
adjacencyMat <- read.csv("BearIsland.csv")

qpressBI <- make.qnm(adjacencyMat)

write.dia(qpressBI, "BearIsland.dia")

#--------------------------------------------------------------------------

# 1. Load in the Bear Island dia to make an edge list
edges <- model.dia("BearIsland.dia")

## Examine unweighted adjacency matrix
A <- adjacency.matrix(edges, labels=TRUE)

# Take a peak at the adjacency matrix
A

#Visualize A using:
adjacency.image(edges) # RW: this doesn't match -> labling on left inverted

# 2. Build a set of stable matricies, 
#We could add additional validation criteria to filter out matricies that don't reproduce a known system behavior 
sim<- system.simulate(n.sims=10000, edges=edges, 
                      sampler = community.sampler(edges), validators = NULL) 
#The sim object contains the inverse community matrcies, their corresponding edge weights, and a few other things.... 

#Look at the proportion of stable matricies when drawing edge weights from uniform distributions

sim$stable / sim$total #

# 3. Interactively expore how the nodes respond to different press scenarios

impact.barplot(sim)

# Look at how the community responds when nodes are pressed one at a time.
imptable<- impact.table(sim)

#Which node perturbations have similar community outcomes?

imp_dist<- dist(t(imptable))
plot(hclust(imp_dist), main="Perturbation similarity (Euclidean distance)", hang=-1)

imp_dist<- dist(imptable)
plot(hclust(imp_dist), main="Node similarity across perturbations (Euclidean distance)", hang=-1)

