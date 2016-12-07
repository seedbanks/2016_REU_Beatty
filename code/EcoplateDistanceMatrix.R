#Residence Time Chemostat Experiment
#Ecoplate Resource Utilization Distance Matrix/Ordination
#Jaylen Beatty, July 27th, 2016. Indiana University Lennon Lab REU
install.packages("ade4")
install.packages("gclus")
install.packages("FD")
library("ade4")
library("vegan")
library("gclus")
library("cluster")
library("FD")

#Additional functions required to visualize 
#distance matrix. Coldiss.R written by
#Francois Gillet and the authors of Numerical Ecology with R.

#Coldiss citation


#Gillet, Francois. (2012) Coldiss.R [Computer program], License: GPL-2.
#Available in Numerical Ecology with R (2011).


source("coldiss.R")#Function must be downloaded and placed in working directory
#creates heatmaps of dissimilarity matrix.

#Import Data set. Contains average absorbance for 31 resources after subtracting out
#mean absorbance of water control for each chemostat community before and after
#treatment, plus the 2 source communities.

plate =  read.csv("48hourAverageCarbonSubstrateUsage.csv", row.names = 1)
#The plate file includes the chemostat community label and the treatment status.
#"X_1" denotes chemostat X's community pre-treatment, "X_2" denotes the community post-treatment.
#The row.names = 1 part of the function is necessary to have the Bray-Curtis Dissimiliarity Matrix
#Properly label the 

#Q-mode= looking at similarities between pairs of samples
#I.E - Chemostat 1 pre vs. Chemostat 1 post
#Bray-curtis similarity matrix on absorbance values
#Dissimiliartiy converted to distances by D = 1 - 5
plate.bc= vegdist(plate[3:32])#vegdist requires only numerical variables. Column 1 and 2 contain labels and treatment status
#negative values may be altering the accuracy of the distance matrix
head(plate.bc)#Reads out dissimilarity values used to calculate the matrix.
absplate.bc= vegdist(abs(plate[3:32]))
head(absplate.bc)

#Drawing dissimilarity matrix using Numerical Ecology authors' code.
coldiss(plate.bc, byrank=TRUE, diag=TRUE)
coldiss(absplate.bc, byrank=TRUE, diag=TRUE)

#Cluster dendrogram to visualize differences from dissimilarity matrix.
clust.res=hclust(plate.bc)
plot(clust.res, hang= -1)
clust.res2= hclust(absplate.bc)
plot(clust.res2, hang=-1)

install.packages("pvclust")
library(pvclust)
absplate= abs(plate[3:32])
tabsplate= t(absplate)
#t is the transpose function. pvclust reads columns, not rows. Data must be transposed 
#or else a cluster dendrogram would be drawn for each resource, instead of each community.
absfit <- pvclust(tabsplate, method.hclust="ward.D",
               method.dist="euclidean")
plot(fit, hang= -1) # dendogram with p values. hang standardizes branch length.
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)
#The absolute values are being considered because of someone the reads are negative
#negative vlues interfere with the cluster analysis. 

tplate= t(plate)
fit= pvclust (tplate, method.hclust="ward.D",
              method.dist="euclidean")
plot(fit, hang=-1)
pvrect(fit, alpha=.95)
