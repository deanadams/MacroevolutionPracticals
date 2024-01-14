## ---- sim_problem -----
set.seed(6679881) ##The n for the largest known Cullen prime. Cullen primes have the form n*(2^n)+1 
library(phytools)

#A quick simulation
mytree<- pbtree(n=50, scale=1) #one way to simulate a tree. Note: in R, one has functions to simulate BD trees, random splits trees, etc. using different functions
plot(mytree)
X<-fastBM(tree=mytree) #simulates a continuous trait on phylogeny under Brownian motion
Y<-fastBM(tree=mytree)
cor(X,Y)

## ---- sim_correction ----

X_pic <- pic(x=X,phy=mytree)
Y_pic <- pic(x=Y,phy=mytree)
cor(X_pic,Y_pic)

## ---- read_data1 ----

library(ape)	
library(geiger)

#Here is a large time-dated molecular phylogeny (a chronogram):
ManderTree<-read.tree("../data/Mander.tre",tree.names=T)
plot(ManderTree)


#Read in data, phylogeny, and match/prune them to one another
ManderTree<-read.tree("../data/Mander.tre",tree.names=T)
plot(ManderTree)
Mander_dat<-read.csv("../data/PlethodonMns.csv", header=TRUE, row.names = 1) #Notice we read in the first column as row.names. These MUST match the names in the phylogeny

Mander_dat[,12]<-(as.numeric(Mander_dat[,12]=="N")) ##Convert the groups from strings to numerics. Treedata doesn't like different types of data.

#Now, prune the tree to match the data and vice-versa:
Pleth_matched<-treedata(phy=ManderTree,data = Mander_dat,sort=T, warnings=FALSE)
plot(Pleth_matched$phy)  

#grabs appropriate column from data matrix after treedata. We'll use these in analyses later
SVL<-Pleth_matched$data[,1] 
HL<-Pleth_matched$data[,3]
BodyW<-Pleth_matched$data[,5]
Fore<-Pleth_matched$data[,6]
Hind<-Pleth_matched$data[,7]
Groups<-as.factor(Pleth_matched$data[,12]);names(Groups)<-rownames(Pleth_matched$data) ##Store the groups and name the vector according to the samples



## ---- read_data2 ----


#Now read in data for species in the genus Hydromantes
Hyd_dat<-read.csv("../data/HydromantesMns.csv", header=TRUE, row.names = 1)
Hyd_dat  


Hydro_matched<-treedata(phy=ManderTree,data = Hyd_dat, warnings=FALSE)
plot(Hydro_matched$phy)  #We have matched data to the tree!


## ---- phylo_naive ----
  
#Phylogenetically-naive analyses
plot(BodyW,Hind)
cor.test(Hind,BodyW)  #correlation of two traits: pretty high
anova(lm(Hind~BodyW))

## ---- PIC ----  

#Generate PICs and test while conditioning on phylogeny
Hind_pic<-pic(x = Hind, phy = Pleth_matched$phy) 
BodyW_pic<-pic(x = BodyW, phy = Pleth_matched$phy) 
cor.test(Hind_pic,BodyW_pic)
plot(BodyW_pic,Hind_pic)
anova(lm(Hind_pic~BodyW_pic+0))   #regression through the origin (b/c order of taxa for contrast is arbitrary [can 'spin' on node])
summary(lm(Hind_pic~BodyW_pic+0))   #coefficients of the model

## ---- PGLS ----

#To perform PGLS in R, we must first estimate the phylogenetic covariance matrix V (C in matrix form):
spc <- Pleth_matched$phy$tip.label
V<-corBrownian(phy=Pleth_matched$phy,form = ~spc)
C <- vcv.phylo(phy = Pleth_matched$phy)

#Now we run the analysis:
library(nlme)

bm_gls<-gls(HL~Hind,correlation = V, data=data.frame(Hind, HL))
summary(bm_gls)    
anova(bm_gls)

## ----phylo_transform ----

library(RRPP)
##lm.rrpp, the function used to create the model, requires the data in a specific format. rrpp.data.frame puts the data in that format. It is essentially a data frame in a specific format.
rdf<-rrpp.data.frame(Hind=as.matrix(Hind),Fore=as.matrix(Fore),Groups=as.matrix(Groups), C=C)
res.PhyT<-lm.rrpp(Hind~Fore, data = rdf, Cov = C, print.progress = FALSE)
anova(res.PhyT)   
res.PhyT$LM$gls.coefficients ## the estimated coefficients of the model

## ---- phylo_ANOVA_wrong ----

aov.phylo(SVL~Groups, phy = Pleth_matched$phy) ##Not the same analysis

## ---- phylo_ANOVA
anova(lm.rrpp(SVL~Groups, data = rdf, Cov = C, print.progress = FALSE))
anova(gls(SVL~Groups,correlation = V, data=data.frame(SVL, Groups)))  #identical to Phylo.transform


## ---- phylo_signal ----

####Phylogenetic Signal
library(geomorph)
phylosig(tree=Pleth_matched$phy, x=SVL, method="K", test=T, nsim=1000)  #phytools
  res<-physignal(A=SVL,phy = Pleth_matched$phy, print.progress = FALSE)   #geomorph
summary(res)  #the same, but the latter method (from geomorph) can accomodate multivariate data
plot(res)

## ---- OLS_cor_sim ----
####Additional Simulation Approaches
#Simulate correlated data
library(MASS)
R<-matrix(0.7, nrow=2,ncol=2); diag(R)<-1 ## This matrix desribes the covariation between traits
R
dat.sim<-mvrnorm(n=10000,Sigma = R,mu = c(0,0))
plot(dat.sim, asp=1)


## ---- Phylo_cor_sim ----
#Simulate BM correlated data on phylogeny
mytree<- pbtree(n=100, scale=1) #one way to simulate a tree
R<-matrix(0.7, nrow=2,ncol=2); diag(R)<-1 ## This matrix describes the covariation between traits
R
dat.BM<-sim.char(phy = mytree,par = R,nsim = 1)[,,1]
plot(dat.BM, asp=1)

## ----Practical Questions ----

#To apply what you have learned, please do the following:

#1: Simulate a phylogeny (or read one into R)
#2: Simulate two continuous traits with some known correlation between them
#3: Estimate the evolutionary association of the traits, conditioned on the phylogeny, using one of the methods you have learned

#Upload your code to the assignment center on Canvas 
