
## ----read_data----
library(geiger, warn.conflicts=F, quietly=T) 
library(phytools, warn.conflicts=F, quietly=T)
library(corHMM, warn.conflicts=F, quietly=T)

tree<-read.tree("../data/tree.64.tre",tree.names=T)
mydata<-read.csv('../data/DiscreteData.csv', row.names=1, header=TRUE)
head(mydata) ##Head can be used to show the first few rows of a data frame 

#Match data with tree
data.pruned<-treedata(phy=tree,data = mydata, warnings=FALSE)
tree<-data.pruned$phy
mydata<-data.pruned$data

## ----plot_data----

## NOTE: see `plotTree.datamatrix` in phytools for other plotting option
plot.phylo(tree,show.tip.label = F) ## setting show.tip.label = F prevents the species names from appearing

## We make a colour key to denote which state gets which color
colorkey1 <- c('red','black') ##trait V1
colorkey2 <- c('green','orange') ##trait V3 
names(colorkey1) <-c(0,1) ##Here we are saying that 0s get red and 1s get black
names(colorkey2) <-c(0,1) ##Here we are saying that 0s get green and 1s get orange
V1cols<-colorkey1[as.character(mydata[,1])] ##putting the trait values in as.character() is important!
V3cols<-colorkey2[as.character(mydata[,3])]


tiplabels(pch = 17,col = V1cols,cex=0.6)
tiplabels(pch = 19,col = V3cols,offset = 0.2,cex=0.6)

## ----single_trait----

plot.phylo(tree,show.tip.label = F)
tiplabels(pch = 17,col = V1cols)

#Set up data. corHMM is picky about how it wants the data formatted
trt1<-cbind(row.names(mydata),mydata[,1])

#Set up initial rate matrices
trait1_model_er <- getStateMat4Dat(trt1,"ER")  #equal transition rates
trait1_model_ard <- getStateMat4Dat(trt1,"ARD")  #All rates different


## ----single_trait_plot
trait1_model_er
trait1_model_ard

plotMKmodel(trait1_model_er$rate.mat,rate.cat = 1)
plotMKmodel(trait1_model_ard$rate.mat,rate.cat = 1)


## ----single_trait_fit----

#Fit models: TRAIT 1
trait1_fit_er <-corHMM(phy = tree, data = trt1, rate.cat = 1, rate.mat = trait1_model_er$rate.mat)
trait1_fit_ard <-corHMM(phy = tree, data = trt1, rate.cat = 1, rate.mat = trait1_model_ard$rate.mat)

##The estimated rates
trait1_fit_er$solution
trait1_fit_ard$solution

#compare models: logL and AIC
trait1_fit_er$loglik
trait1_fit_ard$loglik

##We can formally perform a likelihood ratio test of nested models by doing the following
summary_stat <- -2 * (trait1_fit_er$loglik - trait1_fit_ard$loglik) ## -2 * (model_constrained - model_full)
##This is our p-value
pchisq(summary_stat,df = 1,lower.tail=FALSE ) ##degrees of freedom is the difference in the number of parameters between the two models. In our case 2-1=1

trait1_fit_er$AIC
trait1_fit_ard$AIC


## ---- two_trait ----

##first plot the data
V1cols<-colorkey1[as.character(mydata[,1])] ##putting the trait values in as.character() is important!
V2cols<-colorkey2[as.character(mydata[,2])]
plot(tree,show.tip.label = F)
tiplabels(pch = 17,col = V1cols)
tiplabels(pch = 19,col = V2cols,offset = 0.2)

##we set up the same format of data frame as before but now we include two traits: trait 1 and trait 2
trtset12<-cbind(row.names(mydata),mydata[,1:2])

trait12_model_er <- getStateMat4Dat(trtset12,"ER")  #equal transition rates
trait12_model_ard <- getStateMat4Dat(trtset12,"ARD")  #All rates different

trait12_model_er
trait12_model_ard


## ----two_trait_fixed ----

##coax the full transition matrix by having a data frame with all trait combinations
## First get all the unique states for the traits in question
states1<-unique(mydata[,1])
states2<-unique(mydata[,2])

traits12_expanded<-expand.grid('species',states1,states2) ##Having some string at the beginning is important! Our data needs the first column to be populated with names, this achieves that

##try making matrices again
trait12_model_er <- getStateMat4Dat(traits12_expanded,"ER")  #equal transition rates
trait12_model_ard <- getStateMat4Dat(traits12_expanded,"ARD")  #All rates different

##The full transition matrix. Hot dog!
trait12_model_er
trait12_model_ard

plotMKmodel(trait12_model_er$rate.mat,rate.cat = 1)
plotMKmodel(trait12_model_ard$rate.mat,rate.cat = 1)


## ----two_trait_pagel-----
tr1<-mydata[,1]; names(tr1)<-row.names(mydata)
tr2<-mydata[,2]; names(tr2)<-row.names(mydata)
pagel_fit <- fitPagel(tree,x=tr1,y=tr2)
pagel_fit

## ---- x_dep_y_fit ----
fitPagel(tree,x=tr1,y=tr2,dep.var = "x")

## ---- y_dep_x_fit ----
fitPagel(tree,x=tr1,y=tr2,dep.var = "y")


## ----two_trait_order----

trait_model_er <- getStateMat4Dat(traits12_expanded,"ER")  #Transition matrix with equal transition rates for to two binary traits


##We will use the ER transition matrix as a starting point for our directional model
trait_model_dir<-trait_model_er

##Give some of the transitions their own rate
trait_model_dir$rate.mat[1,2]<-2
trait_model_dir$rate.mat[3,4]<-3

trait_model_dir

## ----two_trait_plot----


#Two trait analysis: Order to transitions matters (changes in trait 2 DEPEND on values of trait 1: Maddison 1990)
#Must define the two models for comparison
plot.phylo(tree,show.tip.label = F)
tiplabels(pie = to.matrix(mydata[,3],sort(unique(mydata[,3]))),piecol=c("red", "black"),cex=.3, offset=0)
tiplabels(pie = to.matrix(mydata[,4],sort(unique(mydata[,4]))),piecol=c("green", "orange"),cex=.3, offset=0.2)

##----two_trait_dir_fit----

traits34<-cbind(row.names(mydata),mydata[,3:4])

trait34_fit_er  <-corHMM(phy = tree, data = traits34, rate.cat = 1, rate.mat = trait_model_er$rate.mat)
trait34_fit_dir <-corHMM(phy = tree, data = traits34, rate.cat = 1, rate.mat = trait_model_dir$rate.mat)

trait34_fit_er$loglik
trait34_fit_dir$loglik



trait34_fit_er$AIC
trait34_fit_dir$AIC ##Strong preference for trait34_fit_dir 

