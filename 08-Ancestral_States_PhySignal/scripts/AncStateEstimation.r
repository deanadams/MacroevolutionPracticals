#Analyses for Tutorial 06: Ancestral State Estimation in R
  #Parts of this tutorial are based on that of L. Revell: http://www.phytools.org/eqg2015/asr.html

## ----read_data----

library(phytools)
library(geiger)
library(geomorph)

######Discrete Anc. State Estimation
## Read data & tree
tree<-read.tree("../data/anole.gp.tre",tree.names=T)
group<-read.csv('../data/anole.gp.csv', row.names=1)
gp<-as.factor(t(group)); names(gp)<-row.names(group)

## ----plot_phylo----

#Plot 
plot(tree,type="fan")
cols<-setNames(palette()[1:length(unique(gp))],sort(unique(gp)))
tiplabels(pie=model.matrix(~gp-1),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

## ----parsimony----

library(castor)
gp.int<-as.numeric(gp); names(gp.int)<-names(gp)  
anc.mp <- asr_max_parsimony(tree,gp.int)

#Plot
plot(tree, type="fan")
cols<-setNames(palette()[1:length(unique(gp.int))],sort(unique(gp.int)))
tiplabels(pie=model.matrix(~as.factor(gp.int)-1),piecol=cols,cex=0.3)
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=anc.mp$ancestral_likelihoods,piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

##---- ml_est ----

library(corHMM)

##Create the model for trait evolution
group_formatted<-cbind(row.names(group),group[,1]) ##put the data in a specific format because the function is picky
group_model_er <- getStateMat4Dat(group_formatted,"ER")  #equal transition rates
plotMKmodel(group_model_er$rate.mat,rate.cat = 1)

##Perform ancestral character estimation
anc.ML<-ace(x = group$x, phy = tree, 
            type = "discrete", method = 'ML',
            model = group_model_er$rate.mat)

anc.ML ## we can see a lot of information just by calling our object

## We can see the probability of each state here, just the first few internal nodes will be shown
round(head(anc.ML$lik.anc),3)



# We can plot the ancestral character estimates with the nodelabels() function
colorkey <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2") ##A colorblind friendly color palette
names(colorkey) <- colnames(anc.ML$lik.anc) ##assign each character state a color
tip_cols<-colorkey[as.character(group$x)] ##putting the trait values in as.character() is important!

plot(tree,type="fan")
nodelabels(pie=anc.ML$lik.anc,piecol=colorkey,cex=0.5)
tiplabels(pch=19,col=tip_cols)
legend(x='bottomleft',legend = names(colorkey),fill=colorkey)


## ---- simmap ----
# simulate single stochastic character map using empirical Bayes method
gp<-setNames(group$x,rownames(group)) ##make.simmap needs a named vector for the data.
tree.smp1<-make.simmap(tree,gp,model=group_model_er$rate.mat)
plot(tree.smp1,colorkey,type="fan")  #one run.  Not overly useful. 

#Must do many times
tree.smp<-make.simmap(tree,gp,model=group_model_er$rate.mat,nsim=100)
anc.smp<-summary(tree.smp,plot=FALSE)
plot(anc.smp, type="fan", fsize=0.8,ftype="i")

rm(list = ls())

## ---- read_data_cont ----

#Data from Mahler et al. 2010. Evolution
tree<-read.tree("../data/anole.svl.tre",tree.names=T)
cont_data<-read.csv('../data/anole.svl.csv', row.names=1)
cont_data_vect<-setNames(cont_data$svl,rownames(cont_data)) ##contMap needs a named vector for the data.

## -----anc_est_cont-----
anc.cont.ML<-fastAnc(tree,cont_data_vect,vars=TRUE,CI=TRUE)
#anc.cont.ML #This would show us the estimate and 95% CI for each internal node
anc.cont.ML$ace  #ancestral estimates

#PLOT as color map
tree.col<-contMap(tree,cont_data_vect,plot=FALSE)  #runs Anc. St. Est. on branches of tree
plot(tree.col,type="fan")

## ----anc_est_cont_ace_ml----

#using 'ace' function in APE
anc.cont.ML2<-ace(x=cont_data_vect,phy=tree, type="continuous", method="ML")
anc.cont.ML2$ace  #with APE
anc.cont.ML$ace   #with phytools: the same estimates


## ----anc_est_cont_ace_gls ----

anc.cont.gls<-ace(x=cont_data_vect,phy=tree, corStruct = corBrownian(1, tree), method="GLS")  #same as ML  (see Schluter et al. 1997)
anc.cont.gls$ace   #GLS: SCP. the same

rm(list = ls())


## ---- known_nodes ----
tree<-pbtree(n=100,scale=1)
## simulate data with a trend
x<-fastBM(tree,internal=TRUE,mu=3)
phenogram(tree,x,ftype="off")

x.tip<-x[match(tree$tip.label,names(x))]
phenogram(tree,x.tip)  #traitgram under BM with no ancestral information

#estimate with no prior 
a<-x[as.character(1:tree$Nnode+Ntip(tree))]
x<-x[tree$tip.label]
## let's see how bad we do if we ignore the trend
plot(a,fastAnc(tree,x),xlab="true values",
     ylab="estimated states under BM")
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red")
title("estimated without prior information")

## incorporate prior knowledge
pm<-setNames(c(1000,rep(0,tree$Nnode)),
             c("sig2",1:tree$Nnode+length(tree$tip.label)))
## the root & two randomly chosen nodes
nn<-as.character(c(length(tree$tip.label)+1,
                   sample(2:tree$Nnode+length(tree$tip.label),2)))
pm[nn]<-a[as.character(nn)]
## prior variance
pv<-setNames(c(1000^2,rep(1000,length(pm)-1)),names(pm))
pv[as.character(nn)]<-1e-100
## run MCMC
mcmc<-anc.Bayes(tree,x,ngen=100000,
                control=list(pr.mean=pm,pr.var=pv,
                             a=pm[as.character(length(tree$tip.label)+1)],
                             y=pm[as.character(2:tree$Nnode+length(tree$tip.label))]))

anc.est<-colMeans(mcmc$mcmc[201:1001,as.character(1:tree$Nnode+length(tree$tip.label))])

plot(a,anc.est,xlab="true values",
     ylab="estimated states using informative prior")
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red")
title("estimated using informative prior")

## ---- PhySignal ----

#Read Data from Mahler et al. 2010. Evolution
tree<-read.tree("../data/anole.svl.tre",tree.names=T)
SVL<-read.csv('../data/anole.svl.csv', row.names=1)
SVL<-setNames(SVL$svl,rownames(SVL)) ##a named vector for the data.

# Run phylogenetic signal
phylosig(tree,SVL,method = "lambda", test = TRUE)
phylosig(tree,SVL,method = "K", test = TRUE)
physignal(SVL, tree) # same in geomorph
geomorph:::physignal.z(SVL, tree)  #effect size


