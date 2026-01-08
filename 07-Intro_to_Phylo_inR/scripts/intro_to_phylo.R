#R script

#load packages
library(ape)

#create a tree (one way)
tree <- rcoal(5)

tree

class(tree)
names(tree)
tree$edge[1,]
tree$edge.length[1]
tree$tip.label[3]

# Reading and writing trees

write.tree(tree, file = "newick.tre")
newick_tree <- read.tree("newick.tre")
write.nexus(tree, file = "nexus.nex")
nexus_tree <- read.nexus("nexus.nex")

### Plotting trees
plot(tree)
plot(tree, edge.color = rainbow(8))

par(mfrow=c(1,2))
plot(tree, type = "cladogram")
plot(tree, type = "radial")
par(mfrow=c(1,1))

# Other tree plot components
plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()
tree$edge


#Other functions (partial list)
getMRCA(tree, c("t1", "t3"))
node.depth.edgelength(tree)
depths <- node.depth.edgelength(tree)
max(depths) - depths

# Remove species (drop.tip)
tree2 <- drop.tip(phy = tree, tip = "t4")
plot(tree)
plot(tree2)


## Estimating trees in R
library(phangorn)

primates <- read.phyDat("primates.dna", format = "interleaved")
primates

### Parsimony 
primates_start_tree <- random.addition(primates)

parsimony(primates_start_tree, primates)
plot(primates_start_tree)

#branches in number of changes
primates_start_tree <- acctran(primates_start_tree, primates)

plot(primates_start_tree)

#optimize under parsimony
primates_opt <- optim.parsimony(primates_start_tree, primates)
primates_opt <- acctran(primates_opt, primates)

#Parsimony score
parsimony(primates_opt, primates)
plot(primates_opt)

### Maximum Likelihood
#JC
fitJC <- pml(primates_opt, data=primates)
fitJC

#optimize
fitJC <- optim.pml(fitJC, optNni=TRUE)
plot(fitJC$tree)

#GTR
fitGTR <- optim.pml(fitJC, model="GTR", optNni=TRUE)
plot(fitGTR$tree)

#Reroot (will throw error if root not monophyletic)
GTR_tree_rooted <- root(fitGTR$tree,c("Bovine","Mouse"))
try <- root(fitGTR$tree,"Bovine") #ok
plot(try)


fitGTR

#compare fit of models
anova(fitJC, fitGTR)

