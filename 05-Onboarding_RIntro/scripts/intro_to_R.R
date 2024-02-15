# R code for intro to R


#Assign a value to a variable
a <- 3
a

b = 4
b

5 -> c
c

#Combine values into a vector (i.e., array) using the `c` function
b <- c(3,4,5)   
b
b[2]	

#Combine objects into a list using the `list` function
l <- list(number = 3, values = c(3.5, 4, 12), message = "many things can go in a list")  
l
l[[3]] 
l$message


# Generate random deviates
x <- runif(1) 
x

a <- rnorm(50,5,1) 
plot(a)   
hist(a)     

# some built-in summary functions
mean(a)
median(a)
var(a)

## Working with Trees using the _ape_ Package

#install and load packages
#install.packages("ape")
library(ape)

### The tree object
tree <- rcoal(5)
tree

class(tree)
names(tree)
tree$edge[1,]
tree$edge.length[1]
tree$tip.label[3]


#read/write trees
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

### Other useful functions

#Most recent common ancestor
getMRCA(tree, c("t1", "t3"))

#node depth 
node.depth.edgelength(tree)

depths <- node.depth.edgelength(tree)
max(depths) - depths

