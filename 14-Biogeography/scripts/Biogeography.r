## ----read_data ----
library(phytools)
library(BioGeoBEARS)
ant.tree <- read.tree("../data/attine-tree-pruned.tre")
print(ant.tree, prinlen = 2)
writeLines(readLines("../data/attine-distribution-data.txt",4))
ant.data <- getranges_from_LagrangePHYLIP(lgdata_fn = 
                "../data/attine-distribution-data.txt")

## ----plot_data ----
tmp <- ant.data@df
tmp[,1:5] <- lapply(tmp[,1:5], factor)
colnames(tmp) <- c("Nearctic","Middle America", 
                   "South America", "Afrotropics",
                   "Australasia")
colors <- setNames(replicate(ncol(tmp),
              setNames(c("white", "darkgray"),0:1),
              simplify = FALSE), colnames(tmp))
object <- plotTree.datamatrix(ant.tree,tmp,fsize=0.5,
            yexp=1.1, header=TRUE, xexp=1.25,colors = colors)                       
legend("topleft", c("species absent","species present"),
       pch=22, pt.bg=c("white","darkgray"), pt.cex=1.3,cex=0.8,bty="n")                       

## ----run_dec ----
max_range_size <- 2
bgb_run <- define_BioGeoBEARS_run(
  num_cores_to_use = 1,
  max_range_size = max_range_size,
  trfn = "../data/attine-tree-pruned.tre",
  return_condlikes_table = TRUE)
bgb_run$geogfn <- "../data/attine-distribution-data.txt"
check_BioGeoBEARS_run(bgb_run)
DEC.fit <- bears_optim_run(bgb_run)
DEC.fit$optim_result

