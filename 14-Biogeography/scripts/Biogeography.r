## ----read_data ----
library(phytools)
library(biogeobears)

## ----PCA_example ----
PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
plot(PCA, main = "PCA", pch = 21, bg = "red")
