
## ----read_data----
library(geiger, warn.conflicts=F, quietly=T) 
library(phytools, warn.conflicts=F, quietly=T)
library(corHMM, warn.conflicts=F, quietly=T)

#Read the bonyfish data
bonyfish.tree <- read.tree(file = "../data/bonyfish.tre")
print(bonyfish.tree, printle = 3)

bonyfish.data <- read.csv(file="../data/bonyfish.csv", row.names = 1,
                          stringsAsFactors = TRUE)
head(bonyfish.data)

## ----plot_bonyfish_data----
object <- plotTree.datamatrix(bonyfish.tree, bonyfish.data, fsize = 0.5, 
                            yexp = 1, header = FALSE, xecp = 1.45,
                            palettes = c("YlOrRd", "PuBuGn"))
leg <- legend(x = "topright", names(object$colors$spawning_mode),
              cex = 0.7, pch = 22, pt.bg = object$colors$spawning_mode,
              pt.cex = 1.5, bty = "n", title = "spawning mode")
leg<-legend(x=leg$rect$left+4.7,y=leg$rect$top-leg$rect$h,
            names(object$colors$paternal_care),cex=0.7,
            pch=22,pt.bg=object$colors$paternal_care,pt.cex=1.5,
            bty="n",title="paternal care")

## ----Association_Test----

spawning_mode<-setNames(bonyfish.data[,1],
                        rownames(bonyfish.data))
paternal_care<-setNames(bonyfish.data[,2],
                        rownames(bonyfish.data))

interdependent<-fitPagel(bonyfish.tree,paternal_care,
                         spawning_mode)
print(interdependent)

plot(interdependent,cex.main=1,cex.sub=0.8,
     cex.traits=0.7,cex.rates=0.7,
     lwd.by.rate=TRUE,max.lwd=6)

## ----Directional_Tests----

dependent_care<-fitPagel(bonyfish.tree,paternal_care,
                         spawning_mode,dep.var="x")

dependent_spawning<-fitPagel(bonyfish.tree,paternal_care,
                             spawning_mode,dep.var="y")

anova(dependent_care,dependent_spawning,interdependent)
