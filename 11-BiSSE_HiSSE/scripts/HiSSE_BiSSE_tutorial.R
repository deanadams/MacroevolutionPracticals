## ----load_libraries----
library(diversitree, warn.conflicts=F, quietly=T) 
library(hisse, warn.conflicts=F, quietly=T)
library(geiger, warn.conflicts=F, quietly=T)
library(phytools, warn.conflicts=F, quietly=T)

## ----read_data----
#Read the grunts data
gt <- read.tree(file = "../data/grunts.phy")
print(gt, printlen = 2)

gd <- read.csv(file="../data/grunts.csv", row.names = 1,
                          stringsAsFactors = TRUE)
head(gd)

## ----plot_tree----
hab <-gd[,1]
names(hab) <- rownames(gd)
plotTree(gt, ftype = "i", fsize = 0.7, offset = 1)
tiplabels(pie=to.matrix(hab,0:1)[gt$tip.label,],
          piecol = c("white", "black"), cex = 0.5)
legend("bottomleft", c("non-reef","reef"),
       pch=21,pt.cex = 1.6, pt.bg = c("white", "black"))

## ----bisse_start----
bisse.model <- make.bisse(gt,hab)
p <- starting.point.bisse(gt)
p

## ----bisse_run----
bisse.mle <- find.mle(bisse.model, p)
bisse.mle$par
bisse.mle$lnLik

## ----bisse_null----
bisse.null <- constrain(bisse.model, lambda1 ~ lambda0,
                mu1 ~ mu0)
bisse.null.mle <- find.mle(bisse.null, p[c(-2,-4)])
bisse.null.mle$par
bisse.null.mle$lnLik

## ----bisse_test----
bisseAnova <- anova(bisse.mle, bisse.null.mle)
bisseAnova
aicw(setNames(bisseAnova$AIC, rownames(bisseAnova)))

## ----hisse_setup----
hd <- data.frame(Genus.species=rownames(gd),
                 x=gd[,"habitat"])
head(hd)

rates.bisse <- TransMatMakerHiSSE(hidden.traits = 0)
rates.hisse <- TransMatMakerHiSSE(hidden.traits = 1)
rates.bisse
rates.hisse

## ----hisse_run----
#NOTE: sann = FALSE and bounded.search = FALSE for tutorial runtime
bisse.hmle <- hisse(gt,hd,turnover = c(1,2),
                eps = c(1,2), hidden.states = FALSE,
                trans.rate = rates.bisse, 
                sann = FALSE, bounded.search = FALSE)
#NOTE: sann = FALSE and bounded.search = FALSE for tutorial runtime
hisse.hmle <- hisse(gt,hd,turnover = c(1,1,2,2),
                    eps = c(1,1,2,2), hidden.states = TRUE,
                    trans.rate = rates.hisse,
                    sann = FALSE, bounded.search = FALSE)
bisse.hmle
hisse.hmle
bisse.hmle$AIC
hisse.hmle$AIC
