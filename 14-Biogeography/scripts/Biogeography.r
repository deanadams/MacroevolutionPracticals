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

## ----res_dec ----
DEC.fit$optim_result

## ----plot_dec ----
layout(matrix(1:2, 1,2), widths=c(0.2,0.8))
par(mar=c(4.1,0.1,1.3,0.1),cex=0.8)
plot_BioGeoBEARS_results(DEC.fit,
    analysis_titletxt = "DEC model",
    plotlegend = TRUE, tipcex = 0.4,statecex = 0.4)

## ----run_decj ----
dstart <- DEC.fit$outputs@params_table["d","est"]
estart <- DEC.fit$outputs@params_table["e","est"]
jstart <- 0.0001

bgb_run$BioGeoBEARS_model_object@params_table["d",
        "init"] <- dstart
bgb_run$BioGeoBEARS_model_object@params_table["d",
        "est"] <- dstart
bgb_run$BioGeoBEARS_model_object@params_table["e",
        "init"] <- estart
bgb_run$BioGeoBEARS_model_object@params_table["e",
        "est"] <- estart
bgb_run$BioGeoBEARS_model_object@params_table["j",
        "type"] <- "free"
bgb_run$BioGeoBEARS_model_object@params_table["j",
        "init"] <- jstart
bgb_run$BioGeoBEARS_model_object@params_table["j",
        "est"] <- jstart
check_BioGeoBEARS_run(bgb_run)
DEC_J.fit <- bears_optim_run(bgb_run)

## ----res_decj ----
DEC_J.fit$optim_result

## ----plot_decj ----
layout(matrix(1:2, 1,2), widths=c(0.2,0.8))
par(mar=c(4.1,0.1,1.3,0.1),cex=0.8)
plot_BioGeoBEARS_results(DEC_J.fit,
      analysis_titletxt = "DEC+J model",
      plotlegend = TRUE, tipcex = 0.4,statecex = 0.4)

## ----model_compare ----
logL.DEC <- get_LnL_from_BioGeoBEARS_results_object(DEC.fit)
logL.DECJ <- get_LnL_from_BioGeoBEARS_results_object(DEC_J.fit)
AIC.table <- AICstats_2models(logL.DECJ,logL.DEC, numparams1 = 3,
                              numparams2 = 2)
print(AIC.table)
