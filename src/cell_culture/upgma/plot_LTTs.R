/## ---------------------------
##
## Script name: plot_LTTs
##
## Purpose of script: plot UPGMA ltt vs BEAST2 outputs
##
## Author: Antoine Zwaans
##
## Date Created: 2023-07-21
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory for Mac 

## Set working directory to where the log and tree files are
setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/single_clock/1000_cells/")  

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)

## ----------------------
## Plot MCC against UPGMA
## ----------------------

#load individual trees
tree_MCC <- ape::read.nexus(file = "MCC_CommonAncestorHeights.tree")
tree_UPGMA <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")

trees <- c(tree_MCC,tree_UPGMA)
names(trees) <- c("MCC","UPGMA")

#simple built-in LTT: 
png(file="LTT_comp.png",
    width=600, height=600)
ape::mltt.plot(trees ) 
dev.off()


## ----------------------------------------
## Plot UPGMA and MCC against BEAST2 sample
## ----------------------------------------

##read in upgma and mcc trees
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")
MCC <- ape::read.nexus(file = "MCC_CommonAncestorHeights.tree")

##adjust tip labels
cell_ids <- read.csv(header = F, file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000_cell_names_.txt")
cell_ids$numeric_label <- 0:999
cell_ids_sorted <- cell_ids[match(upgma$tip.label, cell_ids$V1), ]
upgma$tip.label <- as.character(cell_ids_sorted$numeric_label)

##read in SciPhy trees
trees <- ape::read.nexus(file = "smallest_combined.trees")

##Plot them all 
depth.range <- range(unlist(lapply(trees,ape::branching.times)), unlist(lapply(list(upgma),ape::branching.times)))
max.depth <- sum(abs(depth.range))

png(file="LTT_comp_all_MCC_UPGMA.png",
    width=30, height=12,units = "cm",res = 300)
par(mar = c(4.1, 4.1, 1.1, 1.1))
plot(x=c(0, -1*max.depth), y=c(1, ape::Ntip(trees[[1]])), type="n", bty="n", xaxt="n", xlab="Time (days)", ylab="Number of lineages",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
colors=c("#E1E1F7",rgb(1,0,0,0.5),"#060647")
list.of.both <- list(trees, list(upgma),list(MCC))
for (i in sequence(3)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
axis(1, at=seq(-25, 0, by=5), labels = seq(0, 25, by=5), cex.axis=1.5)
legend(x=-5,y= 200, legend=c("SciPhy", "UPGMA","MCC"), fill=colors)
dev.off()
