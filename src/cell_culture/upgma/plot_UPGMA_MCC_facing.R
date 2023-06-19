## ---------------------------
##
## Script name: plot UPGMA MCC facing
##
## Purpose of script: create cophyloplots to compare UPGMA trees (as 
## constructed with build_UPGMA_subsample.R) and MCC trees from BEAST2 runs, 
## using cophyloplots package, plotted in the order: UPGMA (left) - MCC (right)
##
## Author: Antoine Zwaans
##
## Date Created: 2023-06-14
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## set working directory for Mac 

setwd("~/typewriter_analysis")    # Antoine's working directory 

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)

## adjusting labels to match the numeric labels as provided in xml alignments for BEAST2. 
## This is done by assigning the corresponding index in the original 100 sample
tree_UPGMA <- read.newick("results/analysis_cell_culture_data/upgma/UPGMAtree_100_13.txt")

new_labels <- tree_UPGMA$tip.label 
for(i in 1:length(cell100)) {
  new_labels[which(new_labels == cell100[i])] <- i - 1 
}
tree_UPGMA$tip.label <- new_labels

# importing the MCC tree
tree_MCC <- ape::read.nexus("results/analysis_cell_culture_data/inference_results/single_clock/MCC.tree")

## plot facing phylogenies
tree_cophylo <- cophylo(tree_UPGMA,tree_MCC)

## open file
png("results/analysis_cell_culture_data/upgma/cophylo_UPGMA_MCC.png",width = 1400, height = 1200)
## create the plot
plot(tree_cophylo,link.type="curved",link.lwd=4,link.lty="solid",link.col=make.transparent("blue",0.25),fsize=0.5,use.edge.length=TRUE)
## close file
dev.off()
