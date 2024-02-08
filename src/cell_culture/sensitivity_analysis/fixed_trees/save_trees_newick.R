## ---------------------------
##
## Script name: save trees newick
##
## Purpose of script: save a scaled UPGMA and MCC for fixed analyses in BEAST2 in newick format
##
## Author: Antoine Zwaans
##
## Date Created: 2024-02-08
##
## Copyright (c) Antoine Zwaans, 2024
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory to where trees and log files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/")    
source("~/typewriter_analysis/src/cell_culture/upgma/scripts.R")

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
require(ape)


## ----------------------------------------
## Plot UPGMA and MCC against BEAST2 sample
## ----------------------------------------


##read in upgma tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")

##read in log
log <- read.table("combined.log", header = T)

##read in the mcc tree
MCC <- ape::read.nexus(file = "MCC_on_thinned4000000_check.tree")

#get the median posterior tree height
median_posterior_height <- median(log[,"treeHeight.t.alignment"])

##read in scaled upgma
upgma_scaled <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_scaled$edge.length <- upgma_scaled$edge.length * (median_posterior_height/upgma_height)

##save this for fixed analysis
ape::write.tree(upgma_scaled, file='~/typewriter_analysis/results/analysis_cell_culture_data/sensitivity_analysis/fixed_trees/UPGMAtree_1000_scaled_median_heights.txt')
###save MCC for fixed analysis
ape::write.tree(MCC, file='~/typewriter_analysis/results/analysis_cell_culture_data/sensitivity_analysis/fixed_trees/MCC_median_heights.txt')
