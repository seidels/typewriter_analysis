## ---------------------------
##
## Script name: run_TreeDist_RF 
##
## Purpose of script: Calculate RF pairwise distances and map tree in 12dim
##
## Author: Antoine Zwaans
##
## Date Created: 2023-09-11
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Note: 
##
## script to run TreeDist on euler, relies on sciphy tree files: 
## found in : results/clock_per_target/analysis_cell_culture_data/inference_results/clocks_per_target/1000_cells
## and corresponding UPGMA
## found in :results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt   
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)

#load sciphy trees and UPGMA
sciphy_trees <- ape::read.nexus(file = "thinned4000000.trees")
upgma <- ape::read.tree(file = "UPGMAtree_1000.txt")
cell_ids <- read.csv(header = F, file = "UPGMAtree_1000_cell_names_.txt")
cell_ids$numeric_label <- 0:999
cell_ids_sorted <- cell_ids[match(upgma$tip.label, cell_ids$V1), ]
upgma$tip.label <- as.character(cell_ids_sorted$numeric_label)

#add the upgmas to the sciphy_trees list
all_trees <- sciphy_trees
all_trees[[length(all_trees) + 1]] <- upgma 
names(all_trees)[length(all_trees)] <- "upgma"


#load the MCC tree and add to list
MCC <- ape::read.nexus(file = "MCC_medianHeights.tree")
all_trees[[length(all_trees) + 1]] <- MCC 
names(all_trees)[length(all_trees)] <- "MCC"

#calculate pairwise distances between all trees
distances <- RobinsonFoulds(all_trees)
distances_df <- data.frame(distances)
write_csv(distances_df,"distances_RF.csv")

#use PCoA to map in 12 dimensions 
mapping <- cmdscale(distances, k = 12)
colnames(mapping) <- paste0(rep("A",12),1:12)
mapping_df <- data.frame(mapping)
write_csv(mapping_df,"mapping_df_RF.csv")
