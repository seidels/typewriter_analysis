## ---------------------------
##
## Script name: plot_trees_in_2D_treespace_with_likelihoods
##
## Purpose of script:
##
## Author: Antoine Zwaans
##
## Date Created: 2023-11-03
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## set working directory for Mac 
setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells")    

## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)

##sourcing helper functions
source("~/typewriter_analysis/src/cell_culture/upgma/plot_trees_in_2d_with_likelihoods.R")

## Load tree data and extract sample nbrs
sciphy_trees <- ape::read.nexus(file = "thinned4000000.trees")

sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
log <- read.table("combined.log", header = T)

step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- log$Sample[2] - log$Sample[1]

#resample he log file at the same frequency
subsampled_log <- log[seq(1,length(log$Sample),by=step_tree/step_log),]

#resample the same number
subsampled_log <- subsampled_log[1:min(length(sciphy_trees),length(subsampled_log$Sample)),]

#check that all tree sample nrs and likelihood sample nrs match
which(subsampled_log$Sample != sample_nr_tree)


#extract likelihood values
tree_likelihood <- subsampled_log$likelihood

#if needed check low likelihood sciphy_trees and remove them (potential remnants of burnin)
#plot(tree_likelihood)

#get the upgma tree corresponding to the dataset and relabel the tips to match BEAST tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000_13_SERIOUS.txt")
cell_ids <- read.csv(header = F, file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000_13_cell_names_SERIOUS.txt")
cell_ids$numeric_label <- 0:999
cell_ids_sorted <- cell_ids[match(upgma$tip.label, cell_ids$V1), ]
upgma$tip.label <- as.character(cell_ids_sorted$numeric_label)

#create a upgma scaled by the median posterior tree height.
median_posterior_height <- median(log[,"treeHeight.t.alignment"])
upgma_rescaled <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_rescaled$edge.length <- upgma_rescaled$edge.length * (median_posterior_height/upgma_height)

#get the MCC tree 
MCC <- ape::read.nexus(file = "MCC_on_thinned4000000_check.tree")

#create a list of trees to analyse
all_trees <- sciphy_trees
#append upgma and MCC to trees list
all_trees <-c(all_trees, MCC) 
all_trees <-c(all_trees, upgma) 

#checking whether MCC is acutally a tree in the posterior lis
#for(i in 1:length(all_trees)) {

#  if(all.equal(all_trees[[i]],MCC,use.edge.length = FALSE)) {
#   print(i)
#   print("IDENTICAL PHYLO OBJECT")
# }
#}

#save that list:
#ape::write.tree(all_trees, file='tree_list_100Sciphy_MCC_UPGMA.txt')

## -----------------------------------------------------------
## Place all all_trees in 2d using the Robinson Foulds (RF) metric
## -----------------------------------------------------------
res_rf <- treespace(all_trees, nf = 12, "RF")

# Plot RF scenario
tree_df_rf <- res_rf$pco$li
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,2)
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,3)
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,4)
plot2_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",2,3)
plot2_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",2,4)
plot3_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",3,4)

blank_plot <- ggplot() + theme_void()
combined <- cowplot::plot_grid(plot1_2,plot1_3,plot1_4,blank_plot,plot2_3,plot2_4,blank_plot,blank_plot,plot3_4,nrow = 3,ncol=3)
ggsave("combined_treespace_4dim_RF_MCC_UPGMA_TREESPACE.png", combined, width = 60, height = 60, units = "cm", dpi = 300)

## ---------------------------------------------------------------------
## Place all all_trees in 2d using the weighted Robinson Foulds (wRF) metric
## ---------------------------------------------------------------------
all_trees <- c(all_trees,upgma_rescaled) 

res_wrf <- treespace(all_trees, nf = 12, method = "wRF")

# Plot wRF scenario
tree_df_wrf <- res_wrf$pco$li
plot1_2 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",1,2)
plot1_3 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",1,3)
plot1_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",1,4)
plot2_3 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",2,3)
plot2_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",2,4)
plot3_4 <- plot_tree_distances_branch_lengths_metric(tree_df_wrf, tree_likelihood, "Weigthed Robinson Foulds",3,4)

combined <- cowplot::plot_grid(plot1_2,plot1_3,plot1_4,blank_plot,plot2_3,plot2_4,blank_plot,blank_plot,plot3_4,nrow = 3,ncol=3)
ggsave("combined_treespace_4dim_wRF_MCC_UPGMA.png", combined, width = 60, height = 60, units = "cm", dpi = 300)
