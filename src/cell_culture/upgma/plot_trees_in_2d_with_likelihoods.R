## ---------------------------
##
## Script name: plot_trees_in_2d_with_likelihoods
##
## Purpose of script: Show similarity of UPGMA tree to
## the trees in the posterior set with respect to tree
## topology only (measured by Robinson Foulds) and with
## respect to tree topology and branch lengths (weighted
## Robinson Foulds.
##
## Based on Sophie Seidel's plot_trees_in_2d.R script
##
## Date Created: 2023-07-21
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## Set working directory to where the log and tree files are
setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/single_clock/1000_cells/")     

## ---------------------------

## Load required packages
library(treespace)
library(adegraphics)
library(ggplot2)
library(rgl)
library(phangorn)

## ---------------------------

## Define functions

# Function to plot tree distances colored by likelihood
plot_tree_distances <- function(tree_df,tree_likelihood,metric_name) {
 tree_df <-  tree_df
 plot_title = paste(metric_name,"metric space")
  
  plot <- ggplot(tree_df[1: (length(tree_df$A1) - 2),], aes(x = A1, y = A2)) +
    geom_point(aes(col = tree_likelihood), size = 6, alpha = 0.5) + 
    scale_color_gradient(low = "#E1E1F7", high = "#060647",name = "Likelihood")  +  
    xlab("") + ylab("") + theme_bw(base_family = "")+ 
    theme(legend.position = c(0.8,0.2), text = element_text(size = 14)) +
    geom_point(aes(x=tree_df$A1[length(tree_df$A1)-1],y=tree_df$A1[length(tree_df$A2)-1]),colour="red") +
    geom_text(aes(x=tree_df$A1[length(tree_df$A1)-1]*1.05,y=tree_df$A1[length(tree_df$A2)-1]*1.05,label="UPGMA")) +
    geom_point(aes(x=tree_df$A1[length(tree_df$A1)],y=tree_df$A1[length(tree_df$A2)]),colour="red") +
    geom_text(aes(x=tree_df$A1[length(tree_df$A1)]*(0.95),y=tree_df$A1[length(tree_df$A2)]*(0.95),label="MCC")) + ggtitle(plot_title)
  
  return(plot)
}

## ---------------------------

## Load tree data and extract sample nbrs
trees <- ape::read.nexus(file = "smallest_combined.trees")
sample_nr_tree <- as.numeric(unlist(strsplit(names(trees),"_"))[seq(2,2*length(trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
log <- read.table("combined.log", header = T)

step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- log$Sample[2] - log$Sample[1]

#resample he log file at the same frequency
subsampled_log <- log[seq(1,length(log$Sample),by=step_tree/step_log),]

#resample the same number
subsampled_log <- subsampled_log[1:min(length(trees),length(subsampled_log$Sample)),]

#check that all tree sample nrs and likelihood sample nrs match
subsampled_log$Sample == sample_nr_tree

#extract likelihood values
tree_likelihood <- subsampled_log$likelihood

#get the upgma tree corresponding to the dataset and relabel the tips to match BEAST tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")
cell_ids <- read.csv(header = F, file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000_cell_names.txt")
cell_ids$numeric_label <- 0:999
cell_ids_sorted <- cell_ids[match(upgma$tip.label, cell_ids$V1), ]
upgma$tip.label <- as.character(cell_ids_sorted$numeric_label)

#add the upgma to the trees list
trees[[length(trees) + 1]] <- upgma 
names(trees)[length(trees)] <- "upgma"

#get the MCC tree and add to the trees list
MCC <- ape::read.nexus(file = "MCC_commonAncestorHeights.tree")
trees[[length(trees) + 1]] <- MCC 
names(trees)[length(trees)] <- "MCC"

## -----------------------------------------------------------
## Place all trees in 2d using the Robinson Foulds (RF) metric
## -----------------------------------------------------------

res_rf <- treespace(trees, nf = 2, method="RF")

# Plot RF scenario
tree_df_rf <- res_rf$pco$li
plot <- plot_tree_distances(tree_df_rf, tree_likelihood, "Robinson Foulds")
#remove the legend, the other plit will have it
plot_rf <- plot + theme(legend.position = "none") 
ggsave("2d_likelihood_mds_RF.png", plot_rf, width = 15, height = 15, units = "cm", dpi = 300)

## ---------------------------------------------------------------------
## Place all trees in 2d using the weighted Robinson Foulds (wRF) metric
## ---------------------------------------------------------------------

res_wrf <- treespace(trees, nf = 2, method = "wRF")

# Plot wRF scenario
tree_df_wrf <- res_wrf$pco$li
plot_wrf <- plot_tree_distances(tree_df_wrf, tree_likelihood, "Weighted Robinson Foulds")
ggsave("2d_likelihood_mds_wRFl.png", plot_wrf, width = 15, height = 15, units = "cm", dpi = 300)

#make a combined version of these plots
combined <- cowplot::plot_grid(plot_rf,plot_wrf,nrow = 1)
ggsave("combined_2D_plot.png", combined, width = 30, height = 15, units = "cm", dpi = 300)
