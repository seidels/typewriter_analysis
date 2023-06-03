## ---------------------------
##
## Script name: plot_trees_in_2d
##
## Purpose of script: Show similarity of UPGMA tree to
## the trees in the posterior set with respect to tree
## topology only (measured by Robinson Foulds) and with
## respect to tree topology and branch lengths (weighted
## Robinson Foulds.
##
## Author: Sophie Seidel
##
## Date Created: 2023-04-03
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## Set working directory for Mac
setwd("~/Projects/typewriter_analysis/")      # Sophie's working directory (mac)

## ---------------------------

## Load required packages
library(treespace)
library(adegraphics)
library(ggplot2)
library(rgl)
library(phangorn)

## ---------------------------

## Define functions

# Function to plot tree distances
plot_tree_distances <- function(tree_df, filename) {
  tree_df$label <- as.factor(c(rep("posterior", length(trees) - 1), "UPGMA"))
  plot <- ggplot(tree_df, aes(x = A1, y = A2)) +
    geom_point(size = 6, shape = 1, colour = "gray50") +
    geom_point(aes(col = label), size = 6, alpha = 0.2) +
    xlab("") + ylab("") + theme_bw(base_family = "")+
    theme(legend.position = c(0.16, 0.93),
          legend.title = element_blank(),
          text = element_text(size = 14))
  ggsave(filename, plot, width = 15, height = 15, units = "cm", dpi = 300)
}

## ---------------------------

## Load tree data

upgma <- ape::read.tree(file = "results/analysis_cell_culture_data/upgma/UPGMAtree_100_12.txt")
cell_ids <- read.csv(header = F, file = "~/Desktop/cell_ids_seed1.txt")
cell_ids$numeric_label <- 0:99
cell_ids_sorted <- cell_ids[match(upgma$tip.label, cell_ids$V1), ]
upgma$tip.label <- as.character(cell_ids_sorted$numeric_label)

trees <- ape::read.nexus(file = "results/analysis_cell_culture_data/typewriter_model_real_12_100_seed1.tree.resample5000000.1.trees")
trees[[length(trees) + 1]] <- upgma

## ---------------------------

## Using the Robinson Foulds metric
res_rf <- treespace(trees, nf = 2, method = "RF")

# Plot RF scenario
tree_df_rf <- res_rf$pco$li
plot_tree_distances(tree_df_rf, "src/cell_culture/upgma/plots/upgma_posterior_100cells_mds_RF.png")

## ---------------------------

## Using the weighted Robinson Foulds (wRF) metric
res_wrf <- treespace(trees, nf = 2, method = "wRF")

# Plot wRF scenario
tree_df_wrf <- res_wrf$pco$li
plot_tree_distances(tree_df_wrf,"src/cell_culture/upgma/plots/upgma_posterior_100cells_mds_wRF.png")

