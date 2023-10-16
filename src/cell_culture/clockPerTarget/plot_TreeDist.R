## ---------------------------
##
## Script name: plot_TreeDist 
##
## Purpose of script: Plot outputs of 2D mapping obtained with TreeDist and assess quality
##
## Author: Antoine Zwaans
##
## Date Created: 2023-09-11
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory for Mac 

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/TreeDist")    

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)
library(gdata)

#input RF 2D mapping 
RF_treeDist <- read.csv( file = "mapping_df_RF.csv")

#input original distances
distances <- read.csv( file = "distances_RF.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=651, nrow=651) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances

#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(RF_treeDist)), function(k) {
  newDist <- dist(RF_treeDist[, seq_len(k)])
  MappingQuality(mat, newDist, 10)["TxC"]
}, 0)

#plot the mapping quality
png("mapping_quality_RF.png")
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

#### preprocess log to extract the correct likelihood values to plot ####
sciphy_trees <- ape::read.nexus(file = "~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/thinned4000000.trees")

sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
log <- read.table("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/combined.log", header = T)

step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- log$Sample[2] - log$Sample[1]

#resample the log file at the same frequency
subsampled_log <- log[seq(1,length(log$Sample),by=step_tree/step_log),]

#resample the same number
subsampled_log <- subsampled_log[1:min(length(sciphy_trees),length(subsampled_log$Sample)),]

#check that all tree sample nrs and likelihood sample nrs match
which(subsampled_log$Sample != sample_nr_tree)

#extract likelihood values
tree_likelihood <- subsampled_log$likelihood

#plot with likelihoods as previously
plot <- plot_tree_distances_topology_metric(RF_treeDist, tree_likelihood, "Robinson Foulds")
#remove the legend, the other plit will have it
plot_rf <- plot + theme(legend.position = "none",text=element_text(size = 15),panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 
ggsave("2d_likelihood_mds_RF_medianHeights_RF_treedist.png", plot_rf, width = 15, height = 15, units = "cm", dpi = 1000)

#input PhyloInfo 2D mapping 
PI_treeDist <- read.csv( file = "mapping_df_PI.csv")

#input original distances
distances <- read.csv( file = "distances_PI.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=651, nrow=651) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances

#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(RF_treeDist)), function(k) {
  newDist <- dist(RF_treeDist[, seq_len(k)])
  MappingQuality(mat, newDist, 10)["TxC"]
}, 0)

#plot the mapping quality
png("mapping_quality_PI.png")
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

plot <- plot_tree_distances_topology_metric(PI_treeDist, tree_likelihood, "Phylogenetic information")

#remove the legend, the other plit will have it
plot_pi <- plot + theme(legend.position = c(0.8,0.9),text=element_text(size = 15),panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 
ggsave("2d_likelihood_mds_RF_medianHeights_PI_treedist.png", plot_pi, width = 15, height = 15, units = "cm", dpi = 1000)
combined <- cowplot::plot_grid(plot_rf,plot_pi,nrow = 1) 
ggsave("combined_TreeDist.png", combined, width = 40, height = 20, units = "cm", dpi = 600)
