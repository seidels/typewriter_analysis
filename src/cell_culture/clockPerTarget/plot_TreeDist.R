## ---------------------------
##
## Script name: plot_TreeDist 
##
## Purpose of script: Plot outputs of 2D mapping obtained with TreeDist and assess quality.
## Relies on outputs from run_TreeDist_XX.R, unmapped distances and mapped distance run on cluster.
##
## Author: Antoine Zwaans
##
## Date Created: 2023-09-11
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## -----------------------------


## set working directory for Mac 
setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/TreeDistFinal/")    

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(data.table)
library(ape)
library(TreeDist)
library(gdata)

#########################################################
## Getting likelihood values matching the SciPhy trees ##
#########################################################

##sourcing helper functions
source("~/typewriter_analysis/src/cell_culture/upgma/plot_trees_in_2d_with_likelihoods.R")

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

#############################
## Plotting in RF 2D space ##
############################

#input RF 2D mapping 
RF_treeDist <- read.csv( file = "mapping_df_RF.csv")

#input original distances
distances <- read.csv( file = "distances_RF.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=651, nrow=651) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances
dist <- as.dist(mat)
#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(RF_treeDist)), function(k) {
  newDist <- dist(RF_treeDist[, seq_len(k)])
  MappingQuality(dist, newDist, 10)["TxC"]
}, 0)
plot(txc, xlab = "Dimension")

#plot the mapping quality
png("mapping_quality_RF_MCC.png")
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

tree_df_rf <- data.frame(RF_treeDist)
plot1_2_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_3_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",2,3) + xlab("Coordinate 2") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_4_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",2,4) + xlab("Coordinate 2") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot3_4_RF <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Robinson Foulds",3,4) + xlab("Coordinate 3") + ylab("Coordinate 4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

blank_plot <- ggplot() + theme_void()
combined <- cowplot::plot_grid(plot1_2_RF+ theme(plot.title = element_blank()),plot1_3_RF + theme(plot.title = element_blank()),plot1_4_RF+ theme(plot.title = element_blank()),blank_plot,plot2_3_RF + theme(plot.title = element_blank()),plot2_4_RF + theme(plot.title = element_blank()),blank_plot,blank_plot,plot3_4_RF+ theme(plot.title = element_blank()),nrow = 3,ncol=3)

#plot with titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Robinson Foulds metric space",
    fontface = 'bold',size = 25,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
combined_with_title <- cowplot::plot_grid(
  title, combined,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave("combined_treeDist_4dim_RF_MCC.png", combined_with_title, width = 60, height = 60, units = "cm", dpi = 300)

###########################################
## Plotting in 2D Clustering Information ##
###########################################

CI_treeDist <- read.csv( file = "mapping_df_CI.csv")

#input original distances
distances <- read.csv( file = "distances_CI.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=651, nrow=651) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances
dist <- as.dist(mat)

#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(CI_treeDist)), function(k) {
  newDist <- dist(CI_treeDist[, seq_len(k)])
  MappingQuality(dist, newDist, 10)["TxC"]
}, 0)

#plot the mapping quality
png("mapping_quality_CI.png")
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

#remove the legend, the other plit will have it

tree_df_rf <- CI_treeDist
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",2,3) + xlab("Coordinate 2") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",2,4) + xlab("Coordinate 2") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot3_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",3,4) + xlab("Coordinate 3") + ylab("Coordinate 4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

blank_plot <- ggplot() + theme_void()
combined <- cowplot::plot_grid(plot1_2+ theme(plot.title = element_blank()),plot1_3 + theme(plot.title = element_blank()),plot1_4+ theme(plot.title = element_blank()),blank_plot,plot2_3 + theme(plot.title = element_blank()),plot2_4 + theme(plot.title = element_blank()),blank_plot,blank_plot,plot3_4+ theme(plot.title = element_blank()),nrow = 3,ncol=3)

#plot with titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Clustering Information metric space",
    fontface = 'bold',size = 25,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
combined_with_title <- cowplot::plot_grid(
  title, combined,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave("combined_treeDist_4dim_CI_MCC.png", combined_with_title, width = 60, height = 60, units = "cm", dpi = 300)


###############################################################
## Plotting in 2D Clustering Information and Robinson Foulds ##
###############################################################

combined <- cowplot::plot_grid(plot1_2 + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                               axis.ticks.x=element_blank(), #remove x axis ticks
                                               axis.text.y=element_blank(),  #remove y axis labels
                                               axis.ticks.y=element_blank()),
                               plot1_3 + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                               axis.ticks.x=element_blank(), #remove x axis ticks
                                               axis.text.y=element_blank(),  #remove y axis labels
                                               axis.ticks.y=element_blank())
                               ,plot1_4+ theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                               axis.ticks.x=element_blank(), #remove x axis ticks
                                               axis.text.y=element_blank(),  #remove y axis labels
                                               axis.ticks.y=element_blank()),plot1_2_RF + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                axis.text.y=element_blank(),  #remove y axis labels
                                                                                                axis.ticks.y=element_blank()),plot2_3 + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                              axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                              axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                              axis.ticks.y=element_blank()),plot2_4 + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                                                                            axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                                                                            axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                                                                            axis.ticks.y=element_blank()),plot1_3_RF + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                                                                                                                             axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                                                                                                                             axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                                                                                                                             axis.ticks.y=element_blank()),plot2_3_RF + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                                                                                                                                                                              axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                                                                                                                                                                              axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                                                                                                                                                                              axis.ticks.y=element_blank()),plot3_4 + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                                                                                                                                                                                                                            axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                                                                                                                                                                                                                            axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                                                                                                                                                                                                                            axis.ticks.y=element_blank()),plot1_4_RF + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                                                                                                                                                                                                                                                                             axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                                                                                                                                                                                                                                                                             axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                                                                                                                                                                                                                                                                             axis.ticks.y=element_blank()),plot2_4_RF + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                                                                                                                                                                                                                                                                                                                              axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                                                                                                                                                                                                                                                                                                                              axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                                                                                                                                                                                                                                                                                                                              axis.ticks.y=element_blank()),plot3_4_RF + theme(plot.title = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               axis.ticks.y=element_blank()),nrow = 4,ncol=3)

#adding titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Clustering Information metric space",
    fontface = 'bold',size = 25,
    x = 1,
    hjust = 1
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

title2 <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Robinson Foulds metric space",
    fontface = 'bold',size = 25,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
combined_with_title <- cowplot::plot_grid(
  title, combined,title2,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1,0.05)
)

ggsave("combined_treeDist_4dim_CI_RF_MCC_dual.png", combined_with_title, width = 50, height = 50, units = "cm", dpi = 300)

#############################################
## Plotting in 2D Phylogenetic Information ##
#############################################

PI_treeDist <- read.csv( file = "mapping_df_PI.csv")

#input original distances
distances <- read.csv( file = "distances_PI.csv")

#reformat into lowertriangular matrix
distances <- unlist(distances)
mat <- matrix(NA, ncol=651, nrow=651) 
lowerTriangle(mat, diag=FALSE, byrow=FALSE) <- distances
dist <- as.dist(mat)

#calculate trustworthyness x continuity (mapping quality)
txc <- vapply(seq_len(ncol(PI_treeDist)), function(k) {
  newDist <- dist(PI_treeDist[, seq_len(k)])
  MappingQuality(dist, newDist, 10)["TxC"]
}, 0)

#plot the mapping quality
png("mapping_quality_PI.png")
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)
dev.off()

#remove the legend, the other plit will have it

tree_df_rf <- PI_treeDist
plot1_2 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,2) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(legend.position = "None", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,3) + xlab("Coordinate 1") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot1_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",1,4) + xlab("Coordinate 1") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_3 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",2,3) + xlab("Coordinate 2") + ylab("Coordinate 3") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot2_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",2,4) + xlab("Coordinate 2") + ylab("Coordinate 4") + theme(legend.position = "None",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot3_4 <- plot_tree_distances_topology_metric(tree_df_rf, tree_likelihood, "Clustering Information",3,4) + xlab("Coordinate 3") + ylab("Coordinate 4") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

blank_plot <- ggplot() + theme_void()

blank_plot <- ggplot() + theme_void()
combined <- cowplot::plot_grid(plot1_2+ theme(plot.title = element_blank()),plot1_3 + theme(plot.title = element_blank()),plot1_4+ theme(plot.title = element_blank()),blank_plot,plot2_3 + theme(plot.title = element_blank()),plot2_4 + theme(plot.title = element_blank()),blank_plot,blank_plot,plot3_4+ theme(plot.title = element_blank()),nrow = 3,ncol=3)

#plot with titles
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "Phylogenetic Information metric space",
    fontface = 'bold',size = 25,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
combined_with_title <- cowplot::plot_grid(
  title, combined,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave("combined_treeDist_4dim_PI_MCC.png", combined_with_title, width = 60, height = 60, units = "cm", dpi = 300)

