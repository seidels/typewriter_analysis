## ---------------------------
##
## Script name: plot_ltts
##
## Purpose of script: plot UPGMA, MCC LTTs vs BEAST2 outputs
##
## Author: Antoine Zwaans
##
## Date Created: 2023-07-21
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory to where trees and log files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/single_clock/1000_cells/")    

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)


## ----------------------------------------
## Plot UPGMA and MCC against BEAST2 sample
## ----------------------------------------

##read in SciPhy trees
trees <- ape::read.nexus(file = "smallest_combined.trees")

##read in upgma tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")

##read in the mcc tree
MCC <- ape::read.nexus(file = "MCC_medianHeights.tree")

##read in scaled upgma
upgma_scaled <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/UPGMAtree_1000_medianPosteriorHeight.txt")


#extract ltt coordinates from all trees
all_ltt <- c()
for (i in 1:length(trees)) {
  coords <- ltt.plot.coords(trees[[i]])
  all_ltt <- rbind(all_ltt, cbind(coords,number=rep(i,nrow(coords))))
}

#add upgma and mcc coordinates to the list
ltt_upgma <- ltt.plot.coords(upgma)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma,rep(length(trees)+1,nrow(ltt_upgma))))

ltt_upgma_scaled <- ltt.plot.coords(upgma_scaled)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma_scaled,rep(length(trees)+2,nrow(ltt_upgma_scaled))))

ltt_mcc <- ltt.plot.coords(MCC)                                
all_ltt <- rbind(all_ltt, cbind(ltt_mcc,rep(length(trees)+3,nrow(ltt_mcc))))

all_ltt <- data.frame(all_ltt)

#create another column to colour by type  
all_ltt$type <- c(rep("SciPhy",length(which(all_ltt$number<=(length(trees))))),rep("UPGMA",length(which(all_ltt$number==(length(trees)+1)))),rep("Scaled UPGMA",length(which(all_ltt$number==(length(trees)+2)))),rep("MCC",length(which(all_ltt$number==(length(trees)+3)))))
cols <- c("SciPhy" = "#E1E1F7", "MCC" = "#060647", "UPGMA" = rgb(1,0,0,0.5),"Scaled UPGMA" = "#920B0B")

ltt_all <-  ggplot(all_ltt,aes(x=time,y=N,group=number,colour=type)) + 
            geom_step() + 
            theme_bw() + 
            scale_color_manual(values = cols) + 
            theme(legend.position = c(0.8,0.2),legend.title = element_blank()) +
            ylab("Total lineages") + 
            xlab("Time (days)")

ggsave("LTT_SciPhy_UPGMA_scaledUPGMA_MCC_median_heights.png", ltt_all, width = 30, height = 10, units = "cm", dpi = 300)


