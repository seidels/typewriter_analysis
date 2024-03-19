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

## set working directory to where trees and log_constant files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/")    

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)


## ----------------------------------------
## Plot UPGMA and MCC_constant against BEAST2 sample
## ----------------------------------------

##read in SciPhy trees
trees_constant <- ape::read.nexus(file = "thinned4000000.trees")

##read in upgma tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")

##read in log_constant
log_constant <- read.table("combined.log", header = T)

##read in the MCC_constant tree
MCC_constant <- ape::read.nexus(file = "MCC_on_thinned4000000_check.tree")

#get the median posterior tree height
median_posterior_height_constant <- median(log_constant[,"treeHeight.t.alignment"])

##read in scaled upgma
upgma_scaled_constant <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_scaled_constant$edge.length <- upgma_scaled_constant$edge.length * (median_posterior_height_constant/upgma_height)

#extract ltt coordinates from all trees_constant
all_ltt_constant <- c()
#record at what time point each tree reaches 500 cells and 1000 cells
all_500_times <- c()
all_1000_times <- c()
for (i in 1:length(trees_constant)) {
  coords <- ltt.plot.coords(trees_constant[[i]])
  all_500_times <- c(all_500_times,coords[500,'time'])
  all_1000_times <- c(all_1000_times,coords[1000,'time'])
  all_ltt_constant <- rbind(all_ltt_constant, cbind(coords,number=rep(i,nrow(coords))))
}

#add upgma 
ltt_upgma <- ltt.plot.coords(upgma)
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_upgma,rep(length(trees_constant)+1,nrow(ltt_upgma))))

##add the scaled upgma
ltt_upgma_scaled_constant <- ltt.plot.coords(upgma_scaled_constant)
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_upgma_scaled_constant,rep(length(trees_constant)+2,nrow(ltt_upgma_scaled_constant))))

#MCC_constant
ltt_MCC_constant <- ltt.plot.coords(MCC_constant)   
all_ltt_constant <- rbind(all_ltt_constant, cbind(ltt_MCC_constant,rep(length(trees_constant)+3,nrow(ltt_MCC_constant))))

all_ltt_constant <- data.frame(all_ltt_constant)

#making times start at zero (instead of a negative timescale)
all_ltt_constant$time <- all_ltt_constant$time + 25

#create another column to colour by type  
all_ltt_constant$type <- c(rep("SciPhy constant",length(which(all_ltt_constant$number<=(length(trees_constant))))),rep("UPGMA",length(which(all_ltt_constant$number==(length(trees_constant)+1)))),rep("UPGMA (Scaled constant)",length(which(all_ltt_constant$number==(length(trees_constant)+2)))),rep("SciPhy MCC",length(which(all_ltt_constant$number==(length(trees_constant)+3)))))
cols <- c("SciPhy constant" = "#E1E1F7", "SciPhy MCC" = "black", "UPGMA" = "red","UPGMA (Scaled constant)" = "#E07E5E")

ltt_all <-  ggplot(all_ltt_constant,aes(x=time,y=N,group=number,colour=type)) + 
  geom_step() + 
  theme_bw() + 
  scale_color_manual(values = cols) + 
  theme(legend.position = c(0.8,0.3),legend.title = element_blank(),text = element_text(size = 22),panel.grid.minor = element_blank()) +
  ylab("Total lineages") + 
  xlab("Time (days)")

ggsave("LTT_SciPhy_upgma_scaled_constantUPGMA_MCC_median_heights.png", ltt_all, width = 50, height = 10, units = "cm", dpi = 300)

#when does the LTT upgma reach 500 lineages
25 + ltt_upgma[500,'time']
25 + ltt_upgma_scaled_constant[500,'time']

#when do sciphy trees reach 500
all_500_times <- 25 + all_500_times
all_500_times <- as.mcmc(all_500_times)
HPDinterval(all_500_times)

#when do sciphy trees reach 1000
median(all_1000_times)
all_1000_times <- as.mcmc(all_1000_times)
HPDinterval(all_1000_times)

#loop to examine when LTT mcc and LTT UPGMA overlap/are at certain number of lineage difference
#for(i in 1:length(ltt_mcc[,'time'])) {
#  index_upgma <- which.min(abs(ltt_upgma[,'time'] - ltt_mcc[i,'time']))
#  print(index_upgma)
  
#  if((index_upgma -15) == i) {
#    print("the difference between both LTTs is less than 15")
#    print(25 + ltt_mcc[i,'time'])
#  }
#}

## --------------
## add 13-epochs 
## --------------

## set working directory to where trees and log_constant files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_OU/")    

##read in SciPhy trees
trees <- ape::read.nexus(file = "thinned_5000000.trees")

##read in upgma tree
upgma <- ape::read.tree(file = "~/typewriter_analysis/results/analysis_cell_culture_data/upgma/UPGMAtree_1000.txt")

##read in log
log <- read.table("combined_burnin10.log", header = T)

##read in the mcc tree
MCC <- ape::read.nexus(file = "MCC_median_thinned_5000000.tree")

#get the median posterior tree height
median_posterior_height <- median(log[,"treeHeight.t.alignment"])

##read in scaled upgma
upgma_scaled <- upgma
upgma_height <- tree_height_calc(upgma)
upgma_scaled$edge.length <- upgma_scaled$edge.length * (median_posterior_height/upgma_height)

#extract ltt coordinates from all trees
all_ltt <- c()
#record at what time point each tree reaches 500 cells and 1000 cells
all_500_times <- c()
all_1000_times <- c()
for (i in 1:length(trees)) {
  coords <- ltt.plot.coords(trees[[i]])
  all_500_times <- c(all_500_times,coords[500,'time'])
  all_1000_times <- c(all_1000_times,coords[1000,'time'])
  all_ltt <- rbind(all_ltt, cbind(coords,number=rep(i,nrow(coords))))
}

#add upgma 
ltt_upgma <- ltt.plot.coords(upgma)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma,rep(length(trees)+1,nrow(ltt_upgma))))

##add the scaled upgma
ltt_upgma_scaled <- ltt.plot.coords(upgma_scaled)
all_ltt <- rbind(all_ltt, cbind(ltt_upgma_scaled,rep(length(trees)+2,nrow(ltt_upgma_scaled))))

#MCC
ltt_mcc <- ltt.plot.coords(MCC)   
all_ltt <- rbind(all_ltt, cbind(ltt_mcc,rep(length(trees)+3,nrow(ltt_mcc))))

all_ltt <- data.frame(all_ltt)

#making times start at zero (instead of a negative timescale)
all_ltt$time <- all_ltt$time + 25


all_ltt$type <- c(rep("SciPhy 13-epochs",length(which(all_ltt$number<=(length(trees))))),rep("UPGMA",length(which(all_ltt$number==(length(trees)+1)))),rep("UPGMA (Scaled 13-epochs)",length(which(all_ltt$number==(length(trees)+2)))),rep("SciPhy MCC",length(which(all_ltt$number==(length(trees)+3)))))
cols <- c("SciPhy constant" = "#E1E1F7", "UPGMA" = "red","UPGMA (Scaled constant)" = "#E07E5E","SciPhy MCC" = "black", "SciPhy 13-epochs" = "#56996E", "UPGMA" = "red","UPGMA (Scaled 13-epochs)" = "#E07E5E")
alphas <- c("SciPhy constant" = 0.3, "UPGMA" = 1.0,"UPGMA (Scaled constant)" = 1.0,"SciPhy MCC" = 1.0,"SciPhy 13-epochs" = 0.3, "UPGMA" = 1.0,"UPGMA (Scaled 13-epochs)" = 1.0)

ltt_all <-  ggplot(all_ltt,aes(x=time,y=N,group=number,colour=type,alpha = type)) + 
  geom_step() +
  theme_bw() + 
  scale_color_manual(values = cols) + 
  theme(legend.position = c(0.8,0.3),legend.title = element_blank(),text = element_text(size = 22),panel.grid.minor = element_blank()) +
  ylab("Total lineages") + 
  xlab("Time (days)") + scale_alpha_manual(values=alphas)

ggsave("LTT_SciPhy_UPGMA_scaledUPGMA_MCC_median_heights_OU.png", ltt_all, width = 50, height = 10, units = "cm", dpi = 300)


ltt_all <-  ggplot(all_ltt_constant,aes(x=time,y=N,group=number,colour=type,alpha = type)) + 
  geom_step() + 
  geom_step(data=all_ltt) +
  theme_bw() + 
  scale_color_manual(values = cols) + 
  theme(legend.position = c(0.8,0.3),legend.title = element_blank(),text = element_text(size = 22),panel.grid.minor = element_blank()) +
  ylab("Total lineages") + 
  xlab("Time (days)") + scale_alpha_manual(values=alphas)

ggsave("LTT_SciPhy_UPGMA_scaledUPGMA_MCC_median_heights_constant_vs_OU.png", ltt_all, width = 50, height = 10, units = "cm", dpi = 300)
