## ---------------------------
##
## Script name: plot_inference_results_UPGMA_BDS
##
## Purpose of script: Generating violin plots from phylodyn estimates on fixed UPGMA/MCC tree
##
## Author: Antoine Zwaans & Sophie Seidel
##
## Date Created: 2023-07-19
##
## Copyright (c) Antoine Zwaans, Sophie Seidel, 2023
## Email: antoine.zwaans@bsse.ethz.ch, sophie.seidel@posteo.de
##

## set working directory where the log files are

setwd("~/Projects/typewriter_analysis/results/analysis_cell_culture_data/fixed_tree_analyses/fixed_trees/")

## load up the packages we will need:
library(reshape2)
library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)


# ---------------------------------------------------------
#  plot the growth rate for the run on the fixed UPGMA tree
# ---------------------------------------------------------

## Add the UPGMA based estimates
typewriter_file <- "1000_UPGMA.1707411001885.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA"

growth_rates = melt(bd_rates, id.vars = "tree")

## Add UPGMA scaled based estimates
typewriter_file <- "1000_UPGMA_medianPosteriorHeight.1707410740300.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA \n Scaled"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add UPGMA scaled + estimated branch lengths based estimates
typewriter_file <- "1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_nodeReheight.combined.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "UPGMA + \n SciPhy Branches"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the MCC tree based estimates
typewriter_file <- "1000_MCC.1707412619101.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "SciPhy MCC"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the main analysis estimates
typewriter_file <- "../../inference_results/clock_per_target/1000_cells/combined.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)
bd_rates <- data.frame(growthRate = typewriter[,"birthRate"] - typewriter[,"deathRate"])
bd_rates$tree = "SciPhy posterior"
bd_rates = melt(bd_rates, id.vars = "tree")

growth_rates = rbind(growth_rates, bd_rates)

## Add the Prior
prior = data.frame(growthRate=rlnorm(n = 1000, meanlog = -0.6, sdlog = 1))
prior$tree = "Prior"
prior =melt(prior, id.vars = "tree")

growth_rates = rbind(growth_rates, prior)

# Define order
growth_rates$tree = factor(growth_rates$tree, levels=c("Prior", "SciPhy posterior", "SciPhy MCC", "UPGMA", "UPGMA \n Scaled", "UPGMA + \n SciPhy Branches"))

p_growth_fixed_tree <- ggplot(growth_rates, aes(x=tree, y=value)) +
  theme_bw() +
  geom_violin(draw_quantiles = c(0.5))+
  ylim(0.3, 0.5)+
  xlab("")+
  ylab("Posterior growth rate per day")



p_growth_fixed_tree
ggsave("growth_rate_fixed_trees.png", p_growth_fixed_tree, width = 15, height = 15, units = "cm", dpi = 300, )

