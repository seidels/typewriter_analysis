## ---------------------------
##
## Script name: plot_growth_rates
##
## Purpose of script: Plot estimates from SciPhy on HEK293 cell culture data, skyline clock per target, with color annotations matching the tree.
##
## Author: Sophie Seidel
##
## Date Created: 2024-03-10
##
## Copyright (c) Sophie Seidel 2024
## Email: sophie.seidel@posteo.de

## set working directory where the log files are
setwd("~/Projects/typewriter_analysis/results/analysis_cell_culture_data/fixed_tree_analyses/fixed_trees_OU/")

## load up the packages we will need:
library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)
library(pammtools)
library(reshape2)


timeline_format <- c(seq(0,25,by=2), 25)
timeline_format
get_determ_pop_size = function(median){
  total <- 1
  for(i in 1:12) {

    total <- total*exp(median[i]*2)

  }
  total <- total*exp(median[13])

  return(total)
}

get_median_and_hpd = function(growth){

  HPD <- HPDinterval(growth)

  name = paste0("growthRate.", 1:13)
  median <- as.numeric(sapply( data.frame(growth),median))
  up_bd <- as.numeric(HPD[,"upper"])
  low_bd <- as.numeric(HPD[,"lower"])

  return(data.frame(name=name, median=median, hpd_up = up_bd, hpd_low=low_bd))
}

#load the UPGMA log file
typewriter_file <- "1000_UPGMA.combined.log"
typewriter <- read.table(typewriter_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] = c("growthRate.14", growth_hpd[13, 2:4]) # add last timepoint, st for time 24 - 25 the same value remains
get_determ_pop_size(median = growth_hpd$median) # 0.3 million
growth_hpd$tree = "UPGMA"
growth_combined = growth_hpd

#load the UPGMA scaled log file
typewriter_file <- "1000_UPGMA_medianPosteriorHeight.combined.log"
typewriter <- read.table(typewriter_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] = c("growthRate.14", growth_hpd[13, 2:4]) # add last timepoint, st for time 24 - 25 the same value remains
get_determ_pop_size(median = growth_hpd$median) # 0.5 million
growth_hpd$tree = "UPGMA Scaled"
growth_combined = rbind(growth_combined, growth_hpd)

#load the UPGMA branch lengths log file - not converged
typewriter_file <- "1000_UPGMA_medianPosteriorHeight_estimateBranchLengths_nodeReheight.combined.log"
typewriter <- read.table(typewriter_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] = c("growthRate.14", growth_hpd[13, 2:4]) # add last timepoint, st for time 24 - 25 the same value remains
get_determ_pop_size(median = growth_hpd$median)
head(growth_hpd)
growth_hpd$tree = "UPGMA \n SciPhy Branches"
growth_combined = rbind(growth_combined, growth_hpd)

#load the MCC scaled log file
typewriter_file <- "1000_MCC.combined.log"
typewriter <- read.table(typewriter_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] = c("growthRate.14", growth_hpd[13, 2:4]) # add last timepoint, st for time 24 - 25 the same value remains
get_determ_pop_size(median = growth_hpd$median) # 0.7 million
growth_hpd$tree = "SciPhy MCC"
growth_combined = rbind(growth_combined, growth_hpd)

# #load the SciPhy posterior s
typewriter_file <- "../../inference_results/clock_per_target_OU/combined_burnin10.log"
typewriter <- read.table(typewriter_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:13)] - typewriter_mcmc[,paste0("deathRate.",1:13)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[14, ] = c("growthRate.14", growth_hpd[13, 2:4]) # add last timepoint, st for time 24 - 25 the same value remains
get_determ_pop_size(median = growth_hpd$median)
head(growth_hpd)
growth_hpd$tree = "SciPhy posterior"
growth_combined = rbind(growth_combined, growth_hpd)


growth_combined$t = as.numeric(timeline_format)

#growth_long = melt(growth_combined, id.vars = c("tree", "name", "time", "median", "hpd_up", "hpd_low"))

p_growth = ggplot(growth_combined, aes(x=t, y=median, col=tree)) +
  #facet_grid(rows =  (startsWith(x = tree, prefix = "Sci")))+
  #facet_grid(tree~.)+
  #geom_line() +
  geom_step(size=1)+
  geom_stepribbon(aes(ymin = hpd_low, ymax=hpd_up, fill=tree), linetype="dotted", alpha = 0.1)+
  #geom_ribbon(aes(ymin = hpd_low, ymax=hpd_up, fill=tree), linetype="dotted", alpha = 0.1)+
  theme_bw() +
  theme(legend.title=element_blank(), legend.position = c(0.84, 0.858)) +
  ylab(expression("Posterior growth rate [" * d^-1 * "]"))+
  xlab("Time [d]")
p_growth

ggsave("growth_skyline_fixed_trees.png", p_growth , width = 15, height = 10, units = "cm", dpi = 500)


