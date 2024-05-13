## ---------------------------
##
## Script name: plot_growth_rates_OU
##
##
## Author: Sophie Seidel
##
## Date Created: 2024-05-13
##
## Copyright (c) Sophie Seidel 2024
## Email: sophie.seidel@posteo.de

## set working directory where the log files are
setwd("~/Projects/typewriter_analysis/results/preliminary_gastruloid/ou_skyline/")

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



timeline_3_changes <- c(0, 4, 7, 8, 11)
timeline_2_changes <- c(0, 4, 7.5, 11)


get_median_and_hpd = function(growth){

  HPD <- HPDinterval(growth)

  name = paste0("growthRate.", 1:ncol(growth))
  median <- as.numeric(sapply( data.frame(growth),median))
  up_bd <- as.numeric(HPD[,"upper"])
  low_bd <- as.numeric(HPD[,"lower"])

  return(data.frame(name=name, median=median, hpd_up = up_bd, hpd_low=low_bd))
}

#load the log file with 3 change points in the estimated rates
log_file <- "3-mGASv2-skyline-ou.10burnin.combined.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:4)] - typewriter_mcmc[,paste0("deathRate.",1:4)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[5, ] =  growth_hpd[4, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_3_changes
growth_hpd$tree = "3 change points"

growth_combined = growth_hpd


#load the log file with 2 change points in the estimated rates
log_file <- "4-mGASv2-skyline-ou.10burnin.combined.log"
typewriter <- read.table(log_file, header = T)
typewriter_mcmc <- as.mcmc(typewriter)
growth <- typewriter_mcmc[,paste0("birthRate.",1:3)] - typewriter_mcmc[,paste0("deathRate.",1:3)]
growth_hpd = get_median_and_hpd(growth)
growth_hpd[4, ] =  growth_hpd[3, ]# add last timepoint, st for time - 11 the same value remains
growth_hpd$t = timeline_2_changes
growth_hpd$tree = "2 change points"
growth_combined = rbind(growth_combined, growth_hpd)


growth_combined

p_growth = ggplot(growth_combined, aes(x=t, y=median, col=tree)) +
  #facet_grid(rows =  (startsWith(x = tree, prefix = "Sci")))+
  #facet_grid(tree~.)+
  #geom_line() +
  geom_step(size=1)+
  geom_stepribbon(aes(ymin = hpd_low, ymax=hpd_up, fill=tree), linetype="dotted", alpha = 0.1)+
  #geom_ribbon(aes(ymin = hpd_low, ymax=hpd_up, fill=tree), linetype="dotted", alpha = 0.1)+
  theme_bw() +
  theme(legend.title=element_blank(), legend.position = c(0.27, 0.858)) +
  ylab(expression("Posterior growth rate [" * d^-1 * "]"))+
  xlab("Time [d]")
p_growth

ggsave("growth_skyline_ou.png", p_growth , width = 15, height = 10, units = "cm", dpi = 500)


