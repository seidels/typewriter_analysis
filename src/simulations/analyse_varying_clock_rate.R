## ---------------------------
##
## Script name: analyse_varying_clock_rate.R
##
## Purpose of script: Analyse the data generated by simulating with
## different clock rates.
##
## Author: Sophie Seidel
##
## Date Created: 2022-11-14
##
## Copyright (c) Sophie Seidel, 2022
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## set working directory for Mac

setwd("~/Projects/typewriter_analysis/results/simulation/vary_clock_rates/")      # Sophie's working directory (mac)

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
library(TreeTools)
library(ggplot2)
## ---------------------------

## load up our functions into memory
source("../../../src/useful_scripts_across_categories.R")

## ---------------------------

## set up plot settings
plot_dir = "../../../src/simulations/plots/vary_clock_rates/"

clock_rates = c(0.01, 0.1, 0.3)

seed = 1
for (clock_rate in clock_rates){

  # read and process alignment
  alignment_file = paste0("simulation.clockRate=", clock_rate, ".", seed, ".alignment.nexus")
  dat = ReadCharacters(filepath = alignment_file)
  dat = as.data.frame(dat[ ,c(1,3,5,7,9)]) #remove columns that only contain commas
  colnames(dat) = c("Site1", "Site2", "Site3", "Site4", "Site5")
  dat$nrOfEdits = get_nr_of_edits(dat, unedited_char = "0")

  # plot number of edits
  g = ggplot(dat, aes(x=nrOfEdits)) +
    geom_histogram()+
    theme_minimal()+
    ylim(0, 3700)
  ggsave(filename = paste0(plot_dir, "number_of_edits_clock_", clock_rate, "_", seed, ".pdf"))
}

seed = 4
for (clock_rate in clock_rates){

  # read and process alignment
  alignment_file = paste0("simulation.clockRate=", clock_rate, ".", seed, ".alignment.nexus")
  dat = ReadCharacters(filepath = alignment_file)
  dat = as.data.frame(dat[ ,c(1,3,5,7,9)]) #remove columns that only contain commas
  colnames(dat) = c("Site1", "Site2", "Site3", "Site4", "Site5")
  dat$nrOfEdits = get_nr_of_edits(dat, unedited_char = "0")

  # plot number of edits
  g = ggplot(dat, aes(x=nrOfEdits)) +
    geom_histogram()+
    theme_minimal()+
    ylim(0, 3700)
  ggsave(filename = paste0(plot_dir, "number_of_edits_clock_", clock_rate, "_", seed, ".pdf"))
}


