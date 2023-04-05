## ---------------------------
##
## Script name:
##
## Purpose of script:
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

## set working directory for Mac

setwd("~/Projects/typewriter_analysis/")      # Sophie's working directory (mac)

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(ggplot2)
library(reshape2)
## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R")

## ---------------------------
prior = rlnorm(n = nrow(log), meanlog = -2, sdlog = 0.5)
logfile = "results/analysis_cell_culture_data/inference_results/clock_per_target/typewriter_clockPerSite_13Sites_100Cells_DataSet1.combined1-3.thinned1000000.log"
log = read.delim(logfile, header = T)
log$prior_clock = prior
clocks = colnames(log)[startsWith(x = colnames(log), "clock")]

melted = melt(log[, c(56, 42:54)])

ggplot(melted, aes(x=variable, y=value))+
  geom_violin()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95))

