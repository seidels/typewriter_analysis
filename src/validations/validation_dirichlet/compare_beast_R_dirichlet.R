## ---------------------------
##
## Script name: validate_dirichlet
##
## Purpose of script: Check that the r dirichlet and the beast dirichlet
## distributions are equivalently parameterised and that they are sampling
## from the same distribution
##
## Author: Sophie Seidel
##
## Date Created: 2023-01-10
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

setwd("~/Projects/typewriter_analysis/src/validation_dirichlet/")      # Sophie's working directory (mac)

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
library(ggplot2)

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)

## ---------------------------
beast_log = "../../results/validation_dirichlet/sample_from_dirichlet.1.log"
samples_from_beast = read.delim(file = beast_log )
n_samples = nrow(samples_from_beast)
samples_from_beast = samples_from_beast[(0.1 * n_samples) : n_samples, ]

samples_from_R = as.data.frame(extraDistr::rdirichlet(n=9001, alpha = rep(1.5,13)))
colnames(samples_from_R) = paste0("insertRates.", 1:13)

median(samples_from_beast$insertRates.1)
median(samples_from_R[, 1])

summary(samples_from_beast$insertRates.1)
summary(samples_from_R[, 1])

ggplot(samples_from_beast, aes(x=insertRates.6)) +
  geom_density(col="darkblue")+
  geom_density(data = samples_from_R, aes(x=insertRates.1), col="darkgreen")+
  theme_classic()

