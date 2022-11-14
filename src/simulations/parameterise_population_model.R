## ---------------------------
##
## Script name: parameterise_population_model
##
## Purpose of script: Find and document the population model parameters we will use in the simulation.
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

setwd("~/Projects/typewriter_analysis/src/simulations/")      # Sophie's working directory (mac)

## ---------------------------

# Find the population growth rate
## From the type writer paper:
## ... to isolate a monoclonal cell line that grew from 1 cell to ~1.2 million cells via ~20 doublings over 25 days
N_t = 1200000
N_0 = 1
t = 25 #days

## assume N(t) = N(0) * exp(r*t) and rearrange to compute the growth rate r
growth_rate = log(N_t / N_0) / t
growth_rate

# Find the cell division rate (birth rate)
## From the type writer paper:
## ... to isolate a monoclonal cell line that grew from 1 cell to ~1.2 million cells via ~20 doublings over 25 days
number_of_divisions = 20
cell_division_rate = number_of_divisions / t
round(cell_division_rate, digits = 1)

# Find the apoptosis rate (death rate)
## Assume death rate is the difference between the population growth rate and the cell division rate
cell_death_rate = cell_division_rate -growth_rate
round(cell_death_rate, digits = 1)

# Find the sampling proportion
## From the type writer paper:
## Applying these filters left 3,257 cells, for each of which we recovered intact TAPE-1 sequences...
N_sampled_cells = 3257
sampling_proportion = N_sampled_cells / N_t
round(sampling_proportion, 3)
