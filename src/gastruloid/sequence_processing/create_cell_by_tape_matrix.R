## ---------------------------
##
## Purpose of script: Take the input matrix provided by sanjay and convert it into a
## cell by tape matrix
##
## Author: Sophie Seidel
##
## Date Created: 2024-04-12
##
## Copyright (c) Sophie Seidel, 2024
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## set working directory for Mac

setwd("~/Projects/typewriter_analysis/")

## ---------------------------

## load up our functions into memory

source("src/gastruloid/sequence_processing/functions.R")

## ---------------------------

dat = "~/Projects/typewriter_analysis/data/gastruloid/data_by_sanjay/sanjay_matrix.csv"

dat_sanjay = read.csv(file = dat)

### sanjay dat ###

### rows are cells and the rownames give the cell barcode

### columns are individual Tape sites that were edited.
### For example, if the column name is "AAAGTAAATCTC-0-Site1", it means
### AAAGTAAATCTC was the Tape barcode, 0 is the specific locus of integration
### for that barcode, and we're looking at the first of the 6 possible sites
### that could have been edited for that locus.


cell_by_tape_dat = get_cell_by_tape_matrix(dat_sanjay)

saveRDS(cell_by_tape_dat, file = "~/Projects/typewriter_analysis/src/gastruloid/sequence_processing/cell_by_tape.RDS")

