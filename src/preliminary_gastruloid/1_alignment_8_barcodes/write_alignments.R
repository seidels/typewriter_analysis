## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Sophie Seidel
##
## Date Created: 2023-04-19
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

library(plyr)

## load up our functions into memory

source("src/useful_scripts_across_categories.R")

## ---------------------------

# output files
alignment_file = "results/preliminary_gastruloid/1_alignment_8_barcodes/alignment_filtered_for_8barcodes.xml"
integers_dat_file = "results/preliminary_gastruloid/1_alignment_8_barcodes/edits_to_integer_map.csv"

# input file
filtered_dat_file = "data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"
filtered_dat = readRDS(filtered_dat_file)

# pre-process
integers_dat = convert_edits_to_integer(filtered_dat, number_of_sites = 5)

## Remove because Site 6 is unedited throughout
integers_dat = integers_dat[, ! names(integers_dat) %in% c("Site6")]
write.csv(x = integers_dat, file = integers_dat_file)

# write alignment
for (targetBC in unique(integers_dat$TargetBC)){

  print(targetBC)
  write_targetBC_alignment_to_xml(targets_per_cell_dat = integers_dat,
                                  filename = alignment_file,
                                  targetBC = targetBC)
}

