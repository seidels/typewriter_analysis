## ---------------------------
##
## Script name: write alignments
##
## Purpose of script: write alignments given targets per cell input data
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

setwd("~/typewriter_analysis/")      # Sophie's working directory (mac)

## ---------------------------

library(plyr)

## load up our functions into memory

source("src/useful_scripts_across_categories_with_annotations.R")

## ---------------------------

# output files
alignment_file = "alignment_filtered_for_8barcodes_100.xml"
integers_dat_file = "edits_to_integer_map_100.csv"
conversion_table_file = "conversion_table_100.csv"

# input file
filtered_dat_file = "gastruloid_subsample.RDS"
filtered_dat = readRDS(filtered_dat_file)

#subsample 100 cells:
set.seed(1)
filtered_100_cells <- sample(unique(filtered_dat$Cell),100)
filtered_100 <- c()
for(i in filtered_100_cells) {
  filtered_100 <- rbind(filtered_100,filtered_dat[filtered_dat$Cell==i,])
} 

# pre-process
conversion <- convert_edits_to_integer(filtered_100, number_of_sites = 5)
# get the converted sequences
integers_dat <- conversion$targets_per_cell
# get the mapping of trinucleotide to integer
conversion_table <- conversion$possible_edits

# save this mapping: 
write.csv(x = conversion_table, file = conversion_table_file)

# save the converted sequences
# remove because Site 6 is unedited throughout
integers_dat = integers_dat[, ! names(integers_dat) %in% c("Site6")]
unique(unlist(integers_dat[,3:7]))
write.csv(x = integers_dat, file = integers_dat_file)

# write alignment
for (targetBC in unique(integers_dat$TargetBC)){

  print(targetBC)
  write_targetBC_alignment_to_xml(targets_per_cell_dat = integers_dat,
                                  filename = alignment_file,
                                  targetBC = targetBC)
}

integers_dat
