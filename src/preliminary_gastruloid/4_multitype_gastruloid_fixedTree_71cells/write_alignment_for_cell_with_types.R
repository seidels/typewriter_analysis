## ---------------------------
##
## Purpose of script: Only a subset of the cells in our filtered edit table
## comes with cell type annotations. Here, we filter to only keep those cells
## and write the alignment for them
##
## Author: Sophie Seidel
##
## Date Created: 2023-09-16
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
#scripts

source("src/useful_scripts_across_categories.R")

## ---------------------------
# output files
alignment_file = "results/preliminary_gastruloid/multitype_6types_71cells/alignment.xml"
integers_dat_file = "results/preliminary_gastruloid/multitype_6types_71cells/edits_to_integer_map.csv"

## ---------------------------
# input files

## input cell type annotations
cell_types = read.csv("results/preliminary_gastruloid/multitype_6types/annotations.csv",
                      header = F, col.names = c("cell", "type"))


## input set of cells as filtered by S.Seidel
edit_table_filtered = readRDS("data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS")

## ---------------------------
# subset the edit table to those cells that have a type annotation

edit_table_subset = edit_table_filtered[edit_table_filtered$Cell %in% cell_types$cell, ]

## ---------------------------
# write to file

# pre-process
integers_dat = convert_edits_to_integer(edit_table_subset, number_of_sites = 5)

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

# date trait
taxon_and_date_individual = unname(sapply(edit_table_subset$Cell, function(x){
  paste(x, "=10", sep = "")
}))

# to be pasted in xml
taxon_and_date_merged = paste(taxon_and_date_individual, collapse = ",")
write(taxon_and_date_merged, file = "results/preliminary_gastruloid/multitype_6types_71cells/date_traits.xml")



