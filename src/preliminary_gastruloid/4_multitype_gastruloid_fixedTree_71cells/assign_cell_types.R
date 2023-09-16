## ---------------------------
##
## Script name: assign cell types
##
## Purpose of script: Assign cell types to
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
source("src/preliminary_gastruloid/scripts_across_gastruloid.R")

## ---------------------------
# input cell type annotations
cell_types = read.csv("results/preliminary_gastruloid/multitype_6types_71cells/annotations.csv",
                      header = F, col.names = c("cell", "type"))
cell_types$type_number = sapply(cell_types$type, assign_type_to_number)

# input cell edit file
integers_dat_file = "results/preliminary_gastruloid/multitype_6types_71cells/edits_to_integer_map.csv"
integers_dat = read.csv(integers_dat_file)

## ---------------------------
# subset cell types to match cells in alignment data; and order according to cell in alignment data
cell_types = cell_types[cell_types$cell %in% integers_dat$Cell, ]
cell_types = cell_types[match(unique(integers_dat$Cell), cell_types$cell), ]

## ---------------------------
# write cell type traits to be inserted into the BEAST 2 xml
# date trait
taxon_and_type_individual = unname(apply(cell_types, MARGIN = 1, function(x){
  paste(x[1], x[3], sep = "=")
}))

# to be pasted in xml
taxon_and_type_merged = paste(taxon_and_type_individual, collapse = ",")
write(taxon_and_type_merged, file = "results/preliminary_gastruloid/multitype_6types_71cells/type_traits.xml")


