## ---------------------------
##
## Script name: sequence_process_xml
##
## Purpose of script: Pre-processing and sampling of sequences for xml inputs
##
## Author: Sophie Seidel
##
## Date Created: 2022-12-13
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

## set working directory for Mac and PC

setwd("~/Projects/typewriter_analysis/")

## ---------------------------
# load libs

library(stringr)
library(dplyr)

## ---------------------------
source("src/cell_culture/sequence_processing/functions.R")

## ---------------------------

output_folder = "results/analysis_cell_culture_data/alignments/alignments_2sites_1000cells_13tbcs/"

#load the filtered data:

edit_table_by_5 = readRDS("data/cell_culture/edit_table_filtered.RDS")
## remove tapes 3:5
edit_table_by_5 = edit_table_by_5[, c(1:4,8)]
site_columns = 3:4

# 1st Convert 6-mers to 3-mers
# by removing the last 3 nucleotides that do not contain information
# (they are always GGA); replace <Na> by "NA"

for(site_column in site_columns) {

  edit_table_by_5[, site_column] <- substring(edit_table_by_5[, site_column], 1,3)
  edit_table_by_5[which(is.na(edit_table_by_5[, site_column])), site_column] <- "NA"

}

# 2nd, create a map from the trinucleotides to an integer
insert_to_integer_map <- data.frame(table(unlist(edit_table_by_5[, site_columns])))
insert_to_integer_map = insert_to_integer_map[order(insert_to_integer_map$Var1), ]
insert_to_integer_map = insert_to_integer_map[which(insert_to_integer_map$Var1 != "NA"), ]
insert_to_integer_map = rbind(c("NA", NA), insert_to_integer_map)
insert_to_integer_map$integer = 0: (nrow(insert_to_integer_map)-1)
colnames(insert_to_integer_map)[1] = "insert"

#insert_to_integer_map = get_insert_to_integer_map(edit_table_by_5, columns = site_columns) ## Not used because no stable ordering
write.csv(insert_to_integer_map, file = "data/cell_culture/insert_to_integer_map_2_sites.csv")

# convert these trinucleotides to integer in the edit table
for(site_column in site_columns) {
  edit_table_by_5[, site_column] <- unlist(lapply(edit_table_by_5[, site_column], function(x) {pick_code(x, insert_to_integer_map)}))
}

# concatenate edits at sites, adding commas:
edit_table_by_5$beast_seq <- apply(edit_table_by_5[, site_columns], 1, function(x) {str_flatten(x,collapse = ",")})

#targetBCs as used in the paper:
targetBCs = c("ATGGTAAG","ATTTATAT",
                         "ATTTGGTT", "GCAGGGTG",
                         "GTAAAGAT", "TAGATTTT",
                         "TGCGATTT", "TGGACGAC",
                         "TGGTTTTG", "TTAGATTG",
                         "TTGAGGTG",
                         "TTTCGTGA","TTCACGTA")

# create alignments for different cell numbers, 3 different sub samples each
n_cells = 1000

for (seed in 1:3) {

    print(paste("Getting alignments for seed ", seed))

    #output_folder_ncells = paste0(output_folder, "two_sites_", n_cells, "cells_13tbcs/")

    cell_sample = subsample_dataset(n_cells, edit_table_by_5, seed)

    edit_table_for_sample = edit_table_by_5[ (edit_table_by_5$Cell %in% cell_sample) & (edit_table_by_5$TargetBC %in% targetBCs), ]
    write.csv(x = edit_table_for_sample, file = paste0(output_folder, "edit_table_sample_", seed, ".csv"))

    write_cell_ids_to_file(cell_sample = cell_sample, cell_ids_file = paste0(output_folder, "cell_ids_seed", seed, ".txt"))
    write_alignment_to_xml(cell_sample = cell_sample, dataset = edit_table_by_5, targetBCs = targetBCs, n_cells = n_cells,
                           alignment_file = paste0(output_folder, "alignment_seed", seed, ".txt"), targetbcs_same_length = T)
}

## Checks
edit_table_by_5[which(edit_table_by_5$TargetBC == "TTCACGTA"), ]

