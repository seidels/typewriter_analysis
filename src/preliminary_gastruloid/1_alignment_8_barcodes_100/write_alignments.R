## ---------------------------
##
## Script name: write alignments and conversion table
##
## Purpose of script: subsample and write alignments given targets per cell input data filtered with edit conversion (adapted from write_alignments)
##

## set working directory for Mac

setwd("~/typewriter_analysis/")      # typewriter working director

## ---------------------------

library(plyr)

## load up our functions into memory

source("src/useful_scripts_across_categories.R")

## ---------------------------

# output files
alignment_file = "results/preliminary_gastruloid/1_alignment_8_barcodes_100/alignment_filtered_for_8barcodes_100.xml"
integers_dat_file = "results/preliminary_gastruloid/1_alignment_8_barcodes_100/edits_to_integer_map_100.csv"
conversion_table_file = "results/preliminary_gastruloid/1_alignment_8_barcodes_100/conversion_table_100.csv"

# input file
filtered_dat_file = "data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"
filtered_dat = readRDS(filtered_dat_file)

#subsample 100 cells:
set.seed(1)
filtered_100_cells <- sample(unique(filtered_dat$Cell),100)
filtered_100 <- c()
for(i in filtered_100_cells) {
  filtered_100 <- rbind(filtered_100,filtered_dat[filtered_dat$Cell==i,])
} 

# pre-process
conversion <- convert_edits_to_integer_with_edit_list(filtered_100, number_of_sites = 5)
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

  write_targetBC_alignment_to_xml(targets_per_cell_dat = integers_dat,
                                  filename = alignment_file,
                                  targetBC = targetBC)
}

