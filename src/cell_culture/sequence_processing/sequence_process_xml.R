## ---------------------------
##
## Script name: sequence_process_xml
##
## Purpose of script: Pre-processing and sampling of sequences for xml inputs
##
## Author: Antoine Zwaans & Sophie Seidel
##
## Date Created: 2022-12-13
##
## Copyright (c) Antoine Zwaans & Sophie Seidel, 2022
## Email: antoine.zwaans@bsse.ethz.ch, sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## set working directory for Mac and PC

setwd("~/Projects/typewriter_analysis/")
library(dplyr)
library(ggplot2)
library(tidyr)
library(boot)
library(cowplot)
library(stringr)

## ---------------------------

output_folder = "results/analysis_cell_culture_data/alignments/"

# load the preprocessed data ---------

#edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))
#load the filtered data:

edit_table_by_5 = readRDS("data/edit_table_filtered.RDS")

### process the dataset to a txt format for beast input: trinucleotides are mapped to an integer
### this dataset is subsampled and saved as txt file as a BEAST Alignment with ScarData


# 1st remove all trinucleotides and convert NAs with NA

for(i in 3:7) {
  edit_table_by_5[,i] <- substring(edit_table_by_5[,i], 1,3)
  edit_table_by_5[which(is.na(edit_table_by_5[,i])),i] <- "NA"

}

# 2nd, map these to an integer:

# get this from the frequencies:
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]))) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$codes <- as.character(seq(0,19))
bulk_insert_count <- data.frame(edit=bulk_insert_count$Var1, code=as.character(seq(0,19)))



# convert these trinucleotides to integer format:
for(i in 3:7) {
  edit_table_by_5[,i] <- unlist(lapply(edit_table_by_5[,i], function(x) {pick_code(x,bulk_insert_count)}))

}

# concatenate them, adding commas:
edit_table_by_5$beast_seq <- apply(edit_table_by_5[,3:7],1,function(x) {str_flatten(x,collapse = ",")})

#extract 100/500/1000 cells from all targetBCs (except the one that only has 2 edits, total 12) from paper:
targetBCs = c("ATGGTAAG","ATTTATAT",
                         "ATTTGGTT", "GCAGGGTG",
                         "GTAAAGAT", "TAGATTTT",
                         "TGCGATTT", "TGGACGAC",
                         "TGGTTTTG", "TTAGATTG",
                         "TTGAGGTG",
                         "TTTCGTGA","TTCACGTA")

#targetBCs identified as having truncation to 4, they will be processed as such in the sampling
#TODO IMPLEMENT THIS DIFFERENTLY
#"TGGACGAC" - number 8
#"TGGTTTTG" - number 9
#"TTTCGTGA" - number 12

#here, SAMPLE WITH 3 SEEDS FOR EACH DATASET SIZE.
#this saves a list of cell IDs in the same order as the taxon numbers in the alignments
for (seed in 1:3) {

  for (n_cells in c(100, 500, 1000)){

    output_folder_ncells = paste0(output_folder, "simple_", n_cells, "cells_13tbcs/")

    cell_sample = sample_dataset_for_BEAST(n_cells, edit_table_by_5, seed)

    edit_table_for_sample = edit_table_by_5[ (edit_table_by_5$Cell %in% cell_sample) & (edit_table_by_5$TargetBC %in% targetBCs), ]
    write.csv(x = edit_table_for_sample, file = paste0(output_folder_ncells, "edit_table_sample.csv"))

    write_cell_ids_to_file(cell_sample = cell_sample, cell_ids_file = paste0(output_folder_ncells, "cell_ids_seed", seed, ".txt"))
    write_alignment_to_xml(cell_sample = cell_sample, dataset = edit_table_by_5, targetBCs = targetBCs, n_cells = n_cells,
                           alignment_file = paste0(output_folder_ncells, "alignment_seed", seed, ".txt"))

  }


}




#copy/paste this from the console for insert frequencies
#bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]),useNA = "always")) %>% arrange(desc(Freq))
#bulk_insert_count <- bulk_insert_count$Freq / sum(bulk_insert_count$Freq)
#cat(bulk_insert_count)
