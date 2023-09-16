## ---------------------------
##
## Script name: UPGMA subsample
##
## Purpose of script: create upgma from cell culture data subsampled alignments as provided for
## BEAST2.
##
## Author: Antoine Zwaans & Sophie Seidel
##
## Date Created: 2023-06-14
##
## Copyright (c) Antoine Zwaans & Sophie Seidel, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Notes: Based on J.Choi's script in DNA typewriter repo.
##
##
## ---------------------------

## set working directory for Mac

setwd("~/typewriter_analysis")

## load up the packages we will need:

require(tidyverse)
require(data.table)
library(ape)
library(phangorn)
library(tidyverse)
library(ggdendro)
library(dendextend)
library(phytools)

source("src/cell_culture/upgma/scripts.R")
source("src/useful_scripts_across_categories.R")

# =====================================================================
# Generating the wider cell-by-59EditSites file from cell-by-5EditSites
# =====================================================================


## input set of cells as filtered by S.Seidel
edit_table_filtered = readRDS("data/cell_culture/edit_table_filtered.RDS")

## create a 100 cell sample with seed 1 to match sequence_process_xml.
## TODO use the cell barcodes to ensure data matches instead of "resampling"
## the dataset.

set.seed(1)
cell100 <- sample(unique(edit_table_filtered$Cell),100)

sample_100 <- c()
for(i in cell100) {
  sample_100 <- rbind(sample_100,edit_table_filtered[edit_table_filtered$Cell == i,])
}

sample_100 = mark_unedited_sites(sample_100)
sample_100 = set_non_contracted_sites_to_na(sample_100)

# nSites is 13 targetBCs * 5 sites = 65 sites; minus the contracted sites,
# i.e. 1x 3 contracted (2xTape) and 3x 1 contracted (4xTape)
tree = build_upgma_tree(edit_table = sample_100, nSites = 59)

tree_100_height <- tree_height_calc(tree)

## we want to scale the tree such that is is smaller than 25 in height.
tree$edge.length <- tree$edge.length * (24.999/tree_100_height)

## save UPGMA as txt file
ape::write.tree(tree, file='results/analysis_cell_culture_data/upgma/UPGMAtree_100_13.txt')
write(cell100,"results/analysis_cell_culture_data/upgma/UPGMAtree_100_13_cell_names.txt")


