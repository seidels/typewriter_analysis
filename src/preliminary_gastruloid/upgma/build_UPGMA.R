## ---------------------------
##
## Script name: build UPGMA
##
## Purpose of script: create upgma from gastruloid data
##
## Author: Sophie Seidel
##
## Date Created: 2023-06-14
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes: Based on J.Choi's script in DNA typewriter repo.
##
##
## ---------------------------

## set working directory for Mac

setwd("~/Projects/typewriter_analysis")

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
edit_tab = read.csv("data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_10X_bamExtractV2_t3_collapse.csv")
edit_tab = read.csv("data/preliminary_gastruloid/mGASv2_Lane2_Group1_cell_annotation.csv")

edit_table_filtered = readRDS("data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS")

edit_table_filtered = mark_unedited_sites(edit_table_filtered)

# 8targetBCs times 5 sites
tree = build_upgma_tree(edit_table = edit_table_filtered, nSites = 45)

cell_names = tree$tip.label

tree_height <- tree_height_calc(tree)

## we want to scale the tree such that is is smaller than 25 in height.
tree$edge.length <- tree$edge.length * (9.0/tree_height)
tree$edge.length <- tree$edge.length + 0.01

tree$edge.length

## save UPGMA as txt file
ape::write.tree(tree, file='results/preliminary_gastruloid/upgma/UPGMAtree_780cells_8tBCs.txt')
write(cell_names, "results/preliminary_gastruloid/upgma/UPGMAtree_780cells_8tBCs_cell_names.txt")


