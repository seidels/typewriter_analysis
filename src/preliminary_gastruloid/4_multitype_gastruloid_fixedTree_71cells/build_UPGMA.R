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


# input cell type annotations
cell_types = read.csv("results/preliminary_gastruloid/multitype_6types/annotations.csv",
                      header = F, col.names = c("cell", "type"))


## input set of cells as filtered by S.Seidel
edit_table_filtered = readRDS("data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS")
edit_table_filtered = mark_unedited_sites(edit_table_filtered)

# subset the edit table to those cells that have a type annotation
edit_table_subset = edit_table_filtered[edit_table_filtered$Cell %in% cell_types$cell, ]

# assign cell types to numbers
cell_types$type_number = sapply(cell_types$type, assign_type_to_number)

# 8targetBCs times 5 sites
tree = build_upgma_tree(edit_table = edit_table_subset, nSites = 45)

cell_names = tree$tip.label
cell_number = length(tree$tip.label)

tree_height <- tree_height_calc(tree)

## we want to scale the tree such that is is smaller than 25 in height.
tree$edge.length <- tree$edge.length * (9.0/tree_height)
tree$edge.length <- tree$edge.length + 0.001

sort(tree$edge.length)

## save UPGMA as txt file
ape::write.tree(tree, file= paste0('results/preliminary_gastruloid/upgma/UPGMAtree_', cell_number, 'cells_8tBCs_funnyBranchLengths_forBEASTinit.txt'))
write(cell_names, paste0("results/preliminary_gastruloid/upgma/UPGMAtree_", cell_number, "cells_8tBCs_cell_names.txt"))


assign_type_to_number = function(type){

  if (type == "Prog"){
    return(0)
  }else if (type == " SMD"){
    return(1)
  }else if (type == " NMP"){
    return(2)
  }else if (type == " PhMD"){
    return(3)
  }else if(type == " PxMD"){
    return(4)
  }else if(type == " SC"){
    return(5)
  }else{
    stop("Type not found!")
  }
}
