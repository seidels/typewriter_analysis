## ---------------------------
##
## Script name: cut_tree.R
##
## Purpose of script: is to read a nexus file containing a phylogenetic 
## tree into R using the ape package, and then cut the tree at specific
## time points using the cutPhylo function from the RRphylo package. 
## The resulting cutoff trees are then written to separate nexus files 
## for further analysis 
##
## Author: Sophie Seidel
##
## Date Created: 2023-05-25
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##


## set working directory for Mac

setwd("~/Projects/typewriter_analysis/")      

## ---------------------------

# Load the ape and RRphylo packages
library(ape)
library(RRphylo)

# Define the directory
dir = "results/simulation/simulate_tree/"
  
# Read the nexus file
nexus_file <- "simulate_tree.1.trees"
tree <- read.nexus(paste0(dir, nexus_file))

# Define the time points at which to cut the tree
time_points <- c(5, 10, 15, 20)

# Cut the tree at each time point and write to separate nexus files
for (time in time_points) {
  tree_cut <- cutPhylo(tree, age = time)
  
  # Write the new tree to a nexus file
  output_file <- paste0(dir, "tree_1_cut/", "tree_cut_time_", time, ".newick")
  write.tree(tree_cut, file = output_file)
}

