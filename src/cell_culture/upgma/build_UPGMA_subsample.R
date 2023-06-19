## ---------------------------
##
## Script name: UPGMA subsample
##
## Purpose of script: create upgma from cell culture data subsampled alignments as provided for 
## BEAST2.
##
## Author: Antoine Zwaans
##
## Date Created: 2023-06-14
##
## Copyright (c) Antoine Zwaans, 2023
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
install.packages("ggdendro")
library(ggdendro)
install.packages("dendextend")
library(dendextend)
library(phytools)

# =====================================================================
# Generating the wider cell-by-59EditSites file from cell-by-5EditSites
# =====================================================================


## input set of cells as filtered by S.Seidel
edit_table_filtered = readRDS("data/edit_table_filtered.RDS")

## create a 100 cell sample with seed 1 to match sequence_process_xml.
## TODO use the cell barcodes to ensure data matches instead of "resampling"
## the dataset.

set.seed(1)
cell100 <- sample(unique(edit_table_filtered$Cell),100)

sample_100 <- c()
for(i in cell100) {
  sample_100 <- rbind(sample_100,edit_table_filtered[edit_table_filtered$Cell == i,])
}

sample_100[is.na(sample_100)] <- 'None'
sample_100[sample_100$TargetBC == 'TGGACGAC',7] <- NA
sample_100[sample_100$TargetBC == 'TTTCGTGA',7] <- NA
sample_100[sample_100$TargetBC == 'TGGTTTTG',7] <- NA
sample_100[sample_100$TargetBC == 'TTCACGTA',5:7] <- NA

# edit_cell_table_65 = Ordered cell-by-65EditSites table, including non-existing sites before contraction
edit_cell_table_65 <- select(sample_100, -nUMI) %>%
  pivot_longer(cols = c('Site1','Site2','Site3','Site4','Site5'), names_to = 'Sites', values_to ='Insert') %>%
  pivot_wider(id_cols = Cell, names_from = c(TargetBC,Sites), names_sep = ".", values_from = Insert)
edit_cell_table_65 <- arrange(edit_cell_table_65,Cell) %>%
  select(order(colnames(edit_cell_table_65)))


###########

sub_edit65 <- as.matrix(select(edit_cell_table_65,-Cell))
sub_edit65[is.na(sub_edit65)] <- 'None'
rownames(sub_edit65) <- edit_cell_table_65$Cell

sub_edit59 <- as.matrix(select(edit_cell_table_65,-c('Cell','TGGACGAC.Site5','TTTCGTGA.Site5','TGGTTTTG.Site5',
                                                     'TTCACGTA.Site3','TTCACGTA.Site4','TTCACGTA.Site5')))
rownames(sub_edit59) <- edit_cell_table_65$Cell
cell_list <- edit_cell_table_65$Cell


# =====================================================================
# Generating the phylogenetic tree based on edits
# =====================================================================

# shared_edit_matrix = Counting all shared edits per cell-pair, consistent with the sequential editing on DNA Tape

# Function for calculating shared_edit_matrix
fun_shared_edit_matrix <- function(x) {
  #if (1 == 1){
  sub_edit65 <- x
  sub_edit65[sub_edit65 == 'None'] <- 1:filter(as.data.frame(table(sub_edit65)), sub_edit65 == 'None')$Freq
  ncell <- nrow(sub_edit65)
  cell_list <- sort(rownames(sub_edit65))
  shared_edit_matrix <- matrix(0, ncell,ncell)
  colnames(shared_edit_matrix) <- cell_list
  rownames(shared_edit_matrix) <- cell_list
  for (ii in 1:(ncell)){
    cell1 <- cell_list[ii]
    for (jj in (ii):ncell){
      cell2 <- cell_list[jj]
      for (kk in seq(0,(dim(sub_edit65)[2]-5),5)){
        if (sub_edit65[cell1,(kk+1)] == sub_edit65[cell2,(kk+1)]){
          shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
          if (sub_edit65[cell1,(kk+2)] == sub_edit65[cell2,(kk+2)]){
            shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
            if (sub_edit65[cell1,(kk+3)] == sub_edit65[cell2,(kk+3)]){
              shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
              if (sub_edit65[cell1,(kk+4)] == sub_edit65[cell2,(kk+4)]){
                shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                if (sub_edit65[cell1,(kk+5)] == sub_edit65[cell2,(kk+5)]){
                  shared_edit_matrix[ii,jj] <- shared_edit_matrix[ii,jj] + 1
                }
              }
            }
          }
        }
      }
      shared_edit_matrix[jj,ii] <- shared_edit_matrix[ii,jj]
    }
  }
  return(shared_edit_matrix)
}

shared_edit_matrix <- fun_shared_edit_matrix(sub_edit65) 
shared_edit_matrix <- as.matrix(shared_edit_matrix)
diag(shared_edit_matrix) <- 59

distance_matrix <- 59 - shared_edit_matrix # Phylogenetic distance caludated as (# of possible sites - # of shared sites)
distance_matrix <- as.matrix(distance_matrix)
tree <- as.phylo(hclust(as.dist(distance_matrix), "average")) # tree built using UPGMA



tree_height_calc <- function(tree) { 
  start_edge <- 1
  sum_path <- 0
  while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
    sum_path <- sum_path + tree$edge.length[tree$edge[,2] == start_edge]
    start_edge = tree$edge[tree$edge[,2] == start_edge,1]
    
  }
  return(sum_path)
  
}
tree_100_height <- tree_height_calc(tree)

## we want to scale the tree such that is is smaller than 25 in height.
tree$edge.length <- tree$edge.length * (24.999/tree_100_height)

## save UPGMA as txt file
ape::write.tree(tree, file='results/analysis_cell_culture_data/upgma/UPGMAtree_100_13.txt')
write(cell100,"results/analysis_cell_culture_data/upgma/UPGMAtree_100_13_cell_names.txt")


