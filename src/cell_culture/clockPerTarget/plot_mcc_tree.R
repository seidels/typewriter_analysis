## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Sophie Seidel
##
## Date Created: 2023-08-09
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

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(ggtree)
library(treeio)

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R")

swap_integer_for_edit = function(integer, insert_to_integer_map){

  insert = insert_to_integer_map[insert_to_integer_map$integer == integer, "insert"]
  return(insert)
}

## ---------------------------

# output dir
output_dir = "results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/"

## define input files

tree_file = "./results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/MCC_medianHeights.tree"
alignment = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/alignment_seed1.txt"
cell_ids_file = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/cell_ids_seed1.txt"
edit_file = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/edit_table_sample_1.csv"
edit_to_integer_map = "data/cell_culture/insert_to_integer_map.csv"
colour_file = "src/cell_culture/clockPerTarget/edit_color_legend.csv"


## load input
tree = read.beast(tree_file)
insert_to_integer_map = read.csv(edit_to_integer_map)

cell_ids = read.csv(cell_ids_file, header = F)
cell_ids$tip_label = 0:999 # keep numbering of tips as used in alignment
cell_ids = cell_ids[order(cell_ids$V1), ]


edits = read.csv(edit_file)

edit_colours = read.csv(colour_file)
colour_scale = edit_colours$cols
names(colour_scale) = substring(edit_colours$label, first = 1, last = 3)
names(colour_scale)[1] = NA
colour_scale

# plot the tree
p_new = ggtree(tree) +
  theme_tree2()

## add subplots with targetbc matrices
ctr = 1
targetBCs = unique(edits$TargetBC)
targetbc_edits_list =  vector(mode = "list", length = length(targetBCs))

for (targetBC in targetBCs){

  p = p_new
  print(targetBC)
  edits_tbc = edits[ edits$TargetBC == targetBC, ]

  # add cell ids ;
  ## order cell barcodes alphabetically and the reuse the tip labels from cell ids
  edits_tbc = edits_tbc[order(edits_tbc$Cell), ]
  rownames(edits_tbc) = cell_ids$tip_label
  edits_tbc = edits_tbc[ , 4:8]

  # convert edits from integers to trinucleotides
  edits_tbc_trinucl = data.frame(Site1 = sapply(edits_tbc$Site1, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site2 = sapply(edits_tbc$Site2, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site3 = sapply(edits_tbc$Site3, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site4 = sapply(edits_tbc$Site4, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
                             Site5 = sapply(edits_tbc$Site5, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)})
  )
  rownames(edits_tbc_trinucl) = rownames(edits_tbc)

  targetbc_edits_list[[ctr]] = edits_tbc_trinucl
  ctr = ctr + 1
  #p_new = gheatmap(p, edits_trinucl, width = 0.3, offset = ctr * 15) + theme(legend.position = "top")
}
# p_new
#
# ### Initially, do this just for a single target bc
# edits = edits[edits$TargetBC == "TGGTTTTG", ]
#
# edits = edits[edits$TargetBC == "ATTTGGTT", ]
#
# ## match ids from tree and edit data via the cell ids
# edits = edits[order(edits$Cell), ]
# rownames(edits) = cell_ids$tip_label
# edits = edits[, 4:8]
#
#
# ## convert edits from integers to trinucleotides
# edits_trinucl = data.frame(Site1 = sapply(edits$Site1, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
#                            Site2 = sapply(edits$Site2, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
#                            Site3 = sapply(edits$Site3, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
#                            Site4 = sapply(edits$Site4, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)}),
#                            Site5 = sapply(edits$Site5, function(x) {swap_integer_for_edit(integer = x, insert_to_integer_map = insert_to_integer_map)})
#                            )
# rownames(edits_trinucl) = rownames(edits)
#
# ## assign edit to tree


## plot tree with edit matrix
p_old = ggtree(tree)+  theme_tree2() + xlim_tree(xlim = 25)
p_old

for (i in 1:length(targetBCs)){

  p_new = gheatmap(p_old, targetbc_edits_list[[i]], width = 0.07, colnames = F, offset = 2.0 * (i-1)) + scale_fill_manual(values = colour_scale)
  p_old = p_new
}

p_final = p_new + theme(legend.position = "top")# + guides(fill=guide_legend("insertBC", nrow = 1))
p_final

svg(filename = paste0(output_dir, "mcc_1_with_alignment.svg"), width = 15, height = 10 )
p_final
dev.off()

ggsave(filename = paste0(output_dir, "mcc_1_with_alignment.pdf"), plot = p_final, width = 15, height = 10
        )
