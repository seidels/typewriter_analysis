## ---------------------------
##
## Script name: plot mcc_tree_annotated
##
## Purpose of script: plot tree facing the typewriter sequences, with barcode annotations and colors matching insertion probability plots
##
## Date Created: 2023-10-13
##
## ---------------------------
##
## Notes: based on plot_mcc_tree.R
##
## ---------------------------

## set working directory for Mac

setwd("~/Projects/typewriter_analysis")    # personal working directory

## load up the packages we will need:  (uncomment as required)
library(ggtree)
library(ggplot2)
library(treeio)
#library(ape)
## ---------------------------

swap_integer_for_edit = function(integer, insert_to_integer_map){

  insert = insert_to_integer_map[insert_to_integer_map$integer == integer, "insert"]
  return(insert)
}

## ---------------------------

# output dir
output_dir = "results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/"

## define input files

tree_file = "results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/MCC_medianHeights.tree"
alignment = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/alignment_seed1.txt"
cell_ids_file = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/cell_ids_seed1.txt"
edit_file = "results/analysis_cell_culture_data/alignments/simple_1000cells_13tbcs/edit_table_sample_1.csv"
edit_to_integer_map = "data/cell_culture/insert_to_integer_map.csv"

## load input
tree = read.beast(file = tree_file)
t = get.tree(tree)
t$edge.length[which(t$edge.length <0)] = 0
tree_100_cells = treeio::drop.tip(tree, tip = 1:900)
tree_2 = as.treedata(tree = t)

#that s working now!!
sub =  tree_subset(tree = t, node = 1900,  group_node = T, root_edge = T, levels_back = F)

subtree = tidytree::tree_subset(tree = tree_2, node = 1, root_edge = T, levels_back = 3)
ggtree(tree_100_cells) + geom_tiplab()

dat = tree@data
dat[which(dat$node == 1001), "length"] = 0.4
tree@data = dat

insert_to_integer_map = read.csv(edit_to_integer_map)

cell_ids = read.csv(cell_ids_file, header = F)
cell_ids$tip_label = 0:999 # keep numbering of tips as used in alignment
cell_ids = cell_ids[order(cell_ids$V1), ]
cell_ids = cell_ids[which(cell_ids$tip_label > 900), ]


edits = read.csv(edit_file)

## create targetbc matrices
ctr = 1
targetBCs = unique(edits$TargetBC)
targetbc_edits_list =  vector(mode = "list", length = length(targetBCs))

for (targetBC in targetBCs){

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
}

t = get.tree(tree_100_cells)
## plot tree with edit matrix
basic_tree = ggtree(tree_100_cells,layout="ellipse",root.position = 0.4)


for (i in 1:length(targetBCs)){
  print(i)
  labels <- c(targetBCs[i],rep("",ncol(targetbc_edits_list[[i]])-1))
  annotated_basic_tree = gheatmap(basic_tree, targetbc_edits_list[[i]], hjust=0,colnames_offset_y=-10,colnames_offset_x=0,colnames_position = "bottom",width = 0.07, offset = 2.0 * (i-1),custom_column_labels = labels,font.size=2.5) + coord_cartesian(clip="off")
  basic_tree = annotated_basic_tree
}

basic_tree_bars <- basic_tree +
  theme(legend.position = "top",legend.text = element_text(size=22),
        legend.title = element_blank()) + guides(fill = guide_legend(nrow = 1))

label_segs <- data.frame(xstart=seq(from=25+0.05,to=51-0.1,by=26/13),xend=seq(from=25+0.2,to=51-0.1,by=26/13)+1)

#adding small bars and tape labels
for(i in 1:13) {
basic_tree_bars <- basic_tree_bars + geom_segment(x=label_segs$xstart[i],y=-5,xend=label_segs$xend[i],yend=-5,color="black")
}
#adding a makeshift axis
basic_tree_bars_axis <- basic_tree_bars +
  geom_segment(x=0,y=-5,xend=24.8,yend=-5,color="black") +
  geom_text(x=0,y=-15,label="0",size=7) +
  geom_text(x=5,y=-15,label="5",size=7)  +
  geom_text(x=10,y=-15,label="10",size=7)  +
  geom_text(x=15,y=-15,label="15",size=7)  +
  geom_text(x=20,y=-15,label="20",size=7) +
  coord_cartesian(ylim=c(-20,1000),clip="off")

#combined
ggsave(paste0(output_dir,"MCC_with_alignment_annotated.png"), basic_tree_bars_axis, width = 60, height = 40, units = "cm", dpi = 300)
svg(filename = paste0(output_dir,"MCC_focused_clade_100_cells.svg"), width = 60, height = 40, pointsize = 12 )
basic_tree_bars_axis
dev.off()

tree_subset_internal_2 <- function(tree, node, levels_back = 5, root_edge = TRUE) {

  ## error catching to ensure the tree input is of class phylo
  ## if (class(tree) %in% c("phylo", "treedata")) {
  ##   tree_df <- tidytree::as_tibble(tree)
  ## } else {
  ##   stop("tree must be of class 'phylo'")
  ## }

  ## error catching to ensure the levels_back input is numeric
  ## or can be converted to numeric
  if (!is.numeric(levels_back)) {
    levels_back <- as.numeric(levels_back)
    if (is.na(levels_back)) stop("'levels_back' must be of class numeric")
  }

  tree_df <- tidytree::as_tibble(tree)

  selected_node <- node

  is_tip <- tree_df %>%
    dplyr::mutate(isTip = !.data$node %in% .data$parent) %>%
    dplyr::filter(.data$node == selected_node | .data$label == selected_node) %>%
    dplyr::pull(.data$isTip)

  #if (is_tip & levels_back == 0){
  #  stop("The selected node (", selected_node, ") is a tip. 'levels_back' must be > 0",
  #       call. = FALSE)
  #}

  if (is_tip) {
    group_labels <- tree_df %>%
      dplyr::filter(.data$node == selected_node | .data$label == selected_node) %>%
      dplyr::pull(.data$label)
  } else {
    group_labels <- tree_df %>%
      tidytree::offspring(selected_node) %>%
      dplyr::filter(!.data$node %in% .data$parent) %>%
      dplyr::pull(.data$label)
  }

  ## This pipeline returns the tip labels of all nodes related to
  ## the specified node
  ##
  ## The tail/head combo isolates the base node of the subsetted tree
  ## as the output from ancestor lists the closest parent nodes of a
  ## given node from the bototm up.
  ##
  ## It then finds all of the offspring of that parent node. From there
  ## it filters to include only tip and then pulls the labels.

  if (levels_back == 0) {
    new_root_node <- selected_node
  } else {
    new_root_node <- tidytree::ancestor(tree_df, selected_node) %>%
      tail(levels_back) %>%
      head(1) %>%
      dplyr::pull(.data$node)
  }

  subset_labels <- tidytree::offspring(tree_df, new_root_node) %>%
    dplyr::filter(!.data$node %in% .data$parent) %>%
    dplyr::pull(.data$label)

  ## This finds the nodes associated with the labels pulled
  subset_nodes <- which(tree$tip.label %in% subset_labels)

  root.edge <- NULL
  if (is.null(tree$edge.length)) {
    root_edge <- FALSE
    ## if not branch length info, no need to determine root.edge
  }

  if (root_edge) {
    root.edge <- ancestor(tree_df, new_root_node) %>%
      bind_rows(dplyr::filter(tree_df, node == new_root_node)) %>%
      pull(.data$branch.length) %>%
      sum(na.rm = TRUE)
    if (root.edge == 0)
      root.edge <- NULL
  }

  return(list(
    subset_nodes = subset_nodes,
    new_root_node = new_root_node,
    group_labels = group_labels,
    root.edge = root.edge
  ))
}
