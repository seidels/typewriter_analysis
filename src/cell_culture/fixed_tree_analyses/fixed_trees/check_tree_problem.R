library(ape)

filename = "~/Projects/typewriter_analysis/results/analysis_cell_culture_data/sensitivity_analysis/fixed_trees/UPGMAtree_1000_scaled_median_heights.txt"
output_tree = "~/Projects/typewriter_analysis/results/analysis_cell_culture_data/sensitivity_analysis/fixed_trees/UPGMAtree_1000_scaled_median_heights_positiveBranchLength.txt"

# read tree
tree = read.tree(filename)

# Problem is that 3 edges are of length 0
sort(tree$edge.length)[1:10]

# find the problematic nodes
zero_branch_nodes = which(tree$edge.length == 0)

# add small value to their branches
small_value = 0.001
tree$edge.length[zero_branch_nodes] = small_value

# write changed tree
write.tree(tree, file = output_tree)
