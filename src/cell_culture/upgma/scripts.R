mark_unedited_sites = function(edit_table){
  # set unedited sites to 'None' , s.t. this is treated as a shared site in upgma
  edit_table[is.na(edit_table)] <- 'None'

  return(edit_table)
}

set_non_contracted_sites_to_na = function(edit_table){

  edit_table[edit_table$TargetBC == 'TGGACGAC',7] <- NA
  edit_table[edit_table$TargetBC == 'TTTCGTGA',7] <- NA
  edit_table[edit_table$TargetBC == 'TGGTTTTG',7] <- NA
  edit_table[edit_table$TargetBC == 'TTCACGTA',5:7] <- NA

  return(edit_table)
}


tree_height_calc <- function(tree) {
  start_edge <- 1
  sum_path <- 0
  while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
    sum_path <- sum_path + tree$edge.length[tree$edge[,2] == start_edge]
    start_edge = tree$edge[tree$edge[,2] == start_edge,1]

  }
  return(sum_path)

}
