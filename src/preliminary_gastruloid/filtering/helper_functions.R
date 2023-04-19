## ---------------------------
##
## Script name: helper_functions
##
## Purpose of script: Collect functions to be used by filtering
## scripts in this directory
##
## Author: Sophie Seidel
##
## Date Created: 2023-04-19
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## ---------------------------


cell_has_selected_barcodes = function(cell, dat, selected_barcodes){

  cell_barcodes = dat[dat$Cell == cell, "TargetBC"]

  if (all(selected_barcodes %in% cell_barcodes)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

