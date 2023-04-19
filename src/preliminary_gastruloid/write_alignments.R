## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Sophie Seidel
##
## Date Created: 2023-04-19
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

library(plyr)

## load up our functions into memory

# source("functions/summarise_data.R")

## ---------------------------

alignment_file = "results/preliminary_gastruloid/alignments/1_alignment_filtered_for_8barcodes.xml"
edit_to_integer_file = "results/preliminary_gastruloid/alignments/1_edits_to_integer_map.csv"

filtered_dat_file = "data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_filtered_for_8barcodes.RDS"
filtered_dat = readRDS(filtered_dat_file)

edits_as_integers_dat = convert_edits_to_integer(filtered_dat, number_of_sites = 5)
write.csv(x = edits_as_integers_dat, file = edit_to_integer_file)

write_xml(edits_as_integers_dat, output_file)



convert_edits_to_integer = function(targets_per_cell, number_of_sites=5){

  possible_edits = get_possible_edits(targets_per_cell)

  for (site_i in 1:number_of_sites){

    site_string = paste0("Site", site_i)
    targets_per_cell[, site_string] = mapvalues(x = targets_per_cell[, site_string], from = possible_edits$edits, to = possible_edits$integer)
  }

  return(targets_per_cell)
}

get_possible_edits = function(targets_per_cell){

  columns_with_edits = colnames(targets_per_cell)
  columns_with_edits = columns_with_edits[startsWith(x = columns_with_edits, prefix = "Site")]

  all_edits = unique(unlist(filtered_dat[ , columns_with_edits]))
  all_edits = all_edits[all_edits != "None"]

  edit_to_integer_map = data.frame(edits = all_edits, integer=1:length(all_edits))
  edit_to_integer_map[nrow(edit_to_integer_map) + 1, ] = c("None", 0)

  return(edit_to_integer_map)
}

write_site_alignment_to_xml = function(targets_per_cell, filename, barcode, datName){

  cells = targets_per_cell$Cell

  for (i in 1:length(cells)){
    cell  =  cells[i] #alignment[i, "cell"]
    cellState = alignment[i, "beastState"][[1]][site]

    write(x = c(paste0('<sequence id="', datName, '_cell_', cellNr, '_site', site, '" spec="Sequence" taxon="',
                       i,'" value="', cellState, ',"/>')),
          file = filename, ncolumns = 1, append = TRUE, sep = " ")
  }
}

append_site_alignment = function(alignment, filename, datName, nStates=3){

  for (site in 1:10){
    # write cluster alignment start
    write(x = c(paste0('<data  id="', datName, '_site', site, '" spec="Alignment" name="alignment">
                       <userDataType spec="beast.evolution.datatype.ScarData"
                     nrOfStates="', nStates, '" />')),
          file = filename, ncolumns = 1, append = TRUE, sep = " ")

    #write cluster alignment body
    write_site_alignment_to_xml(alignment, filename, site, datName)

    # write alignment end
    write(x = "</data>", file = filename, append = TRUE, sep = " ")
  }

}

