

cell_has_selected_barcodes = function(cell, dat, selected_barcodes){

  cell_barcodes = dat[dat$Cell == cell, "TargetBC"]

  if (all(selected_barcodes %in% cell_barcodes)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

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

write_targetBC_alignment_to_xml = function(targets_per_cell_dat, targetBC, filename, barcode, datName){

  write_alignment_start(targetBC, filename)

  write_sequence_into_alignment(targets_per_cell_dat, targetBC, filename)

  write_alignment_closure(filename)
}

write_alignment_start = function(targetBC, filename){

  write(x = c(
    paste0('<data  id="preliminary_gastruloid_targetBC_', targetBC,
           '" spec="Alignment" name="alignment" dataType="integer">')),

    file = filename, ncolumns = 1, append = TRUE, sep = " ")
}


write_sequence_into_alignment = function(targets_per_cell_dat, targetBC, filename){

  dat_with_targetBC = targets_per_cell_dat[which(targets_per_cell_dat$TargetBC == targetBC), ]

  cells = unique(dat_with_targetBC$Cell)

  for (i in 1:length(cells)){

    cell  =  cells[i]
    targets_per_cell = dat_with_targetBC[which(dat_with_targetBC$Cell==cell), ]

    targets_concatenated = paste(
      targets_per_cell[startsWith(colnames(targets_per_cell), prefix = "Site")], collapse = ","
    )

    write(x = c(paste0('<sequence id="', (i-1), '_cell_', cell, '_targetBC_', targetBC,
                       '" spec="Sequence" taxon="',
                       (i-1),'" value="', targets_concatenated, '"/>')),
          file = filename, ncolumns = 1, append = TRUE, sep = " ")
  }
}

write_alignment_closure = function(filename){
  write(x = "</data>", file = filename, append = TRUE, sep = " ")
}

get_nr_of_edits = function(edit_table, unedited_char){

  if (!(is.na(unedited_char))){
    edit_table[edit_table == unedited_char] = NA
  }

  ## compare to data
  nr_of_edits = rep(x = 0, nrow(edit_table))

  for (i in 1:nrow(edit_table)){

    if(is.na(edit_table[i, "Site1"])){

      nr_of_edits[i] = 0
      next

    }else if(is.na(edit_table[i, "Site2"])){

      nr_of_edits[i] = 1
      next

    }else if (is.na(edit_table[i, "Site3"])){

      nr_of_edits[i] = 2
      next

    }else if (is.na(edit_table[i, "Site4"])){

      nr_of_edits[i] = 3
      next

    }else if (is.na(edit_table[i, "Site5"])){

      nr_of_edits[i] = 4
      next

    }else{
      nr_of_edits[i] = 5
    }
  }
  return(nr_of_edits)
}
