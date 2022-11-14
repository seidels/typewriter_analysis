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
