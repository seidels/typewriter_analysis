pick_code <- function(edit,code_map) {
  code <- code_map[which(code_map$edit==edit),2]
  return(code)
  
}

sample_dataset_for_BEAST <- function(n_cells, data, seednr=1) {
  
  set.seed(seednr)
  
  cell_sample <- sample(unique(dataset$Cell), n_cells)
  
  return(cell_sample)
}

write_cell_ids_to_file = function(cell_sample, cell_ids_file){
  
  for(i in 1:length(cell_sample)) {
    write(x = cell_sample[i], file = cell_ids_file, append=TRUE)
  }
}


write_alignment_to_xml = function(cell_sample, dataset, targetBCs, n_cells, alignment_file){
  
  for(i in 1:length(targetBCs)) {
    
    targetBC <- targetBCs[i]
    sampled <- collect_sequences_per_cell_targetBC(cell_sample, dataset=dataset, targetBC)
    
    if((targetBC == "TGGACGAC") | (targetBC == "TGGTTTTG") | (targetBC == "TTTCGTGA" )) {
      
      write( paste0("<data  id=\"data_",targetBC,"\" spec=\"Alignment\" name=\"alignment\" >
            <userDataType spec=\"beast.evolution.datatype.ScarData\" nrOfStates=\"20\"/>"),alignment_file, append = TRUE)
      
      for(j in 1:n_cells) {
        write(paste0("            <sequence spec=\"Sequence\" taxon=\"",j-1,"\"  value=\"",str_sub(sampled$beast_seq[j],end =-3),"\"/>"),
              alignment_file,append=TRUE)
        
      }
    } else if (targetBC == "TTCACGTA") {
      write( paste0("<data  id=\"data_",targetBC,"\" spec=\"Alignment\" name=\"alignment\" >
            <userDataType spec=\"beast.evolution.datatype.ScarData\" nrOfStates=\"20\"/>"),alignment_file,append = TRUE)
      for(j in 1:n_cells) {
        write(paste0("            <sequence spec=\"Sequence\" taxon=\"",j-1,"\"  value=\"",str_sub(sampled$beast_seq[j],end =-7),"\"/>"),alignment_file,append=TRUE)
        
      }
    } else{
      write( paste0("<data  id=\"data_",targetBC,"\" spec=\"Alignment\" name=\"alignment\" >
            <userDataType spec=\"beast.evolution.datatype.ScarData\" nrOfStates=\"20\"/>"),alignment_file,append = TRUE)
      for(j in 1:n_cells) {
        write(paste0("            <sequence spec=\"Sequence\" taxon=\"",j-1,"\"  value=\"",sampled$beast_seq[j],"\"/>"),alignment_file,append=TRUE)
        
      }
    }
    
    write("</data>",alignment_file,append=TRUE)
  }
}





# sample based on TargetBC
# to do check: can there be several copies of the same targetBC per cell? NO, I filtered them out ;) S

collect_sequences_per_cell_targetBC <- function(cell_sample, dataset, targetBC) {
  #select all sequences with a specific target BC from a list of cells
  
  sequence_collection <- c()
  
  for(i in cell_sample) {
    
    cell_sequence = dataset[(dataset$Cell == i) & (dataset$TargetBC == targetBC),]
    
    if (nrow(cell_sequence) == 0){
      stop(paste("Cell", i, "does not have targetBC", targetBC))
    }
    
    sequence_collection <- rbind(sequence_collection, cell_sequence)
    
  }
  
  return(sequence_collection)
  
}
