## ---------------------------
##
## Script name: Simulate number of edits Typewriter
##
## Purpose of script: Check what distributions of edits we can expect assuming poisson distributed edits along tree
##
## Author: Antoine Zwaans
##
## Date Created: 2022-09-13
##
## Copyright (c) Antoine Zwaans, 2022
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("/Users/azwaans/typewriter_analysis/")      # Antoine's working directory (MAC)

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

install.packages("TreeSim")

## phylo format 
library(TreeSim)
library(ggtree)
library(Matrix)
library(tidyr)
library(dplyr)

##getting the edit frequencies and edit types for simulations
edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))

editbulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]))) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$Freq <- bulk_insert_count$Freq/sum(bulk_insert_count$Freq)

##this can be used for sequence simulations

# The fraction 0.6 of the extant species is included into the final tree
# (the tree has n species AFTER sampling, extinct and
# non-sampled lineages are not included):

#tree with high branching towards the present
tree_present <- sim.rateshift.taxa(n=500,numbsim=1,c(10,1),c(0,0.3),
                   c(0.6,0.1),c(0,0.3),complete=FALSE)[[1]]

#tree with higher birth rate in the past
tree_past <- sim.rateshift.taxa(n=500,numbsim=1,c(1,100),c(0,0.3),
                            c(0.6,0.1),c(0,0.3),complete=FALSE)[[1]]
#first 500 tips
#tree with a constant birth rate
tree_even <- sim.rateshift.taxa(n=3257,numbsim=1,c(10,10),c(0.3,0.3),
                                c(0.6,0.6),c(0.3,0.3),complete=FALSE)[[1]]

# find tree root age (coulnd not find function)

tree_height_calc <- function(tree) { 
  start_edge <- 1
  sum_path <- 0
  while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
    sum_path <- sum_path + tree$edge.length[tree$edge[,2] == start_edge]
    start_edge = tree$edge[tree$edge[,2] == start_edge,1]
    
  }
 return(sum_path)
  
}
present_height <- tree_height_calc(tree_present)
past_height <- tree_height_calc(tree_past)


#scale the tree heights in tree_past to  tree_present  
tree_past$edge.length <- tree_past$edge.length * (present_height/past_height)

#lambda of 4 is the original value I used
#lambda of 0.28 allows to have ~ 5 edits along the tree
#draw edits numbers along each tree edge == cells, along a tree, assuming a rate lambda
#draw inserts types according to the insert frequencies

draw_edits_along_edges <- function(tree,lambda=0.28,insert_types,insert_freq) {
  num_edit_edges <- sapply(tree$edge.length, function(x) {rpois(1,x*lambda)})
  edit_types_edges <- sapply(num_edit_edges, function(x) {sample(insert_types, x, replace = TRUE, prob = insert_freq)})
  return_list <- list("num_edit"=num_edit_edges,"edits"=edit_types_edges)
  return(return_list)
}

#Get the total number of edits obtained per cell sampled at "present". 
#For each tip we sum the number of edits drawn onedge edge leading to the tip

get_edits_per_sampled_cell <- function(tree,num_edit_edges) {
total_branch_edits <- rep(0,tree$Nnode)
for(i in 1:length(tree$tip)) {
  start_edge <- i
  sum_path <- 0
  while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
    sum_path <- sum_path + num_edit_edges[tree$edge[,2] == start_edge]
    start_edge = tree$edge[tree$edge[,2] == start_edge,1]
    
  }
  total_branch_edits[i] <- sum_path

}
return(data.frame(edits=total_branch_edits))
}

### draw edits on branches and simultaneously length(tree$tip)
### I SIMULATED BACKWARDS IN TIME, THAT'S WRONG" 

get_edits_per_sampled_cell_v2 <- function(tree) {
  
  num_inserts_per_edge <- rep(NA,length(tree$edge.length))
  inserts_per_edge <- list()
  full_sequences <- rep("",length(tree$tip))
  total_branch_edits <- rep(0,tree$Nnode)
  
  for(i in 1:length(tree$tip)) {
    print("i")
    print(i)
    start_edge <- which(tree$edge[,2] == i)
    print(start_edge)
    sum_path <- 0
    seq_path <- ""
    
    while(length(start_edge) != 0) {
      ##if we haven't assigned a number of events on this edge, draw one
      if(is.na((num_inserts_per_edge[start_edge]))) {
        #we draw the number of edits from the step matrix row corresponding from the start length
        print("start edge")
        print(start_edge)
        print(tree$edge.length[start_edge])
        pmatrix <- expm(Qmat*tree$edge.length[start_edge])
        increment <- sample(matrix_step[sum_path+1,],size=1,prob=pmatrix[sum_path+1,])
        inserts <- paste(sample(bulk_insert_count$Var1, increment,replace = TRUE, prob = bulk_insert_count$Freq), collapse = '')

        if(identical(inserts, character(0))) {
          inserts_per_edge[start_edge] <- ""
        }
        else { 
          inserts_per_edge[start_edge] <- inserts
        }
        num_inserts_per_edge[start_edge] <- increment
        sum_path <- sum_path + increment
        seq_path <- paste0(seq_path,inserts) 
        start_edge <- which(tree$edge[,2] == tree$edge[start_edge,1])
        
      }
      ##else just use the value previously drawn
      else {
        sum_path <- sum_path + num_inserts_per_edge[start_edge]
        seq_path <- paste0(seq_path,inserts_per_edge[start_edge])
        start_edge <- which(tree$edge[,2] == tree$edge[start_edge,1])
      }
    }
    total_branch_edits[i] <- sum_path
    full_sequences[i] <- seq_path
    
  }
  return(list("edits" = total_branch_edits, "sequences" = full_sequences, "edges"=num_inserts_per_edge))
  #data.frame(edits=)
}


#####V2 simulations#####
v2_total_branch_edits_past <- get_edits_per_sampled_cell_v2(tree_past)
lengths_v2 <- data.frame(edits=v2_total_branch_edits_past$edits)
density_plot_past <- ggplot(data=lengths_v2,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_past),density_plot_past)


v2_total_branch_edits_present <- get_edits_per_sampled_cell_v2(tree_present)
lengths_v2 <- data.frame(edits=v2_total_branch_edits_present$edits)

##attempt at simulation of 13 barcodes: 
v2_total_branch_edits_present <- get_edits_per_sampled_cell_v2(tree_present)$edits
for(i in 2:13) {
  v2_total_branch_edits_present <-  v2_total_branch_edits_present + get_edits_per_sampled_cell_v2(tree_present)$edits
  
}
lengths_v2 <- data.frame(edits=v2_total_branch_edits_present$edits)
density_plot_present <- ggplot(data=lengths_v2,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_present),density_plot_present)




##Get full barcode sequences at tips based on the ones drawn on each edge. 
## traverse the tree based on edge numbers
get_sequence_per_sampled_cell <- function(tree,edit_types_edges) {
  full_sequences <- vector("list", length(tree$tip))
  for(i in 1:length(tree$tip)) {
    start_edge <- i
    sequence_path <- ""
    while(length(tree$edge[tree$edge[,2] == start_edge,1]) != 0) {
      
      sequence_path <- paste0(sequence_path,unlist(edit_types_edges[tree$edge[,2] == start_edge]))
      start_edge = tree$edge[tree$edge[,2] == start_edge,1]
      
    }
    full_sequences[i] <- sequence_path
    
  }
  return(full_sequences)
}


######## SIMULATING EDIT NUMBERS AND SEQUENCES ########

### Tree with branching towards present ###

#the total number of edits at tips on the tree with high birth towards the present
total_branch_edits_present <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)

#the corresponding barcode sequences
edits_per_edge <- draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$edits
edits_per_edge <- lapply(edits_per_edge, function(x) if(identical(x, character(0))) "" else x)
sim_seqs_present <- get_sequence_per_sampled_cell(tree_present,edits_per_edge)

density_plot_present <- ggplot(data=total_branch_edits_present,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_present),density_plot_present)
ggsave("ditribution_edit_present_branching.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 


### Tree with branching towards past ###

#the total number of edits at tips on the tree with high birth towards the past
total_branch_edits_past <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)

#the corresponding barcode sequences
edits_per_edge <- draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$edits
edits_per_edge <- lapply(edits_per_edge, function(x) if(identical(x, character(0))) "" else x)
sim_seqs_past <- get_sequence_per_sampled_cell(tree_past,edits_per_edge)

density_plot_past <- ggplot(data=total_branch_edits_past,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_past),density_plot_past)
ggsave("ditribution_edit_past_branching.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 



######## GOODNESS OF FIT ########

#test goodness of fit to poisson distribution
observed <- c(as.numeric(table(total_branch_edits$edits)),0)

expected <- sapply(as.numeric(names(table(total_branch_edits))), function(x) dpois(x = x,lambda = mean(total_branch_edits$edits), log = FALSE))
complement <- 1 - sum(expected)
expected_p <- c(expected,complement)

chisq.test(observed,p=expected_p)
chisq.test(observed,p=expected_p,simulate.p.value = TRUE)



#overlay poisson draw
#geom_histogram(data=poisson,aes(x=length,y=..density..),fill="turquoise",alpha=0.2)

#try implementing the matrix exponential:
s1 <- 0.068
s2 <- 0.068
s3 <- 0.068
s4 <- 0.068
s5 <- 0.068
matrix_data <- c(-s1,0,0,0,0,0 ,s1,-s2,0,0,0,0, 0,s2,-s3,0,0,0, 0,0,s3,-s4,0,0, 0,0,0,s4,-s5,0, 0,0,0,0,s5,0)
Qmat <- matrix(matrix_data,nrow=6,ncol=6)
expm(Qmat)

matrix_step <- c(0,0,0,0,0,0,1,0,0,0,0,0,2,1,0,0,0,0,3,2,1,0,0,0,4,3,2,1,0,0,5,4,3,2,1,0)
matrix_step <- matrix(matrix_step,nrow=6,ncol=6)
#to draw an appropriate number of events, we have to use the current number of 


#simulation parameters 
# total process is 1.2 million
#

#start ID = n+1


traverse <- function(edge_index=0,tree,edge_vector1,sequence_vector1,sequence_tips1,total_tips1) {
  
 edge_vector <- edge_vector1
 sequence_vector <- sequence_vector1
 sequence_tips <- sequence_tips1
 total_tips <- total_tips1

  #at the root of the tree, just initialise the first 2 descending nodes with edits 
  if(edge_index==0) {
    
    #find both descending edges
    start_indices <- which(tree$edge[,1] == tree$edge[1,1])

    #append edge 1 with new edits
    pmatrix <- expm(Qmat*tree$edge.length[start_indices[1]])
    increment <- sample(matrix_step[1,],size=1,replace=TRUE,prob=pmatrix[edge_vector[start_indices[1]]+1,])
    edge_vector[start_indices[1]] <-  increment
    sequence_vector[start_indices[1]] <- paste(sample(bulk_insert_count$Var1, increment,replace = TRUE, prob = bulk_insert_count$Freq), collapse = '')
  
    #append edge 2 with new edits
    pmatrix <- expm(Qmat*tree$edge.length[start_indices[2]])
    increment <- sample(matrix_step[1,],size=1,replace=TRUE,prob=pmatrix[edge_vector[start_indices[2]]+1,])
    edge_vector[start_indices[2]] <- increment
    sequence_vector[start_indices[2]] <- paste(sample(bulk_insert_count$Var1, increment,replace = TRUE, prob = bulk_insert_count$Freq), collapse = '')

    #start the recursion 
    traverse(start_indices[1],tree = tree,edge_vector1= edge_vector,sequence_vector1= sequence_vector,sequence_tips1=sequence_tips,total_tips1 = total_tips)
    traverse(start_indices[2],tree = tree,edge_vector1= edge_vector,sequence_vector1 = sequence_vector,sequence_tips1=sequence_tips,total_tips1 = total_tips)
    
  }
  
  new_indices <- which(tree$edge[,1] == tree$edge[edge_index,2])

  #reach a tip, record the sequence  
  if((length(new_indices ) == 0 ) ) {
    
    sequence_tips[tree$edge[edge_index,2]] <- sequence_vector[edge_index]
    total_tips[tree$edge[edge_index,2]] <- edge_vector[edge_index]
  }
  
  else {
    
    #append child edge 1 with new edits
    pmatrix <- expm(Qmat*tree$edge.length[new_indices[1]])
    increment <- sample(matrix_step[edge_vector[edge_index]+1,],size=1,prob=pmatrix[edge_vector[edge_index]+1,])
    edge_vector[new_indices[1]] <- edge_vector[edge_index] + increment
    sequence_vector[new_indices[1]] <- paste0(sequence_vector[edge_index],paste(sample(bulk_insert_count$Var1, increment,replace = TRUE, prob = bulk_insert_count$Freq), collapse = ''))
    
    #append child edge 2 with new edits
    pmatrix <- expm(Qmat*tree$edge.length[new_indices[2]])
    increment <- sample(matrix_step[edge_vector[edge_index]+1,],size=1,prob=pmatrix[edge_vector[edge_index]+1,])
    edge_vector[new_indices[2]] <- edge_vector[edge_index] + increment
    sequence_vector[new_indices[2]] <- paste0(sequence_vector[edge_index],paste(sample(bulk_insert_count$Var1, increment,replace = TRUE, prob = bulk_insert_count$Freq), collapse = ''))
  
    #recursively go to descending nodes
    traverse(new_indices[1],tree = tree,edge_vector1= edge_vector,sequence_vector1=sequence_vector,sequence_tips1=sequence_tips,total_tips1 = total_tips)
    traverse(new_indices[2],tree =tree,edge_vector1= edge_vector,sequence_vector1=sequence_vector,sequence_tips1=sequence_tips,total_tips1 = total_tips)
    
  }
  
  ## modify "by reference" parameters
  eval.parent(substitute(edge_vector1<-edge_vector))
  eval.parent(substitute(sequence_vector1<-sequence_vector))
  eval.parent(substitute(sequence_tips1<-sequence_tips))
  eval.parent(substitute(total_tips1<-total_tips))
  
}
       


s1 <- 0.18
s2 <- 0.18
s3 <- 0.18
s4 <- 0.18
s5 <- 0.18
matrix_data <- c(-s1,0,0,0,0,0 ,s1,-s2,0,0,0,0, 0,s2,-s3,0,0,0, 0,0,s3,-s4,0,0, 0,0,0,s4,-s5,0, 0,0,0,0,s5,0)
Qmat <- matrix(matrix_data,nrow=6,ncol=6)

present_height <- tree_height_calc(tree_present)
past_height <- tree_height_calc(tree_past)


#scale the tree heights in tree_past to  tree_present  
tree_past$edge.length <- tree_past$edge.length * (present_height/past_height)


###for the tree with past branching

total_past <- rep(0,500)
for(i in 1:13) {
  edge_vector <- rep(0,998)
  sequence_vector <- rep(0,998)
  sequence_tips <- rep("",500)
  total_tips <- rep(0,500)
  
  traverse(tree = tree_past,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)
  total_past <- total_past + total_tips
  
}

v2_final_attempt <- data.frame(edits=total_past)
density_plot_past <- ggplot(data=v2_final_attempt,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") + geom_vline(xintercept = mean(v2_final_attempt$edits)) +
  geom_vline(xintercept =  mean(v2_final_attempt$edits) + 1.96*sd(v2_final_attempt$edits),linetype = "dashed") + geom_vline(xintercept = mean(v2_final_attempt$edits) - 1.96*sd(v2_final_attempt$edits),linetype = "dashed") 
  
  #geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_past),density_plot_past)
ggsave("ditribution_edit_past_branching.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 


###for the tree with present branching
total_present <- rep(0,500)
for(i in 1:13) {
  edge_vector <- rep(0,998)
  sequence_vector <- rep(0,998)
  sequence_tips <- rep("",500)
  total_tips <- rep(0,500)
  
  traverse(tree = tree_present,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)
  total_present <- total_present + total_tips
  
}
v2_final_attempt <- data.frame(edits=total_present)
density_plot_past <- ggplot(data=v2_final_attempt,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") + geom_vline(xintercept = mean(v2_final_attempt$edits)) + geom_vline(xintercept =  mean(v2_final_attempt$edits)) 

#geom_density(alpha=.2, fill="#FF6666") 

plot_grid(ggtree(tree_present),density_plot_past)
ggsave("ditribution_edit_present_branching.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 


########SIMULATION OF A FULL TYPERWRITER LIKE DATASET#####


#simulation parameters
#25 days
#3257 cells?
#ratios of editing frequencies for all positions
# difference : c(9.71,26.33,28.18,69.34)

#rescale the tree to 25 days:
even_height <- tree_height_calc(tree_even)
tree_even$edge.length <- tree_even$edge.length * (25/even_height)

#simulate sequences 13 targetBCs with 3 that have only 4 positions 1 with only 2
total_dataset <- rep(0,3257)

#the matrix for normal full barcodes is the following:
s1 <- 0.13
s2 <- 0.13
s3 <- 0.13
s4 <- 0.13
s5 <- 0.13
matrix_data <- c(-s1,0,0,0,0,0 ,s1,-s2,0,0,0,0, 0,s2,-s3,0,0,0, 0,0,s3,-s4,0,0, 0,0,0,s4,-s5,0, 0,0,0,0,s5,0)
Qmat <- matrix(matrix_data,nrow=6,ncol=6)




#barcode 1 
edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_1 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips
#barcode 2 

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_2 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

#barcode 3

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_3 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

#barcode 4

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_4 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips


#barcode 5 

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_5 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

#barcode 6 
edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_6 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

#barcode 7 

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_7 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

# barcode 8 
edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_8 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

# barcode 9 
edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_9 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

# barcode 10 only 4 positions

#we modify the matrix to model that! s5 = 0 

s1 <- 0.13
s2 <- 0.13
s3 <- 0.13
s4 <- 0.13
s5 <- 0.000
matrix_data <- c(-s1,0,0,0,0,0 ,s1,-s2,0,0,0,0, 0,s2,-s3,0,0,0, 0,0,s3,-s4,0,0, 0,0,0,s4,-s5,0, 0,0,0,0,s5,0)
Qmat <- matrix(matrix_data,nrow=6,ncol=6)

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_10 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips


#barcode 11 only 4 positions 

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_11 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

#barcode 12 only 4 positions

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_12 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips



#barcode 13 only 2 positions

#only 2 positions
s1 <- 0.13
s2 <- 0.13
s3 <- 0.00
s4 <- 0.00
s5 <- 0.000
matrix_data <- c(-s1,0,0,0,0,0 ,s1,-s2,0,0,0,0, 0,s2,-s3,0,0,0, 0,0,s3,-s4,0,0, 0,0,0,s4,-s5,0, 0,0,0,0,s5,0)
Qmat <- matrix(matrix_data,nrow=6,ncol=6)

edge_vector <- rep(0,6512)
sequence_vector <- rep(0,6512)
sequence_tips <- rep("",3257)
total_tips <- rep(0,3257)
traverse(tree = tree_even,edge_vector1=edge_vector,sequence_vector1 = sequence_vector,sequence_tips1 = sequence_tips,total_tips1 = total_tips)

targetBC_13 <- data.frame(edits=total_tips,sequences=sequence_tips)
total_dataset <- total_dataset + total_tips

#barcode total number of edits 
total_dataset <- data.frame(edits=total_dataset)
density_plot_past <- ggplot(data=total_dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") + geom_vline(xintercept = mean(total_dataset$edits)) + geom_vline(xintercept = mean(total_dataset$edits) +1.96*sd(total_dataset$edits),linetype = "dashed") + geom_vline(xintercept = mean(total_dataset$edits) - 1.96*sd(total_dataset$edits),linetype = "dashed")
#geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_even),density_plot_past)

ggsave("13_barcodes_total_edit_numbers.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 

#plotting individual barcode number distributions: 

dataset <- data.frame(edits=targetBC_1$edits)
plot_1 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_2$edits)
plot_2 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_3$edits)
plot_3 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_4$edits)
plot_4 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_5$edits)
plot_5 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_6$edits)
plot_6 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_7$edits)
plot_7 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_8$edits)
plot_8 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5)) 

dataset <- data.frame(edits=targetBC_9$edits)
plot_9 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_10$edits)
plot_10 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_11$edits)
plot_11 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_12$edits)
plot_12 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") + scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))

dataset <- data.frame(edits=targetBC_13$edits)
plot_13 <- ggplot(data=dataset,aes(x=edits)) + geom_histogram(alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Number of cells") +  scale_x_continuous(limits = c(-0.5, 5.5),breaks = c(0,1,2,3,4,5))



plot_grid(plot_1,plot_2,plot_3,plot_4,plot_5,plot_6,plot_7,plot_8,plot_9,plot_10,plot_11,plot_12,plot_13, nrow = 3,ncol=5)
ggsave("13_barcodes_individual_edit_numbers.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 



diff <- c(9.71,26.33,28.18,69.34)
diff <- 100 - diff



