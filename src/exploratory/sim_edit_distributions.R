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

##getting the edit frequencies and edit types for simulations
edit_table_by_5 = read.csv("data/cell_culture/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))

bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]))) %>% arrange(desc(Freq))
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
tree_past <- sim.rateshift.taxa(n=500,numbsim=1,c(1,10),c(0.3,0.3),
                            c(0.6,0.1),c(0,0.3),complete=FALSE)[[1]]

#tree with a constant birth rate
tree_even <- sim.rateshift.taxa(n=500,numbsim=1,c(10,10),c(0,0.3),
                                c(0.6,0.1),c(0,0.3),complete=FALSE)[[1]]

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
#lambda of 0.2 allows to have close to ~ 5 edits along the tree
#draw edits numbers along each tree edge == cells, along a tree, assuming a rate lambda
#draw inserts types according to the insert frequencies

draw_edits_along_edges <- function(tree,lambda=0.2,insert_types,insert_freq) {
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
  #Specifying a target length of 5 here. Replace for any other length.  
  total_branch_edits[i] <- min(5,sum_path)

}
return(data.frame(edits=total_branch_edits))
}

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
    #this is to enforce a tape length of at most 5. For longer tapes, adapt accordingly 
    print(sequence_path)
    #if(nchar(sequence_path) <= 5 ) {
    full_sequences[i] <- sequence_path
    #} 
    #else{
     # full_sequences[i] <- substr(sequence_path, 1, 5)
    #}
    
  }
  return(full_sequences)
}


######## SIMULATING EDIT NUMBERS AND SEQUENCES ########

### Tree with branching towards present ###

#the total number of edits at tips on the tree with high birth towards the present
#Simulating 13 barcodes with 5 positions each (default) 
total_branch_edits_present1 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present2 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present3 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present4 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present5 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present6 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present7 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present8 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present9 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present10 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present11 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present12 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_present13 <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)

#total at the tips is the sum for all barcodes
summed_total_present <- total_branch_edits_present1 + 
  total_branch_edits_present2 + 
  total_branch_edits_present3 + 
  total_branch_edits_present4 + 
  total_branch_edits_present5 + 
  total_branch_edits_present6 + 
  total_branch_edits_present7 + 
  total_branch_edits_present8 +
  total_branch_edits_present9 +
  total_branch_edits_present10 +
  total_branch_edits_present11 +
  total_branch_edits_present12 +
  total_branch_edits_present13

  
  


#the corresponding barcode sequences
edits_per_edge <- draw_edits_along_edges(tree_present,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$edits
edits_per_edge <- lapply(edits_per_edge, function(x) if(identical(x, character(0))) "" else x)
sim_seqs_present <- get_sequence_per_sampled_cell(tree_present,edits_per_edge)

density_plot_present <- ggplot(data=summed_total_present,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits per cell") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_present),density_plot_present)
ggsave("distribution_edit_present_branching_tape_len_5.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 


### Tree with branching towards past ###

#the total number of edits at tips on the tree with high birth towards the past
#Simulating 13 barcodes with 5 positions each (default) 
total_branch_edits_past1 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past2 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past3 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past4 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past5 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past6 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past7 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past8 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past9 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past10 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past11 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past12 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)
total_branch_edits_past13 <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$num_edit)

#total at the tips is the sum for all barcodes
summed_total_past <- total_branch_edits_past1 + 
  total_branch_edits_past2 + 
  total_branch_edits_past3 + 
  total_branch_edits_past4 + 
  total_branch_edits_past5 + 
  total_branch_edits_past6 + 
  total_branch_edits_past7 + 
  total_branch_edits_past8 +
  total_branch_edits_past9 +
  total_branch_edits_past10 +
  total_branch_edits_past11 +
  total_branch_edits_past12 +
  total_branch_edits_past13

#the corresponding barcode sequences
edits_per_edge <- draw_edits_along_edges(tree_past,insert_freq = bulk_insert_count$Freq,insert_types = bulk_insert_count$Var1)$edits
edits_per_edge <- lapply(edits_per_edge, function(x) if(identical(x, character(0))) "" else x)
sim_seqs_past <- get_sequence_per_sampled_cell(tree_past,edits_per_edge)

density_plot_past <- ggplot(data=summed_total_past,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits per cell") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
plot_grid(ggtree(tree_past),density_plot_past)
ggsave("distribution_edit_past_branching_len5.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 



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
s1 <- 0.1
s2 <- 0.2
s3 <- 0.3
s4 <- 0.4 
s5 <- 0.5
matrix_data <- c(-s1,0,0,0,0,0 ,s1,-s2,0,0,0,0, 0,s2,-s3,0,0,0, 0,0,s3,-s4,0,0, 0,0,0,s4,-s5,0, 0,0,0,0,s5,0)
Qmat <- matrix(matrix_data,nrow=6,ncol=6)
expm(Qmat)
#to draw an appropriate number of events, we have to use the current number of 
       