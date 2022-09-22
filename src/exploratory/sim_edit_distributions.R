## ---------------------------
##
## Script name: Simulate number of edits Typewriter
##
## Purpose of script: Check what distributions of edits we can expec
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

setwd("/Users/azwaans/")    # Antoine's working directory (PC)

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

install.packages("TreeSim")

## phylo format 
library(TreeSim)

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


#draw edits numbers along each tree edge == cells, along a tree, assuming a rate lambda
draw_edits_along_edges <- function(tree,lambda=4) {
num_edit_edges <- sapply(tree$edge.length, function(x) {rpois(1,x*lambda)})
return(num_edit_edges)
}

#get the total number of edits obtained per cell sampled at "present". 
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


total_branch_edits_present <- get_edits_per_sampled_cell(tree_present,draw_edits_along_edges(tree_present))
density_plot_present <- ggplot(data=total_branch_edits_present,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
#overlay poisson draw
#geom_histogram(data=poisson,aes(x=length,y=..density..),fill="turquoise",alpha=0.2)

plot_grid(ggtree(tree_present),density_plot_present)
ggsave("ditribution_edit_present_branching.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 

total_branch_edits_past <- get_edits_per_sampled_cell(tree_past,draw_edits_along_edges(tree_past))
density_plot_past <- ggplot(data=total_branch_edits_past,aes(x=edits)) + geom_histogram(aes(y=..density..),alpha = 0.2,fill="orange") + xlab("Number of edits") + ylab("Density") +
  geom_density(alpha=.2, fill="#FF6666") 
#overlay poisson draw
#geom_histogram(data=poisson,aes(x=length,y=..density..),fill="turquoise",alpha=0.2)

plot_grid(ggtree(tree_past),density_plot_past)
ggsave("ditribution_edit_past_branching.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm") 


#test goodness of fit to poisson distribution
observed <- c(as.numeric(table(total_branch_edits$edits)),0)

expected <- sapply(as.numeric(names(table(total_branch_edits))), function(x) dpois(x = x,lambda = mean(total_branch_edits$edits), log = FALSE))
complement <- 1 - sum(expected)
expected_p <- c(expected,complement)

chisq.test(observed,p=expected_p)
chisq.test(observed,p=expected_p,simulate.p.value = TRUE)



       