# Analysing the cell culture data

## Data preprocessing

For speed reasons sub-sample data set to 100, 500 and 1000 leaves and
run analyses on this sub-sampled data.


## Analyses

### 1 | 12 targets, 1 strict clock and 19 insertion probabilities to rule them all

Assumption is that all targets evolve under the same substitution
model

Problem: Tree likelihoods do not converge


### 2 | 12 target, 1 strict clock per target

Assumption is that the targets may not evolve under the same clock
(same speed) and that this makes the likelihoods not converge

Problem: Same problem as in 1 seems to appear here

### 3 | 1 target

Assumption is that maybe there is something wrong with using multiple
targets; this could be some problems in the pre-processing of the data
or the assumptions we make on the targets sharing parameters of
sequence evolution. If the analysis converges on a single target, we
have to dig there.

Result: [Running]

### 4 | Fixed tree

Maybe there is something tricky in the tree space that we are trying
to sample. Hence fix the tree and see if then the tree likelihoods do
converge.

Result: [Running]
