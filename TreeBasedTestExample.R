##########################################################################
# An example to illustrate the proposed test
# The data set was created based on a real data set from the study of the upper
# respiratory tract microbiome between smokers (n=28) and nonsmokers (n=32).
# We used the microbiome data from the left side of throat. In this example,
# all the reads were assigned to the leaf nodes. For demonstration purpose,
# we pruned the tree and only kept the leaf nodes with relative abundance
# (average prop > 0.5%, a total of 41 leaf nodes).
# The data set contains the tree ("phylo" class),  the count matrix of the
# leaf nodes (columns are arranged by the leaf node order), and the group 
# indicator
##########################################################################


# Set to relevant directory
# setwd('path-to-data-code')
source('TreeBasedTest.R')
load('Example.RData')

p <- length(tree$tip.label)
edge <- tree$edge
edge.len <- tree$edge.length
nbr <- nrow(edge)
edge.ct <- matrix(0, nbr, n)
# Internal nodes receive zero counts
edge.ct[match(1:p, edge[, 2]), ]  <- t(comm)

cat('Asymptotic P:', KRtest(edge, edge.ct, edge.len, grp.ind), '\n') 

