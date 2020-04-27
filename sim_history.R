library(ape)
library(maps)
library(phytools)

args = commandArgs(trailingOnly = TRUE)

tree = read.tree(args[1])
matrix = as.matrix(read.csv(args[2], header=FALSE, sep=',', quote = '\"'))
dimnames(matrix) <- NULL

features_tree = sim.history(tree, matrix)

write(file=paste(args[3], "/simmap_tree.tre", sep = ""), write.simmap(features_tree, map.order = 'right-to-left'))
write.csv(getStates(features_tree, type = 'tips'), paste(args[3], '/tip_states.csv', sep = ""))