library(ape)
library(maps)
library(phytools)


args = commandArgs(trailingOnly = TRUE)

tree = read.tree(args[1])

x = fastBM(tree, sig2 = as.numeric(args[3]))

write.csv(x, args[2])