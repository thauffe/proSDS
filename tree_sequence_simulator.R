library(ape)
library(maps)
library(phyclust)
library(phytools)

args <- commandArgs(trailingOnly = TRUE)


seqgen.to.DNAbin <- function(simseqs, len ) {
  bases <- matrix(nrow = length(simseqs)-1, ncol=len, byrow=TRUE)
  rownames(bases) <- as.character(2:length(simseqs))
  for (s in 2:length(simseqs)) {
    seqstats <- unlist(strsplit(simseqs[s],"\\s+"))
    seqstats[1] <- seqstats[1]
    rownames(bases)[s-1] <- seqstats[1]
    bases[s-1, ] <- unlist(strsplit(tolower(seqstats[2]),''))
  }
  as.DNAbin(bases)
}


ThereIsATree <- FALSE
while(!ThereIsATree){
  Tree <- pbtree(b = as.numeric(args[5]), d = as.numeric(args[6]), n = as.numeric(args[7]), extant.only = TRUE)
  if(!is.null(Tree)){
    ThereIsATree <- TRUE
  }
}
RootAge <- max(branching.times(Tree))
MaxGenDif <- 0.5/2
Len <- 658 # Target length DNA
ScaleFactor <- (1/RootAge) * MaxGenDif

TreeScaled <- Tree
TreeScaled$edge.length <- TreeScaled$edge.length * ScaleFactor
switch (args[1],
  'GTR' = Seq <- seqgen(opts = paste0("-m", args[1], " -r", args[3], ",1 -f", args[2], " -l", 
                                      Len), 
                        rooted.tree = TreeScaled),
  'GTR+G' = Seq <- seqgen(opts = paste0("-m", args[1], " -r", args[3], ",1 -f", args[2], " -a", strsplit(args[4], ",")[[1]][2], " -l", 
                                       Len), 
                         rooted.tree = TreeScaled),
  'GTR+I+G' = Seq <- seqgen(opts = paste0("-m", args[1], " -r", args[3], ",1 -f", args[2], "-i", strsplit(args[4], ",")[[1]][1]," -a", strsplit(args[4], ",")[[1]][2], " -l", 
                                         Len), 
                           rooted.tree = TreeScaled),
  'HKY' = Seq <- seqgen(opts = paste0("-m", args[1], " -t", args[3], ",1 -f", args[2], " -l", 
                                      Len), 
                        rooted.tree = TreeScaled),
  'HKY+G' = Seq <- seqgen(opts = paste0("-m", args[1], " -t", args[3], ",1 -f", args[2], " -a", strsplit(args[4], ",")[[1]][2], " -l", 
                                      Len), 
                        rooted.tree = TreeScaled),
  'HKY+I+G' = Seq <- seqgen(opts = paste0("-m", args[1], " -t", args[3], ",1 -f", args[2], "-i", strsplit(args[4], ",")[[1]][1], " -a", strsplit(args[4], ",")[[1]][2], " -l", 
                                        Len), 
                          rooted.tree = TreeScaled)
)

DnaBin <- seqgen.to.DNAbin(Seq, len = Len)

# Write to database
write.dna(DnaBin, file = paste0(args[8],"Dna.fa"), format = "fasta")
write.tree(Tree, args[9]) # Unscaled!

