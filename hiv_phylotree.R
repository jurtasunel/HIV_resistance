### This script makes a phylogenetic tree from a multiple sequence alignment fasta.
### The first sequence should be the reference sequence.

# Libraries:
library(ape)
library(TreeTools)

# Load the aligned fasta and the 
fasta_path <- "/home/josemari/Desktop/Jose/Projects/HIV_resistance/Results/HIV_in_samplesmut.fasta"

# Read in the file and get the reference.
aln <- read.dna(fasta_path, format = "fasta")
reference <- aln[1,]
reference_label <- rownames(reference)

# Make a tree using neighbor-joining method.
tree <- njs(dist.dna(aln))
# Root the tree to the reference.
rooted_tree <- RootTree(tree, reference_label)

# Calculate bootstrap values.
bs_replicates = 100
bs <- boot.phylo(rooted_tree, aln, FUN = function(xx) nj(dist.dna(xx)), B = bs_replicates, rooted = TRUE)
# Adjust to percentage for 1000 bootstrap replicates.
bs <- bs/10
# Loop through the bootstraps.
for (i in 1:length(bs)){
  # Replace NA and bs less than 60% by empty string.
  if (is.na(bs[i]) == TRUE | as.numeric(bs[i]) < 60){
    bs[i] <- ""
  # Add the % symbol. 
  } else{
    bs[i] <- paste0(bs[i], "%")}
  }
# Round them, add % symbol and add them to the tree as node labels.
rooted_tree$node.label <- bs

# Get the tip labels.
labels <- rownames(aln)

# Colour the reference blue and the other sequences red.
tip.colours <- c("darkblue", rep("darkred", (length(labels) - 1)))

# Calculate pairwise distances between tips.
distances <- cophenetic(rooted_tree)
# Make a sequence between the minimun and maximun for plotting evenly spaced axis tips.
dist_axes <- round(seq(min(distances), max(distances), by = max(distances)/5), 3)

# Plot the tree with coloured tip labels and showing bootstraps.
tree_plot <- plot(rooted_tree, tip.color = tip.colours, cex = 0.3, align.tip.label = TRUE, show.node.label = TRUE)
# Show genetic distance on x-axis.
axis(side = 1, at = dist_axes, labels = TRUE, tick = TRUE, xlim = range(distances), cex.axis = 0.8)
# Add margin text.
mtext(paste0("*Number of bootstrap replicates: ", bs_replicates), side = 3, cex = 0.8)

### For INT, tree suggest great distance to 22G52143.
### After observing the alignment file, looks like a wrong alignment. Remove sample and redo tree.
aln <- aln[-c(which(rownames(aln) == "22G52143")),]

### For PRO_RT, tree suggest great distance to 22G37125.
### After observing the alignment file, looks like a wrong alignment. Remove sample and redo tree.
aln <- aln[-c(which(rownames(aln) == "22G37125")),]

