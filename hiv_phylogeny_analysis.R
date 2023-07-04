#### This script runs a phylogenetic analysis and produces a results directory containing:
# 1: Phylogenetic tree pdf file.
# 2: A csv file with repoting metadata from the target sequences.
# 2: Two csv files containing the nucleotide differences between the sequences and their genetic distance.
# 3: Multiple sequence aligment fasta file.
# 4: Phylogenetic tree raxml files.
# The script requires an existing ~/Data directory containing four files downloaded from GSAID.
# The tree shows the sequences lineages and a second label:
# By default the second label is the patients age but it can be changed manually to show other attributes.

#### Libraries
library(ggplot2)
library(ggtree)
library(treeio) # Load before seqinr to avoid "read.fasta" functions overlapping.
library(seqinr)
library(ape)
library(phangorn) # Get midpoint of tree to re root it.
library(tidytree)

#### CHANGE the analysis_tag, data_path, file_names, target_IDs and second lable according to the specific analysis.
fasta_path <- "/home/josemari/Desktop/Jose/Projects/HIV_resistance/Data/HIV_int_mafft_aligned.fasta"
fasta_file <- read.fasta(fasta_path, as.string = TRUE, forceDNAtolower = FALSE, set.attributes = FALSE)
# Filter for specific IDs
tree_IDs <- names(fasta_file)
# Get the secondary variable label to plot.
#patient_age <- c(); for (i in 1:length(target_IDs)){patient_age<-c(patient_age, metadata_file$Patient.age[i])}
variable1 <- rep("HIV", length(fasta_file))
variable2 <- rep("darkblue", length(fasta_file))

#raxml_command = "raxmlHPC-PTHREADS-AVX -T 12 -f a -x 123 -p 123 -N 100 -m GTRCAT -k -O -s mafft_aligned.fasta -n raxml_tree -w `pwd`"

# Make the data frame for ploting.
tree_df <- data.frame(tree_IDs, variable1, variable2)
colnames(tree_df) <- c("Accession.ID", "variable1", "variable2")
print(tree_df)

# Read tree and re root it to the reference.
tree <- read.tree(file = "/home/josemari/Desktop/Jose/Projects/HIV_resistance/Data/HIV_int_RAxML_bipartitionsBranchLabels.raxml_tree")
rooted_tree <- midpoint(tree, node.labels = "support")
ggtree(rooted_tree) %<+% tree_df


# Plot tree
plot_tree <- ggtree(rooted_tree) %<+% tree_df +
  geom_tiplab(aes(fill = factor(variable2)),
              colour = "black",
              geom = "label",
              size = 3.5,
              label.padding = unit(0.2, "lines"), # amount of padding around the labels
              label.size = 0.2) + # size of label borders
  xlim(0, 0.0019) +
  guides(fill = guide_legend(title = "variable2")) +
  geom_tippoint(aes(color = factor(variable1)),
                size = 3,
                alpha = 1) +
  #scale_fill_hue() +
  guides(color = guide_legend(title = "variable1")) +
  theme_tree2()
plot_tree
# Save the plot as a png.
ggsave("phylogenetic_tree.png", plot_tree, device = "png", height = 297, width = 210, units = "mm", dpi = 400)

# Read the mafft aligned fasta.
aligned_fasta <- read.fasta("mafft_aligned.fasta", as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
# Split all sequences on individual nt for to do the distance calculation.
for (i in 1:length(aligned_fasta)){
  aligned_fasta[[i]] <- unlist(strsplit(aligned_fasta[[i]], ""))
}
# Count nucleotide differences between sequences and raw genetic distance.
nt_diff <- dist.dna(as.DNAbin(aligned_fasta), model = "N")
gen_distance <- dist.dna(as.DNAbin(aligned_fasta), model = "raw")

# Add extra column with rownames because rownames are lost when saving the "dist" objects as tables.
nt_diff <- cbind(as.matrix(nt_diff),  rownames(as.matrix(nt_diff)))
gen_distance <- cbind(as.matrix(gen_distance), rownames(as.matrix(gen_distance)))
# Write csv files.
write.table(as.matrix(nt_diff), "nt_changes.csv", row.names = FALSE)
write.table(as.matrix(gen_distance), "genetic_distance.csv", row.names = FALSE)

# Move all files to the results directory.
system(paste0("mv genetic_distance.csv nt_changes.csv phylogenetic_tree.png mafft_aligned.fasta ref_appended.fasta *.raxml_tree ", results_dir_name))


# Read as raxml for bootstrap values.
raxml_tree <- read.tree(file = "/home/josemari/Desktop/Jose/Projects/HIV_resistance/Data/HIV_int_RAxML_bipartitionsBranchLabels.raxml_tree")

# Get the node of the reference to re root the tree.
noderoot <-  nodeid(raxml_tree, "NC_001802.1_int_3kto5k")
raxml_rooted_tree <- root(raxml_tree, outgroup = noderoot)

ggtree(raxml_rooted_tree) %<+% tree_df +
  geom_label(aes(label=bootstrap, fill=bootstrap)) +
  geom_tiplab() +
  geom_tippoint(aes(color = factor(variable1)),
                size = 3,
                alpha = 1) +
  guides(color = guide_legend(title = variable1)) +
  #geom_tip(aes(color = factor(variable2))) +
  xlim(0, 0.0008) +
  theme_tree2()
