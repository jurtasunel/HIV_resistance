### This script takes in spredsheets of HIV samples and filters out only the samples with specific mutations. 

#### Libraries
library(seqinr)
library(ggplot2)
library(RColorBrewer)

### Functions:
# Write fasta file from seqinr fasta object. Inputs are fasta object and string for file name.
write_fasta_seqinr <- function(fasta_seqs, file_name){
  # Get all sequences on a vector.
  seqs <- lapply(1:length(fasta_seqs), function(x) fasta_seqs[[x]])
  write.fasta(sequences = seqs, names = names(fasta_seqs), file.out = (file_name))
  # Print location.
  print(paste0("Fasta write to: ", getwd()))
}

# Load spreadsheet data.
data_path = "/home/josemari/Desktop/Jose/Projects/HIV_resistance/Results/HIV_drug_resistance_mutation_INT.csv"
data_df <- read.csv(data_path, stringsAsFactors = FALSE)
# Load fasta file.
fasta_path <- "/home/josemari/Desktop/Jose/Projects/HIV_resistance/Data/HIV_int_mafft_aligned.fasta"
fasta_file <- read.fasta(fasta_path, as.string = TRUE, forceDNAtolower = FALSE, set.attributes = FALSE)

# Set sample IDs as rownames
rownames(data_df) <- data_df[,1]
data_df$Sample.ID <- NULL

# Get the total columns, which is the total number of mutations.
total_vars <- ncol(data_df)
# Make an empty vector to store the IDs of the samples with no variants.
novar_samples <- c()
# Loop through the rows of the data frame.
for (i in 1:nrow(data_df)){
  
  # Check if the number of NA in the row equals to the total column.
  if (length(which(is.na(data_df[i,]))) == total_vars){
    # If it does, it means that there are no variants. Get the sample ID and append it to the vector.
    novar_samples <- c(novar_samples, rownames(data_df)[i])
  }
}

# Remove the samples without mutations from the data frame and the fasta.
data_df <- data_df[!(rownames(data_df) %in% novar_samples), ]
fasta_file <- fasta_file[! names(fasta_file) %in% novar_samples]
# Write out the fasta file with only the samples with mutations.
write_fasta_seqinr(fasta_file, "HIV_drug_resistance_mutation_INT.fasta")

# Make a vector to store all the occurrences of each mutations.
freq <- c()
# Loop through the columns and get the total sum.
for (i in 1:ncol(data_df)){
  # Add all the occurrences of that mutation excluding the NA values.
  freq <- c(freq, sum(data_df[,i], na.rm = TRUE))
}

# Make data frame to plot the frequency of the mutations.
piechart_df <- data.frame(colnames(data_df), freq)
colnames(piechart_df) <- c("Mutations", "Frequency")
piechart_df <- piechart_df[order(-piechart_df$Frequency),]
piechart_df

# Make a summary data grouping all minor variants into "others".
summary_data <- piechart_df[c(1:16),]
othersdata <- piechart_df[c(17:nrow(piechart_df)),]
others <- c("Others", sum(othersdata$Frequency))
summary_data <- rbind(summary_data, others)
# Calculate the percentages of the frequencies.
Percentage <- c()
for (i in 1:nrow(summary_data)){
  i.per <- (as.integer(summary_data$Frequency[i]) / sum(as.integer(summary_data$Frequency))) * 100
  Percentage <- c(Percentage, i.per)
}
summary_data <- cbind(summary_data, Percentage)

# Get the desired number of colours.
number_of_colors <- nrow(summary_data)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(number_of_colors)

# Create the pie chart:
mutations_pie <- ggplot(summary_data, aes(x="", y=as.integer(Percentage), fill=Mutations)) + geom_bar(stat="identity", width=1) + # basic bar plot.
  coord_polar("y", start=0) + geom_text(aes(x = 1.1, label = paste0(round(as.numeric(Percentage), 2), "%")),  # x defines the text distance from the center. 
                                        position = position_stack(vjust = 0.5),
                                        check_overlap = TRUE) + # If labels overlap, don't plot them.
  theme_void() +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_manual(values = mycolors, breaks = as.character(summary_data$Mutations)) +  # Rearrange legend on decreasing freq.
  labs(caption = (paste0("*There are ", nrow(othersdata), " variants with less than 1% representation grouped in 'Others' \n \n *Number of samples: ", nrow(data_df),";  Raw freqquency of the most represented mutation: ", summary_data$Frequency[1]))) +
  theme(plot.margin = unit(c(10,5,10,5), "mm"))
mutations_pie

ggsave("taxon_piechart.pdf", mutations_pie, width = 10, height = 10, dpi = 10, limitsize = FALSE)
