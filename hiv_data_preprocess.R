#### This script reads in 2 fasta files (protease/rev.transcriptase and integrase) and HIV.1 reference (fasta and gff) and does the following:
# Removes 

### Libraries:
library(seqinr)
library(stringr)
library(ggplot2)

### Functions:
# Rename fasta headers to avoid duplicates.
rename_duplicates <- function(fasta_seqs){
  
  # Make empty vector to store new headers.
  new_headers <- c()
  
  # Loop through the fasta and get the headers.
  for (i in 1:length(fasta_seqs)){
    i.header <- names(fasta_seqs)[i]
    # If the header is already on the new_headers vector.
    if (length(grep(i.header, new_headers)) > 0){
      # Construct the new header and append it.
      new_header <- paste0(i.header, "_", length(grep(i.header, new_headers)))
      new_headers <- c(new_headers, new_header)
    } else{ # Append the original header as it is.
      new_headers <- c(new_headers, i.header)
    }
  }
  
  # Rename all headers with the new ones and return the updated fasta.
  for (i in 1:length(fasta_seqs)){names(fasta_seqs)[i] <- new_headers[i]}
  
  return(fasta_seqs)
}
# Remove samples below specific length threshold from seqinr fasta object.
remove_short_seqs <- function(fasta_seqs, length_thr){
  
  # Store the fasta headers of the short sequences.
  short_seqs_ids <- c()
  # Loop through the fasta file.
  for (i in 1:length(fasta_seqs)){
    # If the sequence length is below the threshold, append the header.
    if (nchar(fasta_seqs[[i]]) < length_thr){
      short_seqs_ids <- c(short_seqs_ids, names(fasta_seqs[i]))
    }
  }
  # Remove sequences and return the filtered fasta.
  filtered_fasta <- fasta_seqs[! names(fasta_seqs) %in% short_seqs_ids]
  # Print summary.
  print(paste0("Sequences removed with length below ", length_thr, "bp"))
  print(short_seqs_ids)
  
  return(filtered_fasta)
}
# Get the coverage of each sequence of a fasta object.
get_coverage <- function(fasta_seqs){
  
  # Make vectors to store fasta headers and their coverage.
  ids <- c()
  cov <- c()
  
  # Loop through the fasta.
  for (i in 1:length(fasta_seqs)){
    
    # Get the length of the current sequence.
    i.seq_len <- nchar(fasta_seqs[[i]])
    # Count number of non-canonical nucleotides.
    not_acgt <- str_count(fasta_seqs[[i]], "[^acgt]")
    # Calculate the coverage percentage.
    cov_perc <- 100 - (not_acgt / i.seq_len * 100)
    # Append the header and the coverage.
    ids <- c(ids, names(fasta_seqs[i]))
    cov <- c(cov, cov_perc)
  }
  
  # Make coverage df and return it.
  cov_df <- data.frame(ids, cov)
  colnames(cov_df) <- c("IDs", "Coverage")
  return(cov_df)
}
# Write fasta from seqinr fasta object.
write_fasta_seqinr <- function(fasta_seqs, file_name){
  # Get all sequences on a vector.
  seqs <- lapply(1:length(fasta_seqs), function(x) fasta_seqs[[x]])
  write.fasta(sequences = seqs, names = names(fasta_seqs), file.out = (file_name))
  # Print location.
  print(paste0("Fasta write to: ", getwd()))
}

# Load reference data.
reference_path = "/home/josemari/Desktop/Jose/Reference_sequences/HIV"
reference <- read.fasta(paste0(reference_path, "/NC_001802.1.fasta"), as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
ref_gff <- read.table(paste0(reference_path, "/NC_001802.1.gff3"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
# Load fasta files.
data_path = "/home/josemari/Desktop/Jose/Projects/HIV_resistance/Data"
pro_rt_fasta = read.fasta(paste0(data_path, "/ARIA_seqs/pro_rt.fasta"), as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
int_fasta = read.fasta(paste0(data_path, "/MINR_seqs/int.fasta"), as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)

# Rename duplicates.
pro_rt_nodup <- rename_duplicates(pro_rt_fasta)
int_nodup <- rename_duplicates(int_fasta)

# Remove low quality seqs. Seqs below 75%length of total length removed (pro_rt ~ 1200bp so 900 lenght threshold, int ~ 1kb so 750 thresh).
filtered_pro_rt <- remove_short_seqs(pro_rt_nodup, 900)
filtered_int <- remove_short_seqs(int_nodup, 750)

# Get coverage of sequences and low quality sequences.
pro_rt_cov <- get_coverage(filtered_pro_rt)
pro_rt_lowcov <- pro_rt_cov[which(pro_rt_cov$Coverage < 95),]
int_cov <- get_coverage(filtered_int)
int_lowcov <- int_cov[which(int_cov$Coverage < 95),]

# Plot coverage.
# 1. Protease/rev.transcriptase:
pro_rt_plot <- ggplot(pro_rt_cov, aes(x = IDs, y = Coverage)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_hline(yintercept = 95, color = "red") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(),
        axis.text.y = element_text(size = 9),
        axis.ticks.y = element_line(),
        axis.title.x = element_text(),
        axis.text.x = element_text(angle = 90, size = 2.5)) +
  ggtitle("Quality filter: Coverage > 95%") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  labs(caption = paste0("Sequences excluded from analysis: ", nrow(pro_rt_cov[which(pro_rt_cov$Coverage < 95),])))
pro_rt_plot
# 2. Integrase.
int_plot <- ggplot(int_cov, aes(x = IDs, y = Coverage)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_hline(yintercept = 95, color = "red") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(),
        axis.text.y = element_text(size = 9),
        axis.ticks.y = element_line(),
        axis.title.x = element_text(),
        axis.text.x = element_text(angle = 90, size = 2.5)) +
  ggtitle("Quality filter: Coverage > 95%") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  labs(caption = paste0("Sequences excluded from analysis: ", nrow(int_cov[which(int_cov$Coverage < 95),])))
int_plot

# Remove low coverage samples.
filtered_pro_rt <- filtered_pro_rt[! names(filtered_pro_rt) %in% pro_rt_lowcov$IDs]
filtered_int <- filtered_int[! names(filtered_int) %in% int_lowcov$IDs]

# Get the common barcodes between the pro_rt and the int fastas.
common_barcodes <- intersect(names(filtered_pro_rt), names(filtered_int))

# Write out the filtered fasta files.
write_fasta_seqinr(filtered_pro_rt, "filtered_pro_rt.fasta")
write_fasta_seqinr(filtered_int, "filtered_int.fasta")






