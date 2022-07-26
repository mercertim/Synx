#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of Synx conu ladder
# --------------------------------------------------------------------------

# Set working directory (location of data files and scripts)
setwd("~/Desktop/Synx/")


# Load libraries (may need to be installed)
library('rtracklayer')
library('dplyr')
library('plyr')
library('runner')
library('stringr')
library('spgs')

# Identify unique 31mers in copy number features
# --------------------------------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("synx.bed", genome = 'SynX')

# Add sequence information
synx_seqs <- read.table("synx_gene_sequences.txt")
colnames(synx_seqs) <- c("Feature", "Sequence")
synx_seqs <- synx_seqs %>%
  distinct()

idx <- match(synx_ranges$name, synx_seqs$Feature)
synx_ranges$SEQUENCE <- synx_seqs$Sequence [idx]

# Make list of 31mers for each cn element
KD_k <- function(x){paste(x, collapse = "")}

cnlist.names <- c('1cn', '2cn', '3cn', '4cn')
cnlist <- vector("list", length(cnlist.names))
names(cnlist) <- cnlist.names

for (cn in c('1cn', '2cn', '3cn', '4cn')) {
  cn_seq <- unlist(strsplit(unique(synx_ranges[synx_ranges$name == cn]$SEQUENCE), split = ""))
  cn_31mers <- runner(x = cn_seq, k = 31, f = KD_k)
  cn_31mers <- cn_31mers[str_length(cn_31mers) == 31]
  cnlist[[cn]] <- cn_31mers
}

# Calculate coverage across copy number features
# --------------------------------------------------------------------------

# List datasets to loop script over
data_sets <- c('ont_demo_data_synx')

# Make dataframe of counts (both raw and normalised by mean coverage)
for (ds in data_sets) {
  
  # Make table of 31mers frequency in each sample based on jellyfish
  all_kmer_counts <- read.table("All_31mers_counts.tsv", header = FALSE)
  colnames(all_kmer_counts) <- c("Sequence", "Count")
  
  # Add reverse complement sequence and jellyfish collapses kmers to the lowest lexographically (ie the canonical kmer is the first alphabetically)
  all_kmer_counts$Sequence_RC <- toupper(spgs::reverseComplement(all_kmer_counts$Sequence))
  
  # Identify 31mers found in the copy number elements
  for (cn in c('1cn', '2cn', '3cn', '4cn')) {
    cn_kmers <- cnlist[[cn]]
    all_kmer_counts[,paste0('CN', gsub("cn", "", cn))] <- all_kmer_counts$Sequence %in% cn_kmers | all_kmer_counts$Sequence_RC %in% cn_kmers
  }
  
  all_kmer_counts <- all_kmer_counts %>%
    mutate(
      Copy_number = case_when(
        CN1 == TRUE ~ "CN1",
        CN2 == TRUE ~ "CN2",
        CN3 == TRUE ~ "CN3",
        CN4 == TRUE ~ "CN4"
      )
    )
  
  # Normalise counts by mean across whole SynX
  all_kmer_counts$Count <- as.numeric(all_kmer_counts$Count)
  all_kmer_counts$Normalised_count <- as.numeric(all_kmer_counts$Count/mean(all_kmer_counts$Count))
  all_kmer_counts$Copies <- as.numeric(gsub('CN', '', all_kmer_counts$Copy_number))
  
  # Add to kmer table of raw counts
  if(exists('kmer_table')) {
    idx <- match(kmer_table$Sequence, all_kmer_counts$Sequence)
    kmer_table[,paste0(ds, "_count")] <- all_kmer_counts$Count [idx]
  } else {
    kmer_table <- all_kmer_counts[, c("Sequence", "Sequence_RC", "Copies")]
    kmer_table[,paste0(ds, "_count")] <- all_kmer_counts$Count
  }
  
  # Add to kmer table of normalised counts
  if(exists('kmer_norm_table')) {
    idx <- match(kmer_norm_table$Sequence, all_kmer_counts$Sequence)
    kmer_norm_table[,paste0(ds, "_normalised_count")] <- all_kmer_counts$Normalised_count [idx]
  } else {
    kmer_norm_table <- all_kmer_counts[, c("Sequence", "Sequence_RC", "Copies")]
    kmer_norm_table[,paste0(ds, "_normalised_count")] <- all_kmer_counts$Normalised_count
  }
}

# Subset to kmers in copy number elements
cn_elements_counts <- kmer_table[!is.na(kmer_table$Copies),]
write.csv(cn_elements_counts, "Copy_number_31mers_raw_counts.csv", row.names = FALSE)

cn_elements_norm_counts <- kmer_norm_table[!is.na(kmer_norm_table$Copies),]
write.csv(cn_elements_norm_counts, "Copy_number_31mers_mean_normalised_counts.csv", row.names = FALSE)

