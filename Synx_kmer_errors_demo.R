#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of Synx kmers
# --------------------------------------------------------------------------

# Set working directory (location of data files and scripts)
setwd("~/Desktop/Synx/")

# Load libraries (may need to be installed)
library('runner')
library('stringr')
library('GenomicAlignments')
library('dplyr')

# List datasets to be analysed
data_sets <- c('ont_demo_data_synx', 'illumina_demo_data')

# Store kmer analysis for each dataset in lists
resList <- vector("list", length(data_sets))
names(resList) <- data_sets

for (ds in data_sets) {
  # Calculate error rates for each base using pileup stats
  # --------------------------------------------------------------------------
  
  # Load sequencing error rates
  ONTDNA_pileup <- read.table(paste0(ds, ".bam.bed.tsv"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  
  # Find reference mapping rates
  get_correct <- function(df){
    bases <- c('A', 'T', 'G', 'C')
    ref <- df['REF_NT']
    correct <- bases[bases == ref]
    sum(as.numeric(df[correct]))
  }
  
  ONTDNA_pileup$Depth <- ONTDNA_pileup$Ins + ONTDNA_pileup$Del + ONTDNA_pileup$Coverage
  ONTDNA_pileup$REF_counts <- apply(ONTDNA_pileup, 1, get_correct)
  ONTDNA_pileup$REF_rate <- ONTDNA_pileup$REF_counts / ONTDNA_pileup$Depth
  
  # Find substitution rates
  get_subs <- function(df){
    bases <- c('A', 'T', 'G', 'C')
    ref <- df['REF_NT']
    subs <- bases[!(bases == ref)]
    sum(as.numeric(df[subs]))
  }
  
  ONTDNA_pileup$SUB_counts <- apply(ONTDNA_pileup, 1, get_subs)
  ONTDNA_pileup$SUB_rate <- ONTDNA_pileup$SUB_counts / ONTDNA_pileup$Depth
  
  # Find INDEL rates
  ONTDNA_pileup$Ins_rate <- ONTDNA_pileup$Ins / ONTDNA_pileup$Depth
  ONTDNA_pileup$Del_rate <- ONTDNA_pileup$Del / ONTDNA_pileup$Depth
  ONTDNA_pileup$INDEL_counts <- ONTDNA_pileup$Ins + ONTDNA_pileup$Del
  ONTDNA_pileup$INDEL_rate <- ONTDNA_pileup$INDEL_counts / ONTDNA_pileup$Depth
  
  # Find error rates
  ONTDNA_pileup$Error_counts <- ONTDNA_pileup$SUB_counts + ONTDNA_pileup$INDEL_counts
  ONTDNA_pileup$Error_rate <- ONTDNA_pileup$Error_counts / ONTDNA_pileup$Depth
  
  # Calculate GC content indepenedent and Homopolymer content of kmers
  # --------------------------------------------------------------------------
  
  # Calculate all unique 6mers in Synx
  synx_seq_F <- ONTDNA_pileup$REF_NT
  
  KD_k <- function(x){paste(x, collapse = "")}
  synx_seq_F_6mers <- runner(x = synx_seq_F, k = 6, f = KD_k)
  synx_seq_total_6mers <- synx_seq_F_6mers[str_length(synx_seq_F_6mers) == 6]
  
  # Make ranges of 6mers in synx sequence
  k_total_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_total_6mers), width = 6), strand = "+", SEQUENCE = synx_seq_total_6mers)
  
  # Calculate GC content
  GC_k <- function(x){(str_count(x, pattern = 'C') + str_count(x, pattern = 'G')) / 6 *100}
  k_total_granges$GC_pct <- GC_k(k_total_granges$SEQUENCE)
  
  # Calculate homopolymers in kmers
  hp_max <- data.frame(k_total_granges$SEQUENCE)
  hp_max$HP2 <- grepl("AA|TT|GG|CC", k_total_granges$SEQUENCE)
  hp_max$HP3 <- grepl("AAA|TTT|GGG|CCC", k_total_granges$SEQUENCE)
  hp_max$HP4 <- grepl("AAAA|TTTT|GGGG|CCCC", k_total_granges$SEQUENCE)
  hp_max$HP5 <- grepl("AAAAA|TTTTT|GGGGG|CCCCC", k_total_granges$SEQUENCE)
  hp_max$HP6 <- grepl("AAAAAA|TTTTTT|GGGGGG|CCCCCC", k_total_granges$SEQUENCE)
  k_total_granges$HP_length <- rowSums(hp_max[,2:6]) + 1
  
  k_total_granges$Depth_mean <- mean(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges))
  
  # Calculate error rates for each kmer
  # --------------------------------------------------------------------------
  
  k_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup$Error_counts, k_total_granges@ranges))
  k_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
  k_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
  k_total_granges$Error_upper <- k_total_granges$Error_mean + k_total_granges$Error_sd
  k_total_granges$Error_lower <- k_total_granges$Error_mean - k_total_granges$Error_sd
  
  # Calculate error rates across kmers
  # --------------------------------------------------------------------------
  
  # Calculate error rate per 6mer (taking the mean for repeat kmers)
  kmer_error_freq <- data.frame(k_total_granges@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    dplyr::summarize(Kmer_count = n(),
                     Error_rate_mean = mean(Error_mean, na.rm = TRUE))
  
  # Add to ds2 lists
  resList[[ds]] <- kmer_error_freq
}  

# Make dataframe of mean error rate per kmer for all samples
for(i in names(resList)){
  tmp <- data.frame(resList[[i]])
  
  if(exists('error_df')) {
    idx <- match(error_df$SEQUENCE, tmp$SEQUENCE)
    error_df[,i] <- tmp$Error_rate_mean [idx]
  } else {
    error_df <- tmp[,c("SEQUENCE"), drop = FALSE]
    error_df[,i] <- tmp$Error_rate_mean
  }
}
write.csv(error_df, "6mers_mean_error_rate_per_sample.csv", row.names = FALSE)
