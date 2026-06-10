install.packages("BiocManager")
install.packages("rlang")
BiocManager::install("Biostrings")
BiocManager::install("msa")

# Load libraries
library(Biostrings)
library(msa)
library(dplyr)

# Read alignment
aa_aln <- readAAStringSet("aligned-fasta-files/Aa_FASTA_ALIGNED_AND_SECTIONED.fasta")

# View alignment with colors
aa_aln

# Check alignment width (all should match)
width(aa_aln)

# Are there any sequences with different widths?
any(width(aa_aln) != width(aa_aln)[1])

# Inspect first sequence
as.character(aa_aln[[1]])

# Convert to sequence map
aa_sequence_map <- data.frame(
  Sequence = rep(names(aa_aln), times = width(aa_aln)),
  Position = unlist(lapply(width(aa_aln), seq_len)),
  Residue  = unlist(strsplit(as.character(aa_aln), ""))
)

#Finalized automation
process_alignment <- function(alignment_file) {
  
  aa_aln <- Biostrings::readAAStringSet(alignment_file)
  
  seq_names <- names(aa_aln)
  is_num <- grepl("^\\d+$", seq_names)
  seq_names[is_num] <- paste0("Aa_AgSp1_rm1_", seq_names[is_num])
  names(aa_aln) <- seq_names
  
  n_seq <- length(aa_aln)
  seq_lengths <- width(aa_aln)
  
  aa_sequence_map <- data.frame(
    Sequence   = rep(names(aa_aln), times = seq_lengths),
    SeqIndex   = rep(seq_len(n_seq), times = seq_lengths),
    AlnCol     = unlist(lapply(seq_lengths, seq_len)),  # position within each sequence
    AlnLength  = rep(seq_lengths, times = seq_lengths), # each seq's own length
    Residue    = unlist(strsplit(as.character(aa_aln), "")),
    stringsAsFactors = FALSE
  )
  
  list(
    alignment = aa_aln,
    map       = aa_sequence_map  # no matrix returned since seqs are variable length
  )
}

build_master_from_folder <- function(folder = "aligned-fasta-files",
                                     pattern = "\\.fasta$",
                                     output_name = "MASTER_alignment_map.csv") {
  
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No alignment files found in folder.")
  }
  
  alignments <- list()
  
  master_list <- list()
  
  for (i in seq_along(files)) {
    
    file_path <- files[i]
    file_name <- basename(file_path)
    
    message(sprintf("[%d/%d] Processing %s", i, length(files), file_name))
    
    result <- process_alignment(file_path)
    
    alignments[[file_name]] <- result
    
    df <- result$map
    df$File <- file_name
    
    master_list[[i]] <- df
  }
  
  # Combine all files
  master_map <- do.call(rbind, master_list)
  
  # Create unique master ID for each file — gaps (dashes) get NA
  master_map$MasterID <- NA_integer_
  master_map <- do.call(rbind, lapply(split(master_map, master_map$File), function(df) {
    df$MasterID[df$Residue != "-"] <- seq_len(sum(df$Residue != "-"))
    df
  }))
  
  # Clean column order
  master_map <- master_map[, c(
    "MasterID",
    "File",
    "Sequence",
    "SeqIndex",
    "AlnCol",
    "AlnLength",
    "Residue"
  )]
  
  # Save giant master file
  write.csv(master_map,
            file = file.path(folder, output_name),
            row.names = FALSE)
  
  
  return(list(
    alignments = alignments,
    master_map = master_map
  ))
}

results <- build_master_from_folder(
  folder = "aligned-fasta-files",
  pattern = "\\.fasta$"
)

master_map <- results$master_map
View(master_map)

#reading everything in from the folders nicely
library(dplyr)

TMT_folders <- c("TMT_modifications_raw", "TMT_quantified_raw")

df_all <- lapply(TMT_folders, function(folder) {
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
  
  lapply(files, function(f) {
    df <- read.delim(f, sep = "\t", header = TRUE)
    
    # Skip rogue files (e.g. P_tepidarorium_proteins.txt)
    if (ncol(df) < 7) {
      message("Skipping (unexpected format): ", basename(f))
      return(NULL)
    }
    
    df$source_file   <- basename(f)
    df$source_folder <- basename(folder)
    df
  })
}) |> unlist(recursive = FALSE) |> bind_rows()  # <-- this handles mismatched columns

# Split into two named dataframes
TMT_modifications <- df_all[df_all$source_folder == "TMT_modifications_raw", ]
TMT_quantified    <- df_all[df_all$source_folder == "TMT_quantified_raw", ]

