install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")

# Load libraries
library(Biostrings)
library(msa)

# Read alignment
aa_aln <- readAAStringSet("aligned-fasta-files/Ms_General_Aligned.aln-fasta")

# View alignment with colors
aa_aln

# Check alignment width (all should match)
width(aa_aln)

# Are there any sequences with different widths?
any(width(aa_aln) != width(aa_aln)[1])

# Inspect first sequence
as.character(aa_aln[[1]])

# Convert to matrix
aa_aln_matrix <- as.matrix(aa_aln)
View(aa_aln_matrix)

aa_sequence_map <- data.frame(
  Sequence = rep(rownames(aa_aln_matrix), each = ncol(aa_aln_matrix)),
  AlnCol   = rep(1:ncol(aa_aln_matrix), times = nrow(aa_aln_matrix)),
  GlobalID = 1:(nrow(aa_aln_matrix) * ncol(aa_aln_matrix)),
  Residue  = as.vector(t(aa_aln_matrix))
)

#Finalized automation
process_alignment <- function(alignment_file) {
  
  aa_aln <- Biostrings::readAAStringSet(alignment_file)
  
  seq_names <- names(aa_aln)
  
  # If a name is purely numeric (e.g., "4"), rename it using its position
  is_num <- grepl("^\\d+$", seq_names)
  seq_names[is_num] <- paste0("Aa_AgSp1_rm1_", seq_names[is_num])
  
  names(aa_aln) <- seq_names
  
  widths <- width(aa_aln)
  if (any(widths != widths[1])) warning("Not all sequences have same width!")
  
  aa_aln_matrix <- as.matrix(aa_aln)
  
  n_seq <- nrow(aa_aln_matrix)
  aln_len <- ncol(aa_aln_matrix)
  
  aa_sequence_map <- data.frame(
    Sequence   = rep(names(aa_aln), each = aln_len),
    SeqIndex   = rep(seq_len(n_seq), each = aln_len),
    AlnCol     = rep(seq_len(aln_len), times = n_seq),
    AlnLength  = aln_len,
    Residue    = as.vector(t(aa_aln_matrix)),
    stringsAsFactors = FALSE
  )
  
  list(
    alignment = aa_aln,
    matrix    = aa_aln_matrix,
    map       = aa_sequence_map
  )
}

build_master_from_folder <- function(folder = "aligned-fasta-files",
                                     pattern = "\\.aln-fasta$",
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
  
  # Create unique master ID across EVERYTHING — gaps (dashes) get NA
  master_map$MasterID <- NA_integer_
  master_map$MasterID[master_map$Residue != "-"] <- seq_len(sum(master_map$Residue != "-"))
  
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
  
  message("Master file written to: ", file.path(folder, output_name))
  
  return(list(
    alignments = alignments,
    master_map = master_map
  ))
}

results <- build_master_from_folder(
  folder = "aligned-fasta-files",
  pattern = "\\.aln-fasta$"
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

