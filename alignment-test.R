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


#POTENTIAL AUTOMATION!!

# Install packages (run once)
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("msa")

# Automated alignment-to-dataframe function
library(Biostrings)

# ── Core function (unchanged) ─────────────────────────────────────────────────
process_alignment <- function(alignment_file) {
  
  aa_aln <- readAAStringSet(alignment_file)
  
  widths <- width(aa_aln)
  if (any(widths != widths[1])) {
    warning(paste("Not all sequences have the same width in:", alignment_file))
  }
  
  message(paste("Loaded:", basename(alignment_file), 
                "→", length(aa_aln), "sequences,", widths[1], "positions"))
  
  aa_aln_matrix <- as.matrix(aa_aln)
  
  aa_sequence_map <- data.frame(
    Sequence = rep(rownames(aa_aln_matrix), each = ncol(aa_aln_matrix)),
    AlnCol   = rep(1:ncol(aa_aln_matrix), times = nrow(aa_aln_matrix)),
    GlobalID = 1:(nrow(aa_aln_matrix) * ncol(aa_aln_matrix)),
    Residue  = as.vector(t(aa_aln_matrix)),
    stringsAsFactors = FALSE
  )
  
  return(list(
    alignment = aa_aln,
    matrix    = aa_aln_matrix,
    map       = aa_sequence_map
  ))
}

# ── Batch processing ──────────────────────────────────────────────────────────

# Find all .aln-fasta files in the folder
fasta_files <- list.files(
  path       = "aligned-fasta-files",
  pattern    = "\\.aln-fasta$",       # change pattern if your extension differs
  full.names = TRUE
)

message(paste("Found", length(fasta_files), "alignment files to process.\n"))

# Process every file and store results in a named list
alignments <- list()

for (f in fasta_files) {
  # Use the filename (no path, no extension) as the list key
  key <- tools::file_path_sans_ext(basename(f))
  key <- gsub("\\.aln$", "", key)          # strip extra .aln if present
  
  tryCatch({
    alignments[[key]] <- process_alignment(f)
  }, error = function(e) {
    warning(paste("Skipping", f, "→", e$message))
  })
}

# ── Build master dataframe ────────────────────────────────────────────────────

master_map <- do.call(rbind, lapply(names(alignments), function(key) {
  df <- alignments[[key]]$map
  df$SourceFile <- key          # tag every row with which alignment it came from
  df$GlobalID   <- NULL         # drop per-file GlobalID (will regenerate below)
  df
}))

# Fresh global row ID across all files
master_map$GlobalID <- seq_len(nrow(master_map))

# Reorder columns neatly
master_map <- master_map[, c("GlobalID", "SourceFile", "Sequence", "AlnCol", "Residue")]

message(paste("\nMaster dataframe built:", nrow(master_map), "rows,",
              length(unique(master_map$SourceFile)), "alignment files."))

# ── Save / inspect ────────────────────────────────────────────────────────────

# Save to CSV
write.csv(master_map, "master_alignment_map.csv", row.names = FALSE)
message("Saved → master_alignment_map.csv")

# Quick look
View(master_map)
head(master_map)
table(master_map$SourceFile)   # row counts per file


#CHAT GPT CODE!!!!!
process_alignment2 <- function(alignment_file) {
  
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
    
    result <- process_alignment2(file_path)
    
    alignments[[file_name]] <- result
    
    df <- result$map
    df$File <- file_name
    
    master_list[[i]] <- df
  }
  
  # Combine all files
  master_map <- do.call(rbind, master_list)
  
  # Create unique master ID across EVERYTHING
  master_map$MasterID <- seq_len(nrow(master_map))
  
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