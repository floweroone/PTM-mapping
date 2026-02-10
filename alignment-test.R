install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")

# Load libraries
library(Biostrings)
library(msa)

# Read alignment
aa_aln <- readAAStringSet("Ms_General_Aligned.aln-fasta")

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

# Load libraries
library(Biostrings)
library(msa)

# Automated alignment-to-dataframe function
process_alignment <- function(alignment_file) {
  
  # Load libraries (only loads if not already loaded)
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings not installed. Run: BiocManager::install('Biostrings')")
  }
  library(Biostrings)
  
  # Read alignment
  message("Reading alignment file...")
  aa_aln <- readAAStringSet(alignment_file)
  
  # Check alignment validity
  message("Checking alignment...")
  widths <- width(aa_aln)
  if (any(widths != widths[1])) {
    warning("Not all sequences have the same width!")
  }
  
  message(paste("Alignment loaded:", length(aa_aln), "sequences,", widths[1], "positions"))
  
  # Convert to matrix
  aa_aln_matrix <- as.matrix(aa_aln)
  
  # Create mapping dataframe
  message("Creating mapping dataframe...")
  aa_sequence_map <- data.frame(
    Sequence = rep(rownames(aa_aln_matrix), each = ncol(aa_aln_matrix)),
    AlnCol   = rep(1:ncol(aa_aln_matrix), times = nrow(aa_aln_matrix)),
    GlobalID = 1:(nrow(aa_aln_matrix) * ncol(aa_aln_matrix)),
    Residue  = as.vector(t(aa_aln_matrix)),
    stringsAsFactors = FALSE
  )
  
  message("Done!")
  
  # Return both the alignment object and the dataframe
  return(list(
    alignment = aa_aln,
    matrix = aa_aln_matrix,
    map = aa_sequence_map
  ))
}

# Store all results in a named list
alignments <- list()
alignments$ms_general <- process_alignment("Ms_General_Aligned.aln-fasta")
# Add more files as needed:
# alignments$ms_specific <- process_alignment("Ms_Specific_Aligned.aln-fasta")
# alignments$other_protein <- process_alignment("Other_Protein.aln-fasta")

# Access components for each alignment:
# Original alignment object
alignments$ms_general$alignment

# Matrix form - VIEW IT
View(alignments$ms_general$matrix)

# Mapping dataframe
View(alignments$ms_general$map)

# Or extract to separate variables if preferred:
ms_general_matrix <- alignments$ms_general$matrix
ms_general_map <- alignments$ms_general$map

#USE THE FUNCTION!!!
# Create storage list (change this up)
alignments <- list()

# Process your file (replace with your actual filename)
alignments$ms_general <- process_alignment("Ms_General_Aligned.aln-fasta")