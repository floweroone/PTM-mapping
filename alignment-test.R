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
