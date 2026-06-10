# PTM-mapping
Repository for PTM code 

This R script does three main things: sets up sequence alignment processing, builds a master alignment map, and loads TMT proteomics data.
1. Setup & Alignment Exploration (top section)
Installs and loads two Bioconductor packages — Biostrings (for biological sequence handling) and msa (multiple sequence alignment). It then reads in a single amino acid alignment file and does basic sanity checks: are all sequences the same width? What does the first sequence look like? It also converts the alignment into a matrix and a flat dataframe where every row is one residue at one alignment column for one sequence.
2. The Two Core Functions
process_alignment(file) takes a single .aln-fasta file and:
* Reads the amino acid sequences
* Renames any purely numeric sequence names (a quirk-fix for a specific species, Aa_AgSp1)
* Converts the alignment to a matrix, then to a long-format dataframe with columns for sequence name, alignment column position, and the residue (amino acid or gap) at that position
build_master_from_folder(folder) loops over all .aln-fasta files in a folder, runs process_alignment on each, stacks all the results into one giant master dataframe, assigns a unique integer MasterID to every non-gap residue across all files, and saves the result as a CSV. This is the main workhorse — it turns a folder full of alignment files into a single unified map of every residue.
3. TMT Proteomics Data Loading
Reads tab-delimited .txt files from two folders — TMT_modifications_raw and TMT_quantified_raw — into one combined dataframe, skipping any malformed files. It then splits that combined data back into two separate dataframes: one for peptide modifications data and one for quantification data.
