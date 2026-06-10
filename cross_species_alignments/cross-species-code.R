#cross species alingment attempt

library(Biostrings)
library(msa)

# Read all files
files <- list.files("aligned-fasta-files", pattern = "\\.fasta$", full.names = TRUE)

all_seqs <- do.call(c, lapply(files, function(f) {
  tryCatch(readAAStringSet(f), error = function(e) {
    message("Skipping: ", basename(f)); NULL
  })
}))

# Check all sequence names
head(names(all_seqs), 250)

# Just extract region for GROUPING purposes only
# Keep full name as sequence identifier

get_region <- function(seq_name) {
  region <- sub("_[0-9]+[b]?$", "", seq_name)  # remove trailing number
  region <- sub("^.*?(?:FC511|genomic|AgSp1\\.\\d+(?:_[A-Za-z]+)?)_", "", region)
  region
}

regions <- sapply(names(all_seqs), get_region)

# Check grouping
unique(regions)

dir.create("cross_species_alignments", showWarnings = FALSE)

for (r in unique(regions)) {
  
  # Get all sequences for this region (from all species)
  region_seqs <- all_seqs[regions == r]
  
  if (length(region_seqs) < 2) {
    message("Skipping ", r, " — only 1 sequence")
    next
  }
  
  message("Aligning region: ", r, " (", length(region_seqs), " sequences)")
  
  # Full names are preserved automatically — no renaming needed
  aln <- msa(region_seqs, method = "ClustalOmega")
  
  # Save with full species names intact
  aln_ss <- as(aln, "AAStringSet")
  writeXStringSet(aln_ss,
                  filepath = paste0("cross_species_alignments/", r, "_cross_species.fasta"))
}