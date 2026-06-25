install.packages("BiocManager")
install.packages("rlang")
BiocManager::install("Biostrings")
BiocManager::install("msa")
install.packages("dplyr")

# Load libraries
library(Biostrings)
library(msa)
library(dplyr)

#Finalized automation
process_alignment <- function(alignment_file) {
  
  aa_aln <- Biostrings::readAAStringSet(alignment_file)
  
  seq_names <- names(aa_aln)
  is_num <- grepl("^\\d+$", seq_names)
  seq_names[is_num] <- paste0("Aa_AgSp1_rm1_", seq_names[is_num])
  names(aa_aln) <- seq_names
  
  # Strip annotation tags like [BR1], [NA2] and replace ? with X
  aa_aln <- Biostrings::AAStringSet(gsub("\\[[^\\]]*\\]", "", as.character(aa_aln)))
  aa_aln <- Biostrings::AAStringSet(gsub("\\?", "X", as.character(aa_aln)))
  names(aa_aln) <- seq_names
  
  # Add this inside process_alignment after the gsub, before building the data frame
  bracket_seqs <- names(aa_aln)[grepl("\\[", as.character(aa_aln))]
  if (length(bracket_seqs) > 0) {
    warning("Bracket tags found in: ", paste(bracket_seqs, collapse = ", "))
  }
  
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
                                     pattern = "\\.fasta$") {
  
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No alignment files found in folder.")
  
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
  
  master_map <- do.call(rbind, master_list)
  
  # Extract AGSP
  master_map$AGSP <- sapply(as.character(master_map$Sequence), function(s) {
    m <- regexpr("AgSp[0-9]+(?:\\.[0-9]+)?(?=_|$)", s, perl = TRUE)
    if (m != -1) regmatches(s, m) else "TRINITY"
  })
  
  # Assign MasterID per AGSP x File — ONE block only
  master_map$MasterID <- NA_integer_
  master_map$AgspFile <- paste(master_map$AGSP, master_map$File, sep = "||")
  master_map <- do.call(rbind, lapply(split(master_map, master_map$AgspFile), function(df) {
    if (nrow(df) == 0) return(df)
    df$MasterID[df$Residue != "-"] <- seq_len(sum(df$Residue != "-"))
    df
  }))
  rownames(master_map) <- NULL
  master_map$AgspFile <- NULL
  
  # GeneID and column order
  master_map$GeneID <- sub("_(nt|aa)$", "", as.character(master_map$Sequence))
  master_map <- master_map[, c("MasterID","GeneID","AGSP","Sequence",
                               "SeqIndex","AlnCol","AlnLength","Residue","File")]
  
  return(list(alignments = alignments, master_map = master_map))
}

results <- build_master_from_folder(
  folder = "aligned-fasta-files",
  pattern = "\\.fasta$"
)

master_map <- results$master_map
write.csv(master_map, "aligned-fasta-files/MASTER_alignment_map.csv", row.names = FALSE)
View(master_map)

master_map_residues <- master_map %>%
  filter(!is.na(MasterID))   # gaps carry no real residue number

View(master_map_residues)

#reading everything in from the folders nicely
TMT_quantified <- list.files("TMT_quantified_raw", pattern = "\\.txt$", full.names = TRUE) |>
  lapply(function(f) {
    df <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE)
    
    if (ncol(df) < 7) {
      message("Skipping (unexpected format): ", basename(f))
      return(NULL)
    }
    
    df$source_file <- basename(f)
    df
  }) |>
  bind_rows()
View(TMT_quantified)

# 2. Join TMT_quantified (from your existing folder-reading code) to the alignment map
# Strip to common join key
strip_to_agsp <- function(x) {
  result <- x
  has_agsp <- grepl("AgSp[0-9]+(?:\\.[0-9]+)?", x, perl = TRUE)
  result[has_agsp] <- sub("(.*AgSp[0-9]+(?:\\.[0-9]+)?).*", "\\1", x[has_agsp], perl = TRUE)
  has_trinity <- grepl("TRINITY", x) & !has_agsp
  result[has_trinity] <- sub(".*(TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+).*", "\\1", x[has_trinity])
  result
}

master_map$JoinKey <- strip_to_agsp(master_map$GeneID)
TMT_quantified$JoinKey <- strip_to_agsp(TMT_quantified$Gene)

# Check overlap first
message("Matching genes: ", length(intersect(unique(TMT_quantified$JoinKey), unique(master_map$JoinKey))))

# Join using JoinKey instead of GeneID/Gene
TMT_quantified_joined <- master_map %>%
  left_join(
    TMT_quantified,
    by = c("JoinKey" = "JoinKey",
           "MasterID" = "Number",
           "Residue"  = "Amino Acid")
  )

View(TMT_quantified_joined)

TMT_quantified_joined_reduced <- TMT_quantified_joined %>%
  filter(!is.na(source_file))

View(TMT_quantified_joined_reduced)

write.csv(TMT_quantified_joined_reduced, 
          "TMT_quantified_raw/TMT_quantified_joined_reduced.csv", 
          row.names = FALSE)

# do this check to see if theres anything wrong
# anything that didn't match at all
unmatched <- TMT_quantified %>%
  anti_join(master_map_residues,
            by = c("Gene" = "GeneID", "Number" = "MasterID", "Amino Acid" = "Residue"))
nrow(unmatched)
unique(unmatched$Gene)

# matches on Gene+Number where the amino acid letter disagrees (signals an offset)
letter_mismatches <- TMT_quantified %>%
  inner_join(master_map_residues, by = c("Gene" = "GeneID", "Number" = "MasterID")) %>%
  filter(`Amino Acid` != Residue)

letter_mismatches_clean <- letter_mismatches %>%
  rename(MasterID = Number) %>%
  mutate(AA_difference = paste0("Alignment: ", Residue, " | TMT: ", `Amino Acid`)) %>%
  select(Gene, MasterID, Residue, `Amino Acid`, AA_difference, File, Sequence, source_file) %>%
  arrange(Gene, MasterID)

View(letter_mismatches_clean)
write.csv(letter_mismatches_clean, "letter_mismatches_clean.csv", row.names = FALSE)
