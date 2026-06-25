library(dplyr)
library(Biostrings)

# =========================
# LOAD AND CLEAN PTM DATA
# =========================

ptm <- read.csv("TMT_quantified_raw/TMT_quantified_joined_reduced.csv")
master_map <- read.csv("aligned-fasta-files/MASTER_alignment_map.csv")

# =========================
# CALCULATE LINE_POS
# line_pos = ungapped residue position within each named sequence (FASTA line)
# resets to 1 at the start of each sequence, skipping gaps
# this matches the cross-species alignment coordinate system
# =========================

# Rebuild line_pos_map with Residue included
line_pos_map <- master_map %>%
  filter(Residue != "-") %>%
  group_by(Sequence) %>%
  arrange(MasterID) %>%
  mutate(line_pos = row_number()) %>%
  ungroup() %>%
  select(MasterID, Sequence, Residue, line_pos)

# Now test the join
ptm %>%
  filter(!is.na(Modification)) %>%
  left_join(line_pos_map, by = c("Sequence", "MasterID", "Residue")) %>%
  summarise(
    total = n(),
    has_line_pos = sum(!is.na(line_pos)),
    missing_line_pos = sum(is.na(line_pos))
  )

# Verify all PTM sequences are covered
cat("Total PTM sequences:               ", length(unique(ptm$Sequence)), "\n")
cat("PTM sequences covered by line_pos: ", 
    sum(unique(ptm$Sequence) %in% unique(line_pos_map$Sequence)), "\n")

# =========================
# LOAD AND CLEAN PTM DATA
# Join line_pos onto PTM data using MasterID
# =========================

ptm_clean <- ptm %>%
  filter(!is.na(Modification)) %>%
  left_join(line_pos_map, by = c("Sequence", "MasterID", "Residue")) %>%
  mutate(fraction = case_when(
    Fraction.modified == "Cannot Determine"       ~ 0,
    Fraction.modified == "Only modified detected" ~ 1,
    TRUE ~ suppressWarnings(as.numeric(Fraction.modified))
  ))

# Check any PTM sites that didn't get a line_pos
missing_line_pos <- ptm_clean %>% filter(is.na(line_pos))
if (nrow(missing_line_pos) > 0) {
  message("WARNING: ", nrow(missing_line_pos), " PTM sites have no line_pos — check MasterID join")
  print(missing_line_pos %>% select(Sequence, MasterID, Residue, Modification) %>% distinct())
} else {
  message("All PTM sites successfully joined to line_pos")
}

# =========================
# SUMMARISE PTM SITES
# group by Sequence + line_pos (not AlnCol)
# =========================

# Single modification positions
single_positions <- ptm_clean %>%
  group_by(Sequence, line_pos) %>%
  filter(n_distinct(Modification) == 1) %>%
  summarise(
    Modification = first(Modification),
    fraction     = max(fraction),
    Residue      = first(Residue),
    .groups = "drop"
  )

# Positions with both Phos and Hex
both_positions <- ptm_clean %>%
  group_by(Sequence, line_pos) %>%
  filter(n_distinct(Modification) > 1) %>%
  summarise(
    fraction = max(fraction),
    Residue  = first(Residue),
    .groups = "drop"
  ) %>%
  mutate(Modification = "Both")

# Combine into final PTM table
ptm_final <- bind_rows(single_positions, both_positions)

message("PTM modification counts:")
print(table(ptm_final$Modification))

# =========================
# HELPER FUNCTIONS
# =========================

# Map ungapped residue positions (line_pos) to aligned column positions
map_positions <- function(aligned_seq) {
  chars <- strsplit(as.character(aligned_seq), "")[[1]]
  aln_pos <- which(chars != "-")
  setNames(aln_pos, seq_along(aln_pos))
}

# =========================
# JALVIEW COLOR DEFINITIONS
# Jalview format: <featureType>\t<RRGGBB> (no # prefix)
# =========================

color_lines <- c(
  paste("Phos", "2AB0A0", sep = "\t"),   # teal
  paste("Hex",  "FF9937", sep = "\t"),   # orange
  paste("Both", "FF3737", sep = "\t")    # red
)

# =========================
# PROCESS EACH ALIGNMENT
# =========================

dir.create("jalview_features", showWarnings = FALSE)

aln_files <- list.files("cross_species_alignments",
                        pattern = "\\.fasta$",
                        full.names = TRUE)

# Initialise summary table
summary_rows <- list()

for (aln_file in aln_files) {
  
  aln <- readAAStringSet(aln_file)
  aln_name <- tools::file_path_sans_ext(basename(aln_file))
  seqs_in_aln <- names(aln)
  
  ptm_sub <- ptm_final %>% filter(Sequence %in% seqs_in_aln)
  
  if (nrow(ptm_sub) == 0) {
    message("No PTM data for ", aln_name, " — skipping")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      alignment        = aln_name,
      status           = "No PTM data",
      ptm_sites_found  = 0,
      ptm_sites_mapped = 0,
      feature_file     = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  message("Processing ", aln_name, " (", nrow(ptm_sub), " PTM sites)")
  
  # Build position maps for each sequence in the alignment
  # maps line_pos (ungapped residue number) → aligned column number
  pos_maps <- lapply(setNames(seq_along(aln), names(aln)), function(i) {
    map_positions(aln[[i]])
  })
  
  # Map line_pos to aligned column position
  ptm_mapped <- ptm_sub %>%
    rowwise() %>%
    mutate(
      aln_pos = {
        pmap <- pos_maps[[Sequence]]
        if (!is.null(pmap) && line_pos <= length(pmap)) pmap[[line_pos]] else NA_integer_
      }
    ) %>%
    filter(!is.na(aln_pos))
  
  if (nrow(ptm_mapped) == 0) {
    message("No mappable PTM positions for ", aln_name, " — skipping")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      alignment        = aln_name,
      status           = "PTM data found but no positions mappable",
      ptm_sites_found  = nrow(ptm_sub),
      ptm_sites_mapped = 0,
      feature_file     = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # =========================
  # VALIDATE RESIDUE MATCHES
  # Check that Residue from PTM data matches what is in the alignment
  # at the mapped position — flags any remaining mismatches
  # =========================
  
  ptm_mapped <- ptm_mapped %>%
    rowwise() %>%
    mutate(
      seq_chars     = list(strsplit(as.character(aln[[Sequence]]), "")[[1]]),
      aln_residue   = seq_chars[[aln_pos]],
      residue_match = (aln_residue == Residue)
    ) %>%
    ungroup() %>%
    select(-seq_chars)
  
  mismatches <- ptm_mapped %>% filter(!residue_match)
  if (nrow(mismatches) > 0) {
    message("  WARNING: ", nrow(mismatches), " residue mismatches found in ", aln_name)
    print(mismatches %>% select(Sequence, line_pos, aln_pos, Residue, aln_residue, Modification))
    mismatch_file <- file.path("jalview_features", paste0(aln_name, "_MISMATCHES.csv"))
    write.csv(
      mismatches %>% select(Sequence, line_pos, aln_pos, Residue, aln_residue, Modification, fraction),
      mismatch_file, row.names = FALSE
    )
    message("  Mismatches written to: ", mismatch_file)
  } else {
    message("  All residues validated OK")
  }
  
  # =========================
  # BUILD FEATURE LINES
  # Jalview features file format (tab-delimited):
  #   Col 1: description (tooltip — shows modification and residue)
  #   Col 2: sequence ID
  #   Col 3: sequence index (-1 = match by ID)
  #   Col 4: start position (aligned column)
  #   Col 5: end position (aligned column)
  #   Col 6: featureType (matches color header)
  #   Col 7: score (fraction modified, 0-1)
  # =========================
  
  feature_lines <- apply(ptm_mapped, 1, function(row) {
    paste(
      paste0(row["Modification"], " (", row["Residue"], ")"),
      row["Sequence"],
      "-1",
      trimws(row["aln_pos"]),
      trimws(row["aln_pos"]),
      row["Modification"],
      format(round(as.numeric(row["fraction"]), 2), nsmall = 2),
      sep = "\t"
    )
  })
  
  # Combine: color header, blank line, feature lines
  lines <- c(color_lines, "", feature_lines)
  
  out_file <- file.path("jalview_features", paste0(aln_name, ".features"))
  writeLines(lines, out_file)
  
  message("  Written: ", out_file)
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    alignment        = aln_name,
    status           = "Written",
    ptm_sites_found  = nrow(ptm_sub),
    ptm_sites_mapped = nrow(ptm_mapped),
    feature_file     = out_file,
    stringsAsFactors = FALSE
  )
}

# =========================
# WRITE SUMMARY TABLE
# =========================

summary_table <- bind_rows(summary_rows)

write.csv(summary_table, "jalview_features/alignment_summary.csv", row.names = FALSE)

message("\nSummary:")
message("  Written:                        ", sum(summary_table$status == "Written"))
message("  No PTM data:                    ", sum(summary_table$status == "No PTM data"))
message("  PTM found but not mappable:     ", sum(summary_table$status == "PTM data found but no positions mappable"))
message("  Total alignments processed:     ", nrow(summary_table))
message("\nSummary table written to jalview_features/alignment_summary.csv")
message("\nDone!")