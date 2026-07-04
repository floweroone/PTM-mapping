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

# ---- 1. Parse a single alignment file into a tidy per-residue table ----
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
  
  bracket_seqs <- names(aa_aln)[grepl("\\[", as.character(aa_aln))]
  if (length(bracket_seqs) > 0) {
    warning("Bracket tags found in: ", paste(bracket_seqs, collapse = ", "))
  }
  
  n_seq <- length(aa_aln)
  seq_lengths <- Biostrings::width(aa_aln)
  residues <- unlist(strsplit(as.character(aa_aln), ""))
  
  aa_sequence_map <- data.frame(
    Sequence   = rep(names(aa_aln), times = seq_lengths),
    SeqIndex   = rep(seq_len(n_seq), times = seq_lengths),
    AlnLength  = rep(seq_lengths, times = seq_lengths),
    Residue    = residues,
    stringsAsFactors = FALSE
  )
  
  # AlnCol = ungapped position within each sequence; NA for gap rows
  aa_sequence_map$AlnCol <- ave(
    aa_sequence_map$Residue != "-",
    aa_sequence_map$Sequence,
    FUN = function(x) {
      pos <- cumsum(x)
      pos[!x] <- NA
      pos
    }
  )
  
  list(
    alignment = aa_aln,
    map       = aa_sequence_map
  )
}

# ---- 2. Normalize gene names into a clean join key ----
# THIS MUST BE DEFINED BEFORE build_master_from_folder() IS CALLED.
strip_to_agsp <- function(x) {
  # Normalize: insert missing underscore before Cfrag/Nfrag (e.g. "AgSp1.1Cfrag" -> "AgSp1.1_Cfrag")
  x <- gsub("([0-9])(Cfrag|Nfrag|frag)", "\\1_\\2", x, perl = TRUE)
  
  # ANY name containing a TRINITY_DN#####_c#_g# pattern collapses to just that bare ID,
  # regardless of what prefix precedes it (Nc_AgSp1_mod_, No_AgSp1_mod_, etc.)
  has_trinity_id <- grepl("TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+", x, perl = TRUE)
  x[has_trinity_id] <- sub(".*(TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+).*", "\\1", x[has_trinity_id], perl = TRUE)
  
  simple_tokens <- c("nt","sy","cr","ct","cta","ctb","rr")
  drop_prefixes <- c("glycine","serine","tr")
  keep_prefixes <- c("Cfrag","Nfrag","frag")
  
  is_rm_segment  <- function(seg) grepl("^(rm|cr)[0-9]*[a-z]?$", seg)
  is_bare_number <- function(seg) grepl("^[0-9]+$", seg)
  
  strip_one <- function(s) {
    segs <- strsplit(s, "_")[[1]]
    n <- length(segs)
    if (n < 2) return(s)
    
    for (i in 2:n) {
      seg <- segs[i]
      if (seg %in% simple_tokens) {
        return(paste(segs[1:(i-1)], collapse = "_"))
      }
      if (seg %in% drop_prefixes && i < n && (is_rm_segment(segs[i+1]) || is_bare_number(segs[i+1]))) {
        return(paste(segs[1:(i-1)], collapse = "_"))
      }
      if (seg %in% keep_prefixes && i < n && is_rm_segment(segs[i+1])) {
        return(paste(segs[1:i], collapse = "_"))
      }
    }
    s
  }
  
  result <- vapply(x, strip_one, character(1), USE.NAMES = FALSE)
  result
}

# ---- 3. Build master table across a whole folder of alignment files ----
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
  
  # GeneID and JoinKey — identifies which whole gene each fragment belongs to
  master_map$GeneID <- sub("_(nt|aa)$", "", as.character(master_map$Sequence))
  master_map$JoinKey <- strip_to_agsp(as.character(master_map$GeneID))
  
  # MasterID = ungapped position pooled across ALL fragments of the same gene (JoinKey),
  # counted continuously in file order — NOT reset per individual Sequence/fragment.
  master_map$MasterID <- NA_integer_
  master_map$JoinFileKey <- paste(master_map$JoinKey, master_map$File, sep = "||")
  master_map <- do.call(rbind, lapply(split(master_map, master_map$JoinFileKey, drop = TRUE), function(df) {
    df$MasterID[df$Residue != "-"] <- seq_len(sum(df$Residue != "-"))
    df
  }))
  rownames(master_map) <- NULL
  master_map$JoinFileKey <- NULL
  
  # GeneID and JoinKey
  master_map <- master_map[, c("MasterID","GeneID","JoinKey","Sequence",
                               "SeqIndex","AlnCol","AlnLength","Residue","File")]
  
  return(list(alignments = alignments, master_map = master_map))
}

# ---- 4. Run the alignment pipeline ----
results <- build_master_from_folder(
  folder = "aligned-fasta-files",
  pattern = "\\.fasta$"
)
master_map <- results$master_map
write.csv(master_map, "aligned-fasta-files/MASTER_alignment_map.csv", row.names = FALSE)

master_map_residues <- master_map %>%
  filter(!is.na(MasterID))
View(master_map_residues)

# ---- 5. Read in TMT quantification data ----
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

# ---- 6. Build JoinKey on TMT side and check match quality ----
TMT_quantified$JoinKey <- strip_to_agsp(TMT_quantified$Gene)

message("Matching genes: ", length(intersect(unique(TMT_quantified$JoinKey), unique(master_map$JoinKey))))

setdiff(unique(TMT_quantified$JoinKey), unique(master_map$JoinKey))
setdiff(unique(master_map$JoinKey), unique(TMT_quantified$JoinKey))

# ---- 7. Join TMT data onto the alignment map ----
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

# ---- 8. QC checks ----
unmatched <- TMT_quantified %>%
  anti_join(master_map_residues, by = c("JoinKey" = "JoinKey", "Number" = "MasterID", "Amino Acid" = "Residue"))
nrow(unmatched)
unique(unmatched$Gene)

letter_mismatches <- TMT_quantified %>%
  inner_join(master_map_residues, by = c("JoinKey" = "JoinKey", "Number" = "MasterID")) %>%
  filter(`Amino Acid` != Residue)

letter_mismatches_clean <- letter_mismatches %>%
  rename(MasterID = Number) %>%
  mutate(AA_difference = paste0("Alignment: ", Residue, " | TMT: ", `Amino Acid`)) %>%
  select(Gene, MasterID, Residue, `Amino Acid`, AA_difference, File, Sequence, source_file) %>%
  arrange(Gene, MasterID)

View(letter_mismatches_clean)
write.csv(letter_mismatches_clean, "letter_mismatches_clean.csv", row.names = FALSE)