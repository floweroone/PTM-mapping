#updated cross species code
# =========================
# PACKAGES
# =========================
library(Biostrings)
library(msa)
library(stringr)
library(dplyr)
library(tinytex)

# =========================
# 1. READ ALL FASTA FILES
# =========================
files <- list.files("aligned-fasta-files",
                    pattern = "\\.fasta$",
                    full.names = TRUE)
all_seqs <- do.call(c, lapply(files, function(f) {
  tryCatch(
    readAAStringSet(f),
    error = function(e) {
      message("Skipping: ", basename(f))
      NULL
    }
  )
}))

# =========================
# 2. BUILD ANNOTATION TABLE
# =========================
df <- data.frame(
  name = names(all_seqs),
  stringsAsFactors = FALSE
)
df$species <- str_extract(df$name, "^[^_]+")
df$gene_raw <- str_extract(df$name, "AgSp1\\.?[0-9]*")

df$gene_class <- case_when(
  str_detect(df$gene_raw, "AgSp1\\.1") ~ "AgSp1.1",
  str_detect(df$gene_raw, "AgSp1\\.2") ~ "AgSp1.2",
  str_detect(df$gene_raw, "AgSp1$|AgSp1_") ~ "AgSp1",
  TRUE ~ "other"
)
df$gene_class[str_detect(df$name, "TRINITY")] <- "AgSp1"

df$type <- str_extract(
  ifelse(str_detect(df$name, "AgSp1"), 
         str_extract(df$name, "AgSp1.*"),
         df$name),
  regex("(rr|tr|serine|glycine|nt|sy|cr|ct[ab]?|cp)(?=_|$)", ignore_case = TRUE)
)
df$type <- tolower(df$type)

# Fold cta/ctb sub-variants into the general ct type
df$type[df$type %in% c("cta", "ctb")] <- "ct"

df$type[is.na(df$type)] <- "other"

# =========================
# 3. EXTRACT REPEAT MODULE INFO
# =========================
df$name_norm <- df$name
df$name_norm <- str_replace(df$name_norm, 
                            "_(glycine|serine|tr)_rm_([0-9])", 
                            "_\\1_rm1_\\2")
df$name_norm <- str_replace(df$name_norm,
                            "_(glycine|serine|tr)_rm([^0-9a-z]|$)",
                            "_\\1_rm1\\2")

df$rm_module <- str_extract(df$name_norm, "rm([0-9]+)[a-z]?", group = 1)
df$rm_module <- as.integer(df$rm_module)
df$rm_letter <- str_extract(df$name_norm, "rm[0-9]+([a-z])", group = 1)
df$rm_letter[is.na(df$rm_letter)] <- "a"

df$repeat_num <- str_extract(df$name, "[0-9]+$")
df$repeat_num <- as.integer(df$repeat_num)
df$repeat_num[is.na(df$repeat_num)] <- 1

# =========================
# 4. DEFINE ALIGNMENT GROUPS -- AgSp1_all ONLY
# =========================
df$rm_group <- case_when(
  df$rm_letter == "a" ~ paste0("rm", df$rm_module),
  TRUE ~ paste0("rm", df$rm_module, df$rm_letter)
)
df$rm_group[is.na(df$rm_module)] <- NA

# Only build the "all" (gene-class-agnostic) grouping -- no per-gene_class split
df$group_all <- ifelse(
  is.na(df$rm_module),
  paste("AgSp1_all", df$type, sep = "_"),
  paste("AgSp1_all", df$type, df$rm_group, sep = "_")
)
df$group_all[df$type == "cp"] <- "AgSp1_all_cp"
df$group_all[df$group_all %in% c("AgSp1_all_tr_rm2", "AgSp1_all_tr_rm3", "AgSp1_all_tr_rm4") & 
               df$species == "Varen"] <- "AgSp1_all_tr"

df$group_all_combined <- case_when(
  str_detect(df$type, "tr")      ~ "AgSp1_all_tr_combined",
  str_detect(df$type, "glycine") ~ "AgSp1_all_glycine",
  str_detect(df$type, "serine")  ~ "AgSp1_all_serine",
  TRUE ~ NA_character_
)

split_all          <- split(all_seqs, df$group_all)
split_all_combined <- split(all_seqs, df$group_all_combined)

split_seqs <- c(split_all, split_all_combined)
split_seqs <- split_seqs[!duplicated(names(split_seqs))]

# =========================
# 5. ORDERING SCHEMES
# =========================

# Existing phylogenetic order (abbreviations as used in your file names)
phylo_order <- c("Amarm", "Lcorn", "No", "Ncruc", "Mlaby", "Mhut",
                 "Aarg", "Atri", "Msag", "Mg", "Varen")
# "Nc" and "Ncruc" are the same species (Neoscona crucifera) -- normalize to one abbreviation
df$species[df$species == "Nc"] <- "Ncruc"

# --- Species name -> abbreviation mapping ---
# NOTE: please verify these, especially the ones not in your original phylo_order.
species_abbrev <- c(
  "Neoscona oaxacensis"     = "No",
  "Metepeira labyrinthea"   = "Mlaby",
  "Micrathena gracilis"     = "Mg",
  "Micrathena sagittata"    = "Msag",
  "Leucauge venusta"        = "Lven",     # <- please confirm abbreviation
  "Verrucosa arenata"       = "Varen",
  "Argiope argentata"       = "Aarg",
  "Neoscona crucifera"      = "Ncruc",
  "Araneus pegnia"          = "Apeg",     # <- please confirm abbreviation
  "Argiope aurantia"        = "Aaur",     # <- please confirm abbreviation
  "Larinioides cornutus"    = "Lcorn",
  "Tetragnatha elongata"    = "Telong",
  "Araneus marmoreus"       = "Amarm",
  "Trichonephila clavipes"  = "Tclav",
  "Argiope trifasciata"     = "Atri"
)

# Property tables, in rank order (top = highest rank)
extensibility_species <- c("Neoscona oaxacensis","Metepeira labyrinthea","Micrathena gracilis",
                           "Micrathena sagittata","Leucauge venusta","Verrucosa arenata",
                           "Argiope argentata","Neoscona crucifera","Araneus pegnia",
                           "Argiope aurantia","Larinioides cornutus","Tetragnatha elongata",
                           "Araneus marmoreus","Trichonephila clavipes","Argiope trifasciata")

toughness_species <- c("Verrucosa arenata","Micrathena gracilis","Metepeira labyrinthea",
                       "Micrathena sagittata","Neoscona crucifera","Larinioides cornutus",
                       "Neoscona oaxacensis","Argiope trifasciata","Araneus marmoreus",
                       "Argiope argentata","Argiope aurantia")

em_species <- c("Verrucosa arenata","Micrathena gracilis","Metepeira labyrinthea",
                "Neoscona crucifera","Micrathena sagittata","Neoscona oaxacensis",
                "Larinioides cornutus","Araneus marmoreus","Argiope argentata",
                "Argiope trifasciata","Argiope aurantia")

extensibility_order <- unname(species_abbrev[extensibility_species])
toughness_order     <- unname(species_abbrev[toughness_species])
em_order            <- unname(species_abbrev[em_species])

# Helper: order sequence names by a given species-priority vector
order_by_species_priority <- function(seq_df, priority_vec) {
  order(match(seq_df$species, priority_vec), seq_df$rm_module, seq_df$repeat_num)
}

# Helper: order by sequence similarity via hierarchical clustering on the alignment itself
# Gap-aware similarity ordering: compares only positions where BOTH sequences
# have a real residue (excludes columns where either sequence has a gap).
# Pairs with zero overlapping positions are flagged and reported, but still
# assigned a distance of 1 (maximally dissimilar) so clustering can proceed.
order_by_similarity <- function(aln_set, group_label = NULL) {
  n <- length(aln_set)
  if (n < 3) return(seq_len(n))
  
  seq_names <- names(aln_set)
  seq_chars <- lapply(as.character(aln_set), function(s) strsplit(s, "")[[1]])
  
  dist_mat <- matrix(0, nrow = n, ncol = n)
  no_overlap_pairs <- character(0)
  
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      a <- seq_chars[[i]]
      b <- seq_chars[[j]]
      
      both_present <- (a != "-") & (b != "-")
      n_compared <- sum(both_present)
      
      if (n_compared == 0) {
        d <- 1
        no_overlap_pairs <- c(no_overlap_pairs, 
                              paste0(seq_names[i], " <-> ", seq_names[j]))
      } else {
        n_matches <- sum(a[both_present] == b[both_present])
        d <- 1 - (n_matches / n_compared)
      }
      
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d
    }
  }
  
  # Flag/report any zero-overlap pairs found
  if (length(no_overlap_pairs) > 0) {
    label_txt <- if (!is.null(group_label)) paste0(" [group: ", group_label, "]") else ""
    message("WARNING: ", length(no_overlap_pairs), 
            " sequence pair(s) with ZERO overlapping positions", label_txt, ":")
    for (pr in no_overlap_pairs) message("    ", pr)
  }
  
  d <- as.dist(dist_mat)
  hc <- hclust(d, method = "average")
  hc$order
}

# =========================
# 6. OUTPUT DIRS
# =========================
dir.create("cross_species_alignments", showWarnings = FALSE)
dir.create("msa_prettyprints_ordered", showWarnings = FALSE)

# =========================
# 7. CROSS-SPECIES ALIGNMENT LOOP
# =========================
for (g in names(split_seqs)) {
  
  region_seqs <- split_seqs[[g]]
  
  if (length(region_seqs) < 2) {
    message("Skipping ", g, " — only ", length(region_seqs), " sequence(s)")
    next
  }
  
  message("Aligning group: ", g, " (", length(region_seqs), " sequences)")
  
  aln <- msa(region_seqs, method = "ClustalOmega")
  aln_set <- as(aln, "AAStringSet")
  
  seq_names <- names(aln_set)
  seq_df    <- df[match(seq_names, df$name), ]
  
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", g)
  
  # Save the raw alignment object once (unordered) for pretty-printing later
  saveRDS(aln,
          file = file.path("cross_species_alignments",
                           paste0(safe_name, "_alignment.rds")))
  
  # ---- 5 ordered FASTA outputs ----
  orderings <- list(
    phylo         = order_by_species_priority(seq_df, phylo_order),
    extensibility = order_by_species_priority(seq_df, extensibility_order),
    toughness     = order_by_species_priority(seq_df, toughness_order),
    EM            = order_by_species_priority(seq_df, em_order),
    similarity    = order_by_similarity(aln_set)
  )
  
  for (ord_name in names(orderings)) {
    ord_index <- orderings[[ord_name]]
    aln_set_ordered <- aln_set[ord_index]
    
    writeXStringSet(
      aln_set_ordered,
      filepath = file.path("cross_species_alignments",
                           paste0(safe_name, "_", ord_name, "_alignment.fasta"))
    )
  }
}

# =========================
# 8. PRETTY PRINT (VISUALIZATION ONLY)
# =========================
# =========================
# 8. PRETTY PRINT (VISUALIZATION ONLY)
# =========================
rds_files <- list.files("cross_species_alignments",
                        pattern = "\\.rds$",
                        full.names = TRUE)

# Raised the ceiling -- large alignments are still worth trying to render
MAX_SEQS_FOR_PRETTYPRINT <- 200

skipped_groups <- data.frame(
  file = character(0),
  n_sequences = integer(0),
  stringsAsFactors = FALSE
)

for (f in rds_files) {
  
  aln <- readRDS(f)
  n_seqs <- length(aln@unmasked)
  
  if (n_seqs > MAX_SEQS_FOR_PRETTYPRINT) {
    message("Skipping ", basename(f), " — ", n_seqs, 
            " sequences (over limit of ", MAX_SEQS_FOR_PRETTYPRINT, "), use Jalview")
    skipped_groups <- rbind(skipped_groups, 
                            data.frame(file = basename(f), n_sequences = n_seqs))
    next
  }
  
  aln_set   <- as(aln, "AAStringSet")
  seq_names <- names(aln_set)
  seq_df    <- df[match(seq_names, df$name), ]
  
  base_name <- tools::file_path_sans_ext(basename(f))
  
  orderings <- list(
    phylo         = order_by_species_priority(seq_df, phylo_order),
    extensibility = order_by_species_priority(seq_df, extensibility_order),
    toughness     = order_by_species_priority(seq_df, toughness_order),
    EM            = order_by_species_priority(seq_df, em_order),
    similarity    = order_by_similarity(aln_set, group_label = base_name)
  )
  
  for (ord_name in names(orderings)) {
    aln_reordered <- aln
    aln_reordered@unmasked <- aln@unmasked[orderings[[ord_name]]]
    
    tryCatch({
      msaPrettyPrint(
        aln_reordered,
        output = "pdf",
        file = file.path(
          "msa_prettyprints_ordered",
          paste0(base_name, "_", ord_name, ".pdf")
        ),
        showNames = "left",
        showLogo = "none",
        showConsensus = "none",
        shadingMode = "similar",
        askForOverwrite = FALSE
      )
    }, error = function(e) {
      message("FAILED to render ", base_name, "_", ord_name, ": ", conditionMessage(e))
      skipped_groups <<- rbind(skipped_groups,
                               data.frame(file = paste0(base_name, "_", ord_name),
                                          n_sequences = n_seqs))
    })
  }
}

# Report everything that got skipped or failed
cat("\n===== SKIPPED / FAILED GROUPS =====\n")
if (nrow(skipped_groups) > 0) {
  print(skipped_groups)
} else {
  cat("None — everything rendered successfully.\n")
}

write.csv(skipped_groups, "msa_prettyprints_ordered/skipped_groups_report.csv", row.names = FALSE)

