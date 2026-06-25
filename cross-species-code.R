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

# clean gene classification
df$gene_class <- case_when(
  str_detect(df$gene_raw, "AgSp1\\.1") ~ "AgSp1.1",
  str_detect(df$gene_raw, "AgSp1\\.2") ~ "AgSp1.2",
  str_detect(df$gene_raw, "AgSp1$|AgSp1_") ~ "AgSp1",
  TRUE ~ "other"
)

# fix Mg TRINITY sequences
df$gene_class[str_detect(df$name, "TRINITY")] <- "AgSp1"

# strip species and gene prefix, then extract type from remainder
# extract type, searching from AgSp1 onward OR from start for TRINITY sequences
df$type <- str_extract(
  ifelse(str_detect(df$name, "AgSp1"), 
         str_extract(df$name, "AgSp1.*"),
         df$name),
  regex("(rr|tr|serine|glycine|nt|sy|cr|ct|cp)(?=_|$)", ignore_case = TRUE)
)
df$type <- tolower(df$type)

# fix Cfrag → cp
df$type[str_detect(df$name, "Cfrag")] <- "cp"
df$type[df$name == "Lcorn_51676_AgSp1.1_Cfrag_rm2_1"] <- "cp"

df$type[is.na(df$type)] <- "other"

# =========================
# 3. EXTRACT REPEAT MODULE INFO
# =========================
# normalize rm_ and rm_N → rm1 and rm1_N
df$name_norm <- df$name

# rm_N → rm1_N (e.g. serine_rm_3 → serine_rm1_3)
df$name_norm <- str_replace(df$name_norm, 
                            "_(glycine|serine|tr)_rm_([0-9])", 
                            "_\\1_rm1_\\2")

# rm without any number → rm1 (e.g. tr_rm → tr_rm1)
df$name_norm <- str_replace(df$name_norm,
                            "_(glycine|serine|tr)_rm([^0-9a-z]|$)",
                            "_\\1_rm1\\2")

# extract module number and submodule letter
df$rm_module <- str_extract(df$name_norm, "rm([0-9]+)[a-z]?", group = 1)
df$rm_module <- as.integer(df$rm_module)
df$rm_letter <- str_extract(df$name_norm, "rm[0-9]+([a-z])", group = 1)
df$rm_letter[is.na(df$rm_letter)] <- "a"  # no letter = standard = a

# extract trailing repeat number
df$repeat_num <- str_extract(df$name, "[0-9]+$")
df$repeat_num <- as.integer(df$repeat_num)
df$repeat_num[is.na(df$repeat_num)] <- 1  # no number = only one = 1

# =========================
# 4. DEFINE ALIGNMENT GROUPS
# =========================
df$rm_group <- case_when(
  df$rm_letter == "a" ~ paste0("rm", df$rm_module),
  TRUE ~ paste0("rm", df$rm_module, df$rm_letter)
)
df$rm_group[is.na(df$rm_module)] <- NA

df$group <- ifelse(
  is.na(df$rm_module),
  paste(df$gene_class, df$type, sep = "_"),
  paste(df$gene_class, df$type, df$rm_group, sep = "_")
)

df$group_all <- ifelse(
  is.na(df$rm_module),
  paste("AgSp1_all", df$type, sep = "_"),
  paste("AgSp1_all", df$type, df$rm_group, sep = "_")
)

# fix Cfrag → cp group (must come after both group assignments)
df$group[str_detect(df$name, "Cfrag")] <- paste(df$gene_class[str_detect(df$name, "Cfrag")], "cp", sep = "_")
df$group_all[str_detect(df$name, "Cfrag")] <- "AgSp1_all_cp"

# ensure AgSp1_all_cp exists
df$group_all[df$type == "cp"] <- "AgSp1_all_cp"

# move singleton Varen AgSp1.2 tr_rm2/3/4 into general tr group
df$group[df$group %in% c("AgSp1.2_tr_rm2", "AgSp1.2_tr_rm3", "AgSp1.2_tr_rm4")] <- "AgSp1.2_tr"
df$group_all[df$group_all %in% c("AgSp1_all_tr_rm2", "AgSp1_all_tr_rm3", "AgSp1_all_tr_rm4") & 
               df$species == "Varen" & df$gene_class == "AgSp1.2"] <- "AgSp1_all_tr"
# create combined tr/glycine/serine groups (all rm modules together)
df$group_combined <- case_when(
  str_detect(df$type, "tr")      ~ paste(df$gene_class, "tr", sep = "_"),
  str_detect(df$type, "glycine") ~ paste(df$gene_class, "glycine", sep = "_"),
  str_detect(df$type, "serine")  ~ paste(df$gene_class, "serine", sep = "_"),
  TRUE ~ NA_character_
)

df$group_all_combined <- case_when(
  str_detect(df$type, "tr")      ~ "AgSp1_all_tr",
  str_detect(df$type, "glycine") ~ "AgSp1_all_glycine",
  str_detect(df$type, "serine")  ~ "AgSp1_all_serine",
  TRUE ~ NA_character_
)

# rebuild all splits
split_primary      <- split(all_seqs, df$group)
split_all          <- split(all_seqs, df$group_all)
split_combined     <- split(all_seqs, df$group_combined)
split_all_combined <- split(all_seqs, df$group_all_combined)

# merge all groups together
split_seqs <- c(split_primary, split_all, split_combined, split_all_combined)
split_seqs <- split_seqs[!duplicated(names(split_seqs))]

# =========================
# 5. OUTPUT DIRS
# =========================
dir.create("cross_species_alignments", showWarnings = FALSE)
dir.create("msa_prettyprints", showWarnings = FALSE)
dir.create("msa_prettyprints_ordered", showWarnings = FALSE)

# =========================
# 6. CROSS-SPECIES ALIGNMENT LOOP
# =========================
phylo_order <- c("Amarm",  # Araneus marmoreus
                 "Lcorn",  # Lariniodes cornutus
                 "No",     # Neoscona oaxacensis
                 "Ncruc",  # Neoscona crucifera
                 "Mlaby",  # Metepeira labyrinthea
                 "Mhut",   # Mastophora hutchinsoni
                 "Aarg",   # Argiope argentata
                 "Atri",   # Argiope trifasiata
                 "Msag",   # Micrathena sagittata
                 "Mg",     # Micrathena gracilis
                 "Varen")  # Verrucosa arenata

for (g in names(split_seqs)) {
  
  region_seqs <- split_seqs[[g]]
  
  if (length(region_seqs) < 2) {
    message("Skipping ", g, " — only ", length(region_seqs), " sequence(s)")
    next
  }
  
  message("Aligning group: ", g,
          " (", length(region_seqs), " sequences)")
  
  aln <- msa(region_seqs, method = "ClustalOmega")
  
  # get sort keys for this group
  seq_names  <- names(as(aln, "AAStringSet"))
  seq_df     <- df[match(seq_names, df$name), ]
  
  sort_key   <- order(
    seq_df$rm_module,
    seq_df$repeat_num,
    match(seq_df$species, phylo_order)
  )
  
  # reorder fasta
  aln_set         <- as(aln, "AAStringSet")
  aln_set_ordered <- aln_set[sort_key]
  
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", g)
  
  # save msa object (unordered, for pretty print)
  saveRDS(aln,
          file = file.path("cross_species_alignments",
                           paste0(safe_name, "_alignment.rds")))
  
  # save fasta (ordered)
  writeXStringSet(
    aln_set_ordered,
    filepath = file.path("cross_species_alignments",
                         paste0(safe_name, "_alignment.fasta"))
  )
}

# =========================
# 7. PRETTY PRINT (VISUALIZATION ONLY)
# =========================
rds_files <- list.files("cross_species_alignments",
                        pattern = "\\.rds$",
                        full.names = TRUE)

for (f in rds_files) {
  
  aln      <- readRDS(f)
  
  if (length(aln@unmasked) > 100) {
    message("Skipping ", basename(f), " — too many sequences, use Jalview")
    next
  }
  
  aln_set  <- as(aln, "AAStringSet")
  seq_names <- names(aln_set)
  seq_df   <- df[match(seq_names, df$name), ]
  
  sort_key <- order(
    seq_df$rm_module,
    seq_df$repeat_num,
    match(seq_df$species, phylo_order)
  )
  
  # reorder inside msa object
  aln_reordered <- aln
  aln_reordered@unmasked <- aln@unmasked[sort_key]
  
  msaPrettyPrint(
    aln_reordered,
    output = "pdf",
    file = file.path(
      "msa_prettyprints_ordered",
      paste0(tools::file_path_sans_ext(basename(f)), ".pdf")
    ),
    showNames = "left",
    showLogo = "none",
    shadingMode = "similar",
    askForOverwrite = FALSE
  )
}