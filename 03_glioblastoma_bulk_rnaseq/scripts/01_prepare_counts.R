suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

message("Reading GEO metadata...")
gse <- getGEO("GSE147352", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

# Identify the tissue-type text robustly from GEO metadata.
source_text <- if ("source_name_ch1" %in% colnames(pheno)) {
  pheno$source_name_ch1
} else {
  rep("", nrow(pheno))
}

# Some GEO records also encode tissue in characteristics fields.
char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
char_text <- if (length(char_cols) > 0) {
  apply(pheno[, char_cols, drop = FALSE], 1, paste, collapse = "; ")
} else {
  rep("", nrow(pheno))
}

combined_text <- paste(source_text, char_text)

condition <- case_when(
  str_detect(str_to_lower(combined_text), "glioblastoma") ~ "glioblastoma",
  str_detect(str_to_lower(combined_text), "normal brain") ~ "normal",
  TRUE ~ NA_character_
)

metadata <- data.frame(
  gsm = rownames(pheno),
  title = pheno$title,
  condition = condition,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(condition))

metadata$condition <- factor(
  metadata$condition,
  levels = c("normal", "glioblastoma")
)

write.table(
  metadata,
  "data/processed/sample_metadata.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Selected samples:")
print(table(metadata$condition))

count_files <- list.files(
  "data/raw/GSE147352",
  pattern = "^GSM.*htseq_count.*\\.txt\\.gz$",
  full.names = TRUE
)

if (length(count_files) == 0) {
  stop(
    "No HTSeq count files found. Run scripts/00_download_geo_counts.sh first."
  )
}

extract_gsm <- function(path) {
  str_extract(basename(path), "^GSM[0-9]+")
}

file_map <- data.frame(
  gsm = vapply(count_files, extract_gsm, character(1)),
  path = count_files,
  stringsAsFactors = FALSE
)

file_map <- inner_join(metadata, file_map, by = "gsm")

if (nrow(file_map) != nrow(metadata)) {
  warning(
    "Metadata samples: ", nrow(metadata),
    "; matched raw-count files: ", nrow(file_map)
  )
}

read_htseq <- function(path, gsm) {
  x <- read_tsv(
    path,
    col_names = c("gene_id", "count"),
    show_col_types = FALSE
  )
  x <- x %>% filter(!str_starts(gene_id, "__"))
  colnames(x)[2] <- gsm
  x
}

message("Building count matrix...")
count_list <- Map(
  read_htseq,
  file_map$path,
  file_map$gsm
)

merged <- Reduce(
  function(x, y) full_join(x, y, by = "gene_id"),
  count_list
)

merged[is.na(merged)] <- 0

count_matrix <- merged %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

storage.mode(count_matrix) <- "integer"

# Reorder columns to metadata
metadata <- metadata[metadata$gsm %in% colnames(count_matrix), ]
count_matrix <- count_matrix[, metadata$gsm, drop = FALSE]

write.csv(
  count_matrix,
  "data/processed/raw_count_matrix.csv",
  quote = FALSE
)

write.table(
  metadata,
  "data/processed/sample_metadata.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Count matrix: ",
  nrow(count_matrix), " genes x ",
  ncol(count_matrix), " samples"
)
