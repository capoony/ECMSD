# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(rlang)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript process_files.R <OutputDir> <TaxonColumn> [<Prefix>]")
}

Output <- args[1]
Taxon <- args[2]
Prefix <- if (length(args) >= 3 && nchar(args[3]) > 0) paste0(args[3], "_") else ""

# Create output directory if it doesn't exist
if (!dir.exists(file.path(Output, "mapping"))) {
  dir.create(file.path(Output, "mapping"), recursive = TRUE)
}

# Read summary data (use fread for fast I/O on large files)
summary_file <- file.path(Output, "mapping", paste0(Prefix, "Mito_summary.txt"))

if (!file.exists(summary_file)) {
  stop(paste("Error: Mito summary file", summary_file, "does not exist."))
}

data <- as.data.frame(data.table::fread(summary_file, sep = "\t", header = TRUE))

# Keep only primary alignments (MappingQuality > 0; secondaries always have MQ=0
# in minimap2) and deduplicate so each read is counted exactly once.
data <- data %>%
  filter(MappingQuality > 0) %>%
  distinct(SeqID, .keep_all = TRUE)

# ---------------------------------------------------------------------------
# Resolve the requested rank.
#
# Mito_summary.txt holds one column per rank and writes "NA" wherever a
# reference is not resolved that far — most references stop at species, so
# grouping straight on the `subspecies` column dumps the bulk of the reads into
# a single NA category in every plot. Mirror the fallback LinkTaxonomy.py
# applies to the ref_summary label: walk up to the next coarser rank that is
# populated for that read.
# ---------------------------------------------------------------------------
rank_chain <- c(
  "subspecies", "species", "genus", "family", "order",
  "class", "phylum", "kingdom", "domain"
)

# fread turns the literal "NA" into a real NA, but be robust to either form
blank_to_na <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & trimws(x) %in% c("", "NA")] <- NA_character_
  x
}

# tolerate `-x Subspecies` etc.; the summary columns are all lower case
if (!(Taxon %in% names(data)) && tolower(Taxon) %in% names(data)) {
  Taxon <- tolower(Taxon)
}
if (!(Taxon %in% names(data))) {
  stop(paste0(
    "Taxonomic rank '", Taxon, "' is not a column of ", summary_file,
    ". Available ranks: ", paste(intersect(rank_chain, names(data)), collapse = ", ")
  ))
}

resolved <- blank_to_na(data[[Taxon]])
n_missing <- sum(is.na(resolved))

if (Taxon %in% rank_chain) {
  # ranks coarser than the requested one, from finest to coarsest
  coarser <- rank_chain[(match(Taxon, rank_chain) + 1L):length(rank_chain)]
  for (r in intersect(coarser, names(data))) {
    if (!anyNA(resolved)) break
    fill <- is.na(resolved)
    resolved[fill] <- blank_to_na(data[[r]])[fill]
  }
  n_filled <- n_missing - sum(is.na(resolved))
  if (n_filled > 0) {
    message(sprintf(
      "%d of %d reads have no %s rank; labelled with their next coarser rank instead.",
      n_filled, nrow(data), Taxon
    ))
  }
}

resolved[is.na(resolved)] <- "Unknown"
data[[Taxon]] <- resolved

# Summarize read counts by taxon and length
data_sub <- data %>%
  select(!!sym(Taxon), Length) %>%
  group_by(!!sym(Taxon), Length) %>%
  summarise(TotalReads = n(), .groups = "drop") %>%
  arrange(desc(TotalReads))

# Save summarized table
write.table(
  data_sub,
  file = file.path(Output, "mapping", paste0(Prefix, "Mito_summary_", Taxon, ".txt")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

print("Identifying top 10 taxa by total reads...")

# Identify top 10 taxa by total reads
top_taxa_names <- data_sub %>%
  group_by(!!sym(Taxon)) %>%
  summarise(Total = sum(TotalReads), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10) %>%
  pull(!!sym(Taxon))

top_taxa <- data_sub %>%
  filter((!!sym(Taxon)) %in% top_taxa_names) %>%
  mutate(!!Taxon := factor(!!sym(Taxon), levels = top_taxa_names))

if (length(top_taxa_names) == 0 || nrow(top_taxa) == 0) {
  stop("No taxa found — check input data or filtering conditions.")
}

# Plot histogram of read lengths for top 10 taxa
p1 <- ggplot(top_taxa, aes(x = Length, y = TotalReads, color = !!sym(Taxon))) +
  geom_col(width = 1) +
  facet_wrap(as.formula(paste("~", Taxon)), scales = "free_y") +
  labs(title = "Top 10 Taxa by Read Length", x = "Read Length", y = "Total Reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(Output, "mapping", paste0(Prefix, "Mito_summary_", Taxon, "_ReadLengths.png")),
  plot = p1, width = 10, height = 6, dpi = 300
)

print("Calculating proportions per taxon...")

# Calculate proportions per taxon
data_sub2 <- data_sub %>%
  group_by(!!sym(Taxon)) %>%
  summarise(TotalReads = sum(TotalReads), .groups = "drop") %>%
  mutate(Proportion = TotalReads / sum(TotalReads)) %>%
  arrange(desc(TotalReads)) %>%
  mutate(!!Taxon := factor(!!sym(Taxon), levels = !!sym(Taxon)))

# Save proportions table
write.table(
  data_sub2,
  file = file.path(Output, "mapping", paste0(Prefix, "Mito_summary_", Taxon, "_proportions.txt")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Plot proportions
p2 <- ggplot(data_sub2, aes(x = !!sym(Taxon), y = Proportion)) +
  geom_col() +
  labs(title = "Proportion of Reads per Taxon", x = Taxon, y = "Proportion of Reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(Output, "mapping", paste0(Prefix, "Mito_summary_", Taxon, "_Proportions.png")),
  plot = p2, width = 10, height = 6, dpi = 150
)

print("Plots saved. Processing complete.")
