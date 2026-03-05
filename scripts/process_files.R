#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(rlang)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_mito_summary.R <MitoSummaryPath> <OutputDir> <TaxonColumn> [<Prefix>]")
}

MitoSummaryPath <- args[1]
Output <- args[2]
Taxon <- args[3]
Prefix <- ifelse(length(args) >= 4, paste0(args[4], "_"), "")

# Create output directory if it doesn't exist
if (!dir.exists(file.path(Output, "mapping"))) {
  dir.create(file.path(Output, "mapping"), recursive = TRUE)
}

#Check if MitoSummaryPath exists
if (!file.exists(MitoSummaryPath)) {
  stop(paste("Error: Mito summary file", MitoSummaryPath, "does not exist."))
}

# Read summary data
#summary_file <- file.path(Output, "mapping", "Mito_summary.txt")
data <- read.table(MitoSummaryPath, header = TRUE, sep = "\t")

print(paste("Processing Mito summary data from:", MitoSummaryPath))

# Summarize read counts by taxon and length
data_sub <- data %>%
  select(!!sym(Taxon), Length) %>%
  group_by(!!sym(Taxon), Length) %>%
  summarise(TotalReads = n(), .groups = "drop") %>%
  arrange(desc(TotalReads))

# Save summarized table
write.table(
  data_sub,
  file = file.path(Output, "mapping", paste0(Prefix,"Mito_summary_", Taxon, ".txt")),
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

# Plot histogram of read lengths for top 10 taxa make sure that the barwidth is always 1
p1 <- ggplot(top_taxa, aes(x = Length, y = TotalReads, color = !!sym(Taxon))) +
  geom_col(width = 1) +
  facet_wrap(as.formula(paste("~", Taxon)), scales = "free_y") +
  labs(title = "Top 10 Taxa by Read Length", x = "Read Length", y = "Total Reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(Output, "mapping", paste0(Prefix,"Mito_summary_", Taxon, "_ReadLengths.png")),
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
  file = file.path(Output, "mapping", paste0(Prefix,"Mito_summary_", Taxon, "_proportions.txt")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Plot proportions
p2 <- ggplot(data_sub2, aes(x = !!sym(Taxon), y = Proportion)) +
  geom_col() +
  labs(title = "Proportion of Reads per Taxon", x = Taxon, y = "Proportion of Reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(Output, "mapping", paste0(Prefix,"Mito_summary_", Taxon, "_Proportions.png")),
  plot = p2, width = 10, height = 6, dpi = 150
)

print("Plots saved. Processing complete.")