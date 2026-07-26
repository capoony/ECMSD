suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

# ---------------------------------------------------------------------------
# Arguments:
#   1  Output directory  (contains mapping/Mito.paf.gz and
#                         mapping/Mito_summary.ref_summary.txt)
#   2  Path to the PAF file (Mito.paf.gz)
#   3  Path to the ref_summary file (Mito_summary.ref_summary.txt)
#   4  Top-N references to plot (default: 100)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_paf_alignments.R <output_dir> <paf_file> <ref_summary_file> [top_n]")
}

out_dir <- args[1]
paf_file <- args[2]
ref_summary <- args[3]
top_n <- if (length(args) >= 4) as.integer(args[4]) else 25L

plot_dir <- file.path(out_dir, "mapping", "alignment_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load reference ranking; drop rows where taxon_name is NA
# ---------------------------------------------------------------------------
ref_tbl <- read.table(ref_summary,
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, quote = ""
)
ref_tbl <- ref_tbl[!is.na(ref_tbl$taxon_name) & ref_tbl$taxon_name != "NA", ]

# Take the top-N references by read count
top_refs <- head(ref_tbl, top_n)

cat(sprintf("Plotting alignments for %d references...\n", nrow(top_refs)))

# ---------------------------------------------------------------------------
# Load full PAF once using data.table::fread for fast I/O, reading only the
# 12 mandatory PAF columns (skip all optional tag fields).
# Pre-filter to the top-N references via an awk pipe to avoid loading the
# entire (potentially multi-GB uncompressed) file into R memory.
# ---------------------------------------------------------------------------
cat("Reading PAF file (fast path via data.table + awk pre-filter)...\n")

paf_cols <- c(
  "qname", "qlen", "qstart", "qend", "strand",
  "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"
)

# Write the set of needed reference names to a temp file for awk lookup
refs_tmp <- tempfile(fileext = ".txt")
writeLines(top_refs$ref_name, refs_tmp)

# Build awk command: load ref names into a hash, then print matching PAF lines
# (column 6 = tname), selecting only the first 12 fields
awk_cmd <- sprintf(
  "gzip -dc %s | awk -F'\\t' 'NR==FNR{refs[$1]=1; next} ($6 in refs){for(i=1;i<=12;i++) printf \"%%s%%s\", $i, (i<12?\"\\t\":\"\\n\")}' %s -",
  shQuote(paf_file), shQuote(refs_tmp)
)

paf_dt <- fread(
  cmd        = awk_cmd,
  sep        = "\t",
  header     = FALSE,
  col.names  = paf_cols,
  data.table = TRUE
)
unlink(refs_tmp)

paf_by_tname <- split(paf_dt, paf_dt$tname)

cat(sprintf(
  "PAF subsetted to %d alignments across %d references.\n",
  nrow(paf_dt), length(paf_by_tname)
))

# ---------------------------------------------------------------------------
# Helper: "nice" breaks for a log1p-scaled depth axis. Default continuous
# breaks are computed on the untransformed values (e.g. evenly-spaced
# 0, 1000, 2000, ...), which then bunch up once log1p-compressed. Using
# powers of ten instead keeps breaks evenly spaced in log space.
# ---------------------------------------------------------------------------
log_depth_breaks <- function(max_val) {
  if (!is.finite(max_val) || max_val <= 0) {
    return(c(0, 1))
  }
  c(0, 10^(0:ceiling(log10(max_val))))
}

# ---------------------------------------------------------------------------
# Helper: compute per-position sequencing depth from a pafr data frame
# Returns a data.frame(pos, depth)
# ---------------------------------------------------------------------------
compute_depth <- function(paf_sub) {
  if (nrow(paf_sub) == 0) {
    return(data.frame(pos = integer(0), depth = integer(0)))
  }
  ref_len <- max(paf_sub$tlen)
  # PAF is 0-based start, 1-based end; convert to 1-based closed intervals
  starts <- pmax(1L, as.integer(paf_sub$tstart) + 1L)
  ends <- pmin(ref_len, as.integer(paf_sub$tend))
  valid <- starts <= ends
  starts <- starts[valid]
  ends <- ends[valid]
  # Vectorised difference array via tabulate(), then cumsum() — no R-level loop
  inc <- tabulate(starts, nbins = ref_len + 1L)
  dec <- tabulate(ends + 1L, nbins = ref_len + 1L)
  depth <- cumsum(inc - dec)[seq_len(ref_len)]
  data.frame(pos = seq_len(ref_len), depth = depth)
}

# ---------------------------------------------------------------------------
# Plot one PDF per reference
# ---------------------------------------------------------------------------
for (i in seq_len(nrow(top_refs))) {
  ref_name <- top_refs$ref_name[i]
  taxid <- top_refs$taxid[i]
  taxon_name <- top_refs$taxon_name[i]
  species_name <- top_refs$species_name[i]
  read_count <- top_refs$read_count[i]

  # O(1) lookup from pre-split list; also drop any rows with NA tname
  paf_sub <- paf_by_tname[[ref_name]]
  if (!is.null(paf_sub)) paf_sub <- paf_sub[!is.na(paf_sub$tname), ]

  if (is.null(paf_sub) || nrow(paf_sub) == 0) {
    cat(sprintf(
      "  [%d/%d] No alignments found for %s – skipping\n",
      i, nrow(top_refs), ref_name
    ))
    next
  }

  # Safe filename: taxid + taxon_name + species epithet only (to avoid repeating genus)
  safe_name <- gsub("[^A-Za-z0-9_]", "_", taxon_name)
  # species_name uses underscores (e.g. "Drosophila_melanogaster"); strip first word
  epithet <- sub("^[^_ ]+[_ ]", "", species_name)
  if (nchar(trimws(epithet)) == 0) epithet <- species_name # fallback if no separator
  safe_epithet <- gsub("[^A-Za-z0-9_]", "_", epithet)
  png_path <- file.path(
    plot_dir,
    sprintf("%04d_taxid%s_%s_%s.png", i, taxid, safe_name, safe_epithet)
  )

  tryCatch(
    {
      ref_len_val <- max(paf_sub$tlen)

      # Per-position depth (also used below to compute breadth-of-coverage %)
      depth_df <- compute_depth(paf_sub)
      coverage_pct <- 100 * sum(depth_df$depth > 0) / ref_len_val

      # --- Panel 1: alignment coverage (breadth) — rectangle per alignment ---
      p_cov <- ggplot(paf_sub) +
        geom_rect(aes(xmin = tstart, xmax = tend, ymin = 0, ymax = 1),
          fill = "black", alpha = 0.3, colour = NA
        ) +
        scale_x_continuous(
          limits = c(0, ref_len_val),
          labels = function(x) paste0(round(x / 1e3, 1), " kb")
        ) +
        labs(
          title = sprintf("%s  [taxID: %s]", gsub("_", " ", species_name), taxid),
          subtitle = sprintf(
            "Rank %d  |  %d mapped reads  |  %.1f%% coverage",
            i, read_count, coverage_pct
          ),
          x = "Reference position (bp)", y = "Coverage"
        ) +
        theme_bw(base_size = 20) +
        theme(
          plot.title = element_text(face = "italic"),
          plot.subtitle = element_text(size = 16, colour = "grey40"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )

      # --- Panel 2: sequencing depth (log scale so low-depth regions
      # remain visible next to high peaks; log1p handles depth == 0) ---
      p_depth <- ggplot(depth_df, aes(x = pos, y = depth)) +
        geom_area(fill = "black", alpha = 0.7) +
        scale_x_continuous(
          limits = c(0, ref_len_val),
          labels = function(x) paste0(round(x / 1e3, 1), " kb")
        ) +
        scale_y_continuous(
          trans = "log1p",
          breaks = log_depth_breaks(max(depth_df$depth, na.rm = TRUE)),
          labels = scales::label_comma()
        ) +
        labs(x = "Reference position (bp)", y = "Depth (×, log scale)") +
        theme_bw(base_size = 20)

      png(png_path, width = 1800, height = 1050, res = 150)
      g1 <- ggplotGrob(p_cov)
      g2 <- ggplotGrob(p_depth)
      # Equalise column widths so plot boxes align on the x-axis
      max_widths <- grid::unit.pmax(g1$widths, g2$widths)
      g1$widths <- max_widths
      g2$widths <- max_widths
      grid.arrange(g1, g2, ncol = 1, heights = c(1, 1.2))
      dev.off()

      cat(sprintf("  [%d/%d] %s -> %s\n", i, nrow(top_refs), taxon_name, png_path))
    },
    error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat(sprintf("  [%d/%d] ERROR for %s: %s\n", i, nrow(top_refs), taxon_name, e$message))
    }
  )
}

cat(sprintf("Done. PDFs written to: %s\n", plot_dir))
