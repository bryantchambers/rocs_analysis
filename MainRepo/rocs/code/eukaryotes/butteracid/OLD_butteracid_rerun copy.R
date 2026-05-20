#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(scales)
  library(data.table)
  library(stringr)
  library(rentrez)
  library(taxonomizr)
  library(tools)
  library(future)
  library(furrr)
  library(tibble)
  library(DescTools)
})

setwd("/projects/caeg/people/ngm902/apps/repos/rocs")
source("code/eukaryotes/butteracid/postprocessing_functions.R")


#####################################################
# Define variables from config
#####################################################
TAX_PATH_NCBI <- "/datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/"
EUK_DIR       <- "taxonomic-profiling/eukaryotes/unicorn_rerun"

# Build SAMPLE_LIST from files in EUK_DIR matching "{SAMPLE}_collapsed.sort.filtered.agg.stat.gz"
files <- list.files(path = EUK_DIR, pattern = "_collapsed\\.sort\\.filtered\\.agg\\.stat\\.gz$", full.names = FALSE)
if (length(files) == 0) stop(sprintf("No files found in %s matching pattern '*_collapsed.sort.filtered.agg.stat.gz'", EUK_DIR))
samples <- sub("_collapsed\\.sort\\.filtered\\.agg\\.stat\\.gz$", "", files)
SAMPLE_LIST <- file.path(EUK_DIR, "sample_list.tsv")
write.table(data.frame(library_id = samples), SAMPLE_LIST, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
log_msg(paste("Created sample list:", SAMPLE_LIST, "with", length(samples), "samples"))

#####################################################
# Logging helper
#####################################################

log_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}


#####################################################
# Cache helpers
#####################################################
cache_read <- function(path, read_fun, ...) {
  if (file.exists(path) && file.info(path)$size > 0) {
    log_msg(paste("Cache hit:", path))
    return(read_fun(path, ...))
  }
  NULL
}

cache_write <- function(df, path, write_fun, ...) {
  write_fun(df, path, ...)
  log_msg(paste("Wrote:", path))
  df
}

#####################################################
# Helpers
#####################################################
extract_library_id <- function(file) {
  b <- basename(file)
  # Strip everything from .sort. or .filtered. onward
  b <- sub("\\.(sort|filtered)\\..*$", "", b)
  b
}

read_with_id <- function(file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    message(paste("Skipping missing or empty file:", file))
    return(NULL)
  }

  lib_id <- extract_library_id(file)

  fread(file, colClasses = "character", data.table = FALSE) |>
    as_tibble() |>
    mutate(library_id = lib_id)
}

read_or_process <- function(output_file, path = NULL, pattern = NULL) {
  if (file.exists(output_file) && file.info(output_file)$size > 0) {
    log_msg(paste("Reading existing file:", output_file))
    return(fread(output_file, sep = "\t", data.table = FALSE))
  } else {
    if (is.null(path) || is.null(pattern)) stop("Need path and pattern to process raw files")
    files <- list.files(path = path, pattern = pattern, full.names = TRUE, recursive = TRUE)
    if (length(files) == 0) stop(paste("No files found for pattern:", pattern))
    log_msg(paste("Found", length(files), "files for pattern:", pattern))
    df <- map_dfr(files, read_with_id)
    fwrite(df, output_file, sep = "\t", quote = FALSE)
    log_msg(paste("Processed and saved:", output_file))
    df
  }
}

clean_taxstats_columns <- function(df, rank_name) {
  rank_col <- switch(
    rank_name,
    species = "species",
    genus   = "genus",
    family  = "family",
    stop("rank_name must be one of: species, genus, family")
  )

  df |>
    (\(x) setNames(x, sub("^#taxid$", "taxid", names(x))))() |>
    rename(
      !!rank_col := name,
      n_references = num_accessions,
      reference_length = total_length,
      total_reads_from_unicorn = num_reads
    ) |>
    mutate(
      library_id = as.character(library_id),
      taxid = as.character(taxid),
      !!rank_col := as.character(.data[[rank_col]])
    )
}

add_lineage_paths <- function(df) {
  df |>
    mutate(
      lineage_full = paste(kingdom, phylum, class, order, family, genus, species, sep = ";"),
      lineage_to_order = paste(kingdom, phylum, class, order, sep = ";"),
      lineage_to_family = paste(kingdom, phylum, class, order, family, sep = ";"),
      lineage_to_genus = paste(kingdom, phylum, class, order, family, genus, sep = ";")
    )
}

#####################################################
# Start
#####################################################
log_msg("Starting R postprocessing script...")
log_msg(paste("Using input directory:", EUK_DIR))

OUTPUT_DIR <- "taxonomic-profiling/eukaryotes/processed_data_may2026"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

merged_species_file <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_species_level.tsv")
merged_genus_file   <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_genus_level.tsv")
merged_family_file  <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_family_level.tsv")

if (file.exists(merged_species_file) &&
    file.exists(merged_genus_file) &&
    file.exists(merged_family_file)) {
  log_msg("Merged outputs already exist; nothing to do. Exiting.")
  quit(save = "no", status = 0)
}

#####################################################
# Load taxonomy DB (only builds sqlite once)
#####################################################
sqlite_file <- "nameNode.sqlite"

if (!file.exists(sqlite_file)) {
  names_file <- file.path(TAX_PATH_NCBI, "taxdump", "names.dmp")
  nodes_file <- file.path(TAX_PATH_NCBI, "taxdump", "nodes.dmp")

  if (file.exists(names_file) && file.exists(nodes_file)) {
    log_msg("Importing taxonomy into SQLite...")
    read.names.sql(names_file, sqlFile = sqlite_file, overwrite = TRUE)
    read.nodes.sql(nodes_file, sqlFile = sqlite_file, overwrite = TRUE)
    log_msg("Taxonomy import complete.")
  } else {
    stop("Warning: taxonomy files missing!")
  }
} else {
  log_msg("SQLite taxonomy already exists. Skipping import.")
}

#####################################################
# Load sample list
#####################################################
sample_list_file <- SAMPLE_LIST
if (!file.exists(sample_list_file)) stop(paste("Sample list not found at", sample_list_file))
sample_list <- read.table(sample_list_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lib_ids <- sample_list[[1]]
log_msg(paste("Loaded sample list with", length(lib_ids), "library IDs"))

#####################################################
# Read new unicorn rank-level taxstats directly
#####################################################
stats_species_file <- file.path(OUTPUT_DIR, "unicorn_stats_species_level.tsv")
stats_genus_file   <- file.path(OUTPUT_DIR, "unicorn_stats_genus_level.tsv")
stats_family_file  <- file.path(OUTPUT_DIR, "unicorn_stats_family_level.tsv")

stats_species <- cache_read(stats_species_file, fread, sep = "\t", data.table = FALSE)
stats_genus   <- cache_read(stats_genus_file,   fread, sep = "\t", data.table = FALSE)
stats_family  <- cache_read(stats_family_file,  fread, sep = "\t", data.table = FALSE)

if (is.null(stats_species)) {
  log_msg("Species unicorn cache missing; reading per-library species taxstats ...")
  raw_species <- read_or_process(
    file.path(OUTPUT_DIR, "species_taxstats_all.tsv"),
    EUK_DIR,
    "\\.filtered\\.spec\\.taxstats$"
  )
  stats_species <- clean_taxstats_columns(raw_species, "species")
  stats_species <- cache_write(stats_species, stats_species_file, fwrite, sep = "\t", quote = FALSE)
} else {
  log_msg(paste("Loaded unicorn species stats rows:", nrow(stats_species)))
}

if (is.null(stats_genus)) {
  log_msg("Genus unicorn cache missing; reading per-library genus taxstats ...")
  raw_genus <- read_or_process(
    file.path(OUTPUT_DIR, "genus_taxstats_all.tsv"),
    EUK_DIR,
    "\\.filtered\\.genus\\.taxstats$"
  )
  stats_genus <- clean_taxstats_columns(raw_genus, "genus")
  stats_genus <- cache_write(stats_genus, stats_genus_file, fwrite, sep = "\t", quote = FALSE)
} else {
  log_msg(paste("Loaded unicorn genus stats rows:", nrow(stats_genus)))
}

if (is.null(stats_family)) {
  log_msg("Family unicorn cache missing; reading per-library family taxstats ...")
  raw_family <- read_or_process(
    file.path(OUTPUT_DIR, "family_taxstats_all.tsv"),
    EUK_DIR,
    "\\.filtered\\.family\\.taxstats$"
  )
  stats_family <- clean_taxstats_columns(raw_family, "family")
  stats_family <- cache_write(stats_family, stats_family_file, fwrite, sep = "\t", quote = FALSE)
} else {
  log_msg(paste("Loaded unicorn family stats rows:", nrow(stats_family)))
}


#####################################################
# metaDMG caches: keep extracting correct rank rows
# and compute CCC stats
#####################################################
metaDMG_species_file <- file.path(OUTPUT_DIR, "metaDMG_species_level.tsv")
metaDMG_genus_file   <- file.path(OUTPUT_DIR, "metaDMG_genus_level.tsv")
metaDMG_family_file  <- file.path(OUTPUT_DIR, "metaDMG_family_level.tsv")

metaDMG_species <- cache_read(metaDMG_species_file, fread, sep = "\t", data.table = FALSE)
metaDMG_genus   <- cache_read(metaDMG_genus_file,   fread, sep = "\t", data.table = FALSE)
metaDMG_family  <- cache_read(metaDMG_family_file,  fread, sep = "\t", data.table = FALSE)

if (is.null(metaDMG_species) || is.null(metaDMG_genus) || is.null(metaDMG_family)) {
  log_msg("metaDMG rank caches missing; loading/processing lca_all.tsv ...")

  lca_all_df <- read_or_process(
    file.path(OUTPUT_DIR, "lca_all.tsv"),
    EUK_DIR,
    "\\.sort\\.filtered\\.agg\\.stat\\.gz$"
  )
  log_msg(paste("LCA rows:", nrow(lca_all_df)))

  log_msg("Computing CCC damage-fit statistics for Eukaryota rows ...")
  lca_all_df <- add_metadmg_ccc_stats(
    df      = lca_all_df,
    samples = unique(lca_all_df$library_id),
    ci      = "asymptotic",
    nperm   = 100,
    nproc   = 14
  )
  log_msg("CCC statistics joined onto lca_all_df.")

  if (is.null(metaDMG_species)) {
    metaDMG_species <- lca_all_df |>
      filter(rank == "species") |>
      mutate(
        library_id = as.character(library_id),
        taxid = as.character(taxid),
        species = as.character(name)
      ) |>
      rename(total_n_reads_from_mD = nreads)

    metaDMG_species <- cache_write(metaDMG_species, metaDMG_species_file, fwrite, sep = "\t", quote = FALSE)
  }

  if (is.null(metaDMG_genus)) {
    metaDMG_genus <- lca_all_df |>
      filter(rank == "genus") |>
      mutate(
        library_id = as.character(library_id),
        taxid = as.character(taxid),
        genus = as.character(name)
      ) |>
      rename(total_n_reads_from_mD = nreads)

    metaDMG_genus <- cache_write(metaDMG_genus, metaDMG_genus_file, fwrite, sep = "\t", quote = FALSE)
  }

  if (is.null(metaDMG_family)) {
    metaDMG_family <- lca_all_df |>
      filter(rank == "family") |>
      mutate(
        library_id = as.character(library_id),
        taxid = as.character(taxid),
        family = as.character(name)
      ) |>
      rename(total_n_reads_from_mD = nreads)

    metaDMG_family <- cache_write(metaDMG_family, metaDMG_family_file, fwrite, sep = "\t", quote = FALSE)
  }
} else {
  log_msg(paste("Loaded metaDMG species rows:", nrow(metaDMG_species)))
  log_msg(paste("Loaded metaDMG genus rows:", nrow(metaDMG_genus)))
  log_msg(paste("Loaded metaDMG family rows:", nrow(metaDMG_family)))
}

#####################################################
# Taxonomy cache for all taxids
#####################################################
tax_data_file <- file.path(OUTPUT_DIR, "tax_data_all_taxids.tsv")
tax_data <- cache_read(tax_data_file, fread, sep = "\t", data.table = FALSE)

if (is.null(tax_data)) {
  log_msg("Taxonomy cache missing; querying SQLite taxonomy for union of taxids...")

  all_taxids <- unique(c(
    stats_species$taxid,
    stats_genus$taxid,
    stats_family$taxid,
    metaDMG_species$taxid,
    metaDMG_genus$taxid,
    metaDMG_family$taxid
  ))

  all_taxids <- all_taxids[!is.na(all_taxids) & nzchar(all_taxids)]

  tax_data <- getTaxonomy(
    all_taxids,
    sqlite_file,
    desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  ) |>
    as_tibble() |>
    mutate(taxid = all_taxids) |>
    add_lineage_paths()

  tax_data <- cache_write(tax_data, tax_data_file, fwrite, sep = "\t", quote = FALSE)
} else {
  tax_data <- as_tibble(tax_data)
  log_msg(paste("Loaded taxonomy rows:", nrow(tax_data)))
}

tax_data <- tax_data |>
  mutate(taxid = readr::parse_integer(as.character(taxid)))

#####################################################
# Add taxonomy to unicorn tables
#####################################################
stats_species_tax <- stats_species |>
  mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  left_join(tax_data, by = "taxid")

stats_genus_tax <- stats_genus |>
  mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  left_join(tax_data, by = "taxid")

stats_family_tax <- stats_family |>
  mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  left_join(tax_data, by = "taxid")

#####################################################
# Add taxonomy to metaDMG tables
#####################################################
metaDMG_species_tax <- metaDMG_species |>
  mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  left_join(tax_data, by = "taxid")

metaDMG_genus_tax <- metaDMG_genus |>
  mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  left_join(tax_data, by = "taxid")

metaDMG_family_tax <- metaDMG_family |>
  mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  left_join(tax_data, by = "taxid")

#####################################################
# Merge unicorn + metaDMG
#####################################################
bf_md_species <- stats_species_tax |>
  mutate(
    library_id = as.character(library_id),
    taxid = as.integer(taxid)
  ) |>
  inner_join(
    metaDMG_species_tax |>
      mutate(
        library_id = as.character(library_id),
        taxid = as.integer(taxid)
      ),
    by = c("library_id", "taxid"),
    suffix = c("_unicorn", "_metaDMG")
  )

bf_md_genus <- stats_genus_tax |>
  mutate(
    library_id = as.character(library_id),
    taxid = as.integer(taxid)
  ) |>
  inner_join(
    metaDMG_genus_tax |>
      mutate(
        library_id = as.character(library_id),
        taxid = as.integer(taxid)
      ),
    by = c("library_id", "taxid"),
    suffix = c("_unicorn", "_metaDMG")
  )

bf_md_family <- stats_family_tax |>
  mutate(
    library_id = as.character(library_id),
    taxid = as.integer(taxid)
  ) |>
  inner_join(
    metaDMG_family_tax |>
      mutate(
        library_id = as.character(library_id),
        taxid = as.integer(taxid)
      ),
    by = c("library_id", "taxid"),
    suffix = c("_unicorn", "_metaDMG")
  )

#####################################################
# Write merged outputs
#####################################################
if (!file.exists(merged_species_file) || file.info(merged_species_file)$size == 0) {
  fwrite(bf_md_species, merged_species_file, sep = "\t", quote = FALSE)
  log_msg(paste("Wrote species merged file:", merged_species_file))
} else {
  log_msg(paste("Species merged file already exists:", merged_species_file))
}

if (!file.exists(merged_genus_file) || file.info(merged_genus_file)$size == 0) {
  fwrite(bf_md_genus, merged_genus_file, sep = "\t", quote = FALSE)
  log_msg(paste("Wrote genus merged file:", merged_genus_file))
} else {
  log_msg(paste("Genus merged file already exists:", merged_genus_file))
}

if (!file.exists(merged_family_file) || file.info(merged_family_file)$size == 0) {
  fwrite(bf_md_family, merged_family_file, sep = "\t", quote = FALSE)
  log_msg(paste("Wrote family merged file:", merged_family_file))
} else {
  log_msg(paste("Family merged file already exists:", merged_family_file))
}

log_msg("Finished R postprocessing for unicorn + metaDMG analysis.")
