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

#####################################################
# Logging helper
#####################################################
log_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

#####################################################
# Working directory and shared functions
#####################################################
setwd("/projects/caeg/people/ngm902/apps/repos/rocs")
source("code/eukaryotes/butteracid/postprocessing_functions.R")

#####################################################
# Define variables from config
#####################################################
TAX_PATH_NCBI <- "/datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/"
EUK_DIR <- "taxonomic-profiling/eukaryotes/unicorn_rerun_new_metadmg"
OUTPUT_DIR <- "taxonomic-profiling/eukaryotes/processed_data_june2026_incremental"
PER_SAMPLE_DIR <- file.path(OUTPUT_DIR, "per_sample")
PER_SAMPLE_SPECIES_DIR <- file.path(PER_SAMPLE_DIR, "species")
PER_SAMPLE_GENUS_DIR <- file.path(PER_SAMPLE_DIR, "genus")
PER_SAMPLE_FAMILY_DIR <- file.path(PER_SAMPLE_DIR, "family")
STATE_DIR <- file.path(OUTPUT_DIR, "state")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PER_SAMPLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PER_SAMPLE_SPECIES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PER_SAMPLE_GENUS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PER_SAMPLE_FAMILY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STATE_DIR, recursive = TRUE, showWarnings = FALSE)

merged_species_file <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_species_level.tsv")
merged_genus_file <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_genus_level.tsv")
merged_family_file <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_family_level.tsv")
processed_samples_file <- file.path(STATE_DIR, "processed_samples.tsv")
sqlite_file <- file.path(STATE_DIR, "nameNode.sqlite")
# SAMPLE_WORKERS <- suppressWarnings(as.integer(Sys.getenv("BUTTERACID_SAMPLE_WORKERS", "1")))
# if (is.na(SAMPLE_WORKERS) || SAMPLE_WORKERS < 1) SAMPLE_WORKERS <- 1L
SAMPLE_WORKERS <- 10


# In the incremental workflow, CCC runs on a single sample at a time, so
# nested parallelism adds overhead and can hang on HPC. Keep CCC single-threaded.
# CCC_NPROC <- suppressWarnings(as.integer(Sys.getenv("BUTTERACID_CCC_NPROC", "1")))
# if (is.na(CCC_NPROC) || CCC_NPROC < 1) CCC_NPROC <- 1L
# if (SAMPLE_WORKERS > 1 && CCC_NPROC > 1) {
#   log_msg("Forcing CCC_NPROC to 1 because sample-level parallelism is enabled")
#   CCC_NPROC <- 1L
# }
CCC_NPROC <- 24


#####################################################
# Cache helpers
#####################################################
cache_read <- function(path, read_fun, ...) {
  if (file.exists(path) && file.info(path)$size > 0) {
    log_msg("Cache hit: ", path)
    return(read_fun(path, ...))
  }
  NULL
}

cache_write <- function(df, path, write_fun, ...) {
  write_fun(df, path, ...)
  log_msg("Wrote: ", path)
  df
}

#####################################################
# File and table helpers
#####################################################
extract_library_id <- function(file) {
  b <- basename(file)
  b <- sub("\\.(sort|filtered)\\..*$", "", b)
  normalize_library_id(b)
}

normalize_library_id <- function(x) {
  vapply(x, function(id) {
    parts <- strsplit(as.character(id), "_", fixed = TRUE)[[1]]
    if (length(parts) >= 2 && identical(parts[1], parts[2])) {
      parts[1]
    } else {
      as.character(id)
    }
  }, character(1), USE.NAMES = FALSE)
}

file_nonempty <- function(path) {
  file.exists(path) && file.info(path)$size > 0
}

all_files_nonempty <- function(paths) {
  all(vapply(paths, file_nonempty, logical(1)))
}

sample_input_paths <- function(sample_id) {
  list(
    species = file.path(EUK_DIR, paste0(sample_id, "_collapsed.filtered.spec.taxstats")),
    genus = file.path(EUK_DIR, paste0(sample_id, "_collapsed.filtered.genus.taxstats")),
    family = file.path(EUK_DIR, paste0(sample_id, "_collapsed.filtered.family.taxstats")),
    lca = file.path(EUK_DIR, paste0(sample_id, "_collapsed.sort.filtered.agg.stat.gz"))
  )
}

sample_output_paths <- function(sample_id) {
  list(
    species = file.path(PER_SAMPLE_SPECIES_DIR, paste0(sample_id, ".species.tsv")),
    genus = file.path(PER_SAMPLE_GENUS_DIR, paste0(sample_id, ".genus.tsv")),
    family = file.path(PER_SAMPLE_FAMILY_DIR, paste0(sample_id, ".family.tsv"))
  )
}

sample_inputs_ready <- function(sample_id) {
  all_files_nonempty(unlist(sample_input_paths(sample_id), use.names = FALSE))
}

sample_outputs_ready <- function(sample_id) {
  all_files_nonempty(unlist(sample_output_paths(sample_id), use.names = FALSE))
}

read_with_id <- function(file) {
  if (!file_nonempty(file)) {
    log_msg("Skipping missing or empty file: ", file)
    return(NULL)
  }

  lib_id <- extract_library_id(file)

  fread(file, colClasses = "character", data.table = FALSE) |>
    as_tibble() |>
    mutate(library_id = lib_id)
}

clean_taxstats_columns <- function(df, rank_name) {
  rank_col <- switch(
    rank_name,
    species = "species",
    genus = "genus",
    family = "family",
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

read_taxstats_rank <- function(path, rank_name) {
  df <- read_with_id(path)
  if (is.null(df) || nrow(df) == 0) {
    return(tibble())
  }
  clean_taxstats_columns(df, rank_name)
}

read_metadmg_rank_tables <- function(lca_file, sample_id) {
  lca_df <- read_with_id(lca_file)
  if (is.null(lca_df) || nrow(lca_df) == 0) {
    return(list(
      species = tibble(),
      genus = tibble(),
      family = tibble()
    ))
  }

  lca_df <- add_metadmg_ccc_stats(
    df = lca_df,
    samples = normalize_library_id(sample_id),
    ci = "asymptotic",
    nperm = 100,
    nproc = CCC_NPROC
  )

  make_rank_df <- function(rank_name, value_col) {
    lca_df |>
      filter(rank == rank_name) |>
      mutate(
        library_id = as.character(library_id),
        taxid = as.character(taxid),
        !!rank_name := as.character(name)
      ) |>
      rename(total_n_reads_from_mD = nreads)
  }

  list(
    species = make_rank_df("species", "species"),
    genus = make_rank_df("genus", "genus"),
    family = make_rank_df("family", "family")
  )
}

get_taxonomy_for_taxids <- function(taxids) {
  taxids <- unique(as.character(taxids))
  taxids <- taxids[!is.na(taxids) & nzchar(taxids)]

  if (length(taxids) == 0) {
    return(tibble(
      taxid = integer(),
      kingdom = character(),
      phylum = character(),
      class = character(),
      order = character(),
      family = character(),
      genus = character(),
      species = character(),
      lineage_full = character(),
      lineage_to_order = character(),
      lineage_to_family = character(),
      lineage_to_genus = character()
    ))
  }

  getTaxonomy(
    taxids,
    sqlite_file,
    desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  ) |>
    as_tibble() |>
    mutate(taxid = taxids) |>
    add_lineage_paths() |>
    mutate(taxid = readr::parse_integer(as.character(taxid)))
}

join_taxonomy <- function(df, tax_data) {
  if (nrow(df) == 0) {
    return(df)
  }

  df |>
    mutate(taxid = readr::parse_integer(as.character(taxid))) |>
    left_join(tax_data, by = "taxid")
}

merge_unicorn_metadmg <- function(unicorn_df, metadmg_df) {
  if (nrow(unicorn_df) == 0 || nrow(metadmg_df) == 0) {
    return(tibble())
  }

  unicorn_df |>
    mutate(
      library_id = as.character(library_id),
      taxid = as.integer(taxid)
    ) |>
    inner_join(
      metadmg_df |>
        mutate(
          library_id = as.character(library_id),
          taxid = as.integer(taxid)
        ),
      by = c("library_id", "taxid"),
      suffix = c("_unicorn", "_metaDMG")
    )
}

write_processed_samples <- function(processed_samples) {
  df <- tibble(library_id = sort(unique(processed_samples)))
  fwrite(df, processed_samples_file, sep = "\t", quote = FALSE)
  log_msg("Updated processed sample registry: ", processed_samples_file)
}

process_sample <- function(sample_id) {
  input_paths <- sample_input_paths(sample_id)
  output_paths <- sample_output_paths(sample_id)

  if (!sample_inputs_ready(sample_id)) {
    log_msg("Skipping sample with missing inputs: ", sample_id)
    return(FALSE)
  }

  if (sample_outputs_ready(sample_id)) {
    log_msg("Skipping already processed sample: ", sample_id)
    return(TRUE)
  }

  log_msg("Processing sample: ", sample_id)

  stats_species <- read_taxstats_rank(input_paths$species, "species")
  stats_genus <- read_taxstats_rank(input_paths$genus, "genus")
  stats_family <- read_taxstats_rank(input_paths$family, "family")

  metadmg_tables <- read_metadmg_rank_tables(input_paths$lca, sample_id)

  all_taxids <- unique(c(
    stats_species$taxid,
    stats_genus$taxid,
    stats_family$taxid,
    metadmg_tables$species$taxid,
    metadmg_tables$genus$taxid,
    metadmg_tables$family$taxid
  ))

  tax_data <- get_taxonomy_for_taxids(all_taxids)

  species_merged <- merge_unicorn_metadmg(
    join_taxonomy(stats_species, tax_data),
    join_taxonomy(metadmg_tables$species, tax_data)
  )
  genus_merged <- merge_unicorn_metadmg(
    join_taxonomy(stats_genus, tax_data),
    join_taxonomy(metadmg_tables$genus, tax_data)
  )
  family_merged <- merge_unicorn_metadmg(
    join_taxonomy(stats_family, tax_data),
    join_taxonomy(metadmg_tables$family, tax_data)
  )

  cache_write(species_merged, output_paths$species, fwrite, sep = "\t", quote = FALSE)
  cache_write(genus_merged, output_paths$genus, fwrite, sep = "\t", quote = FALSE)
  cache_write(family_merged, output_paths$family, fwrite, sep = "\t", quote = FALSE)

  TRUE
}

combine_rank_outputs <- function(rank_dir, output_file) {
  files <- list.files(rank_dir, pattern = "\\.tsv$", full.names = TRUE)

  if (length(files) == 0) {
    log_msg("No per-sample files found in ", rank_dir, ". Writing empty table to ", output_file)
    fwrite(tibble(), output_file, sep = "\t", quote = FALSE)
    return(invisible(tibble()))
  }

  combined <- map_dfr(
    files,
    ~ fread(.x, sep = "\t", data.table = FALSE) |> as_tibble()
  )

  fwrite(combined, output_file, sep = "\t", quote = FALSE)
  log_msg("Wrote combined file: ", output_file, " (rows: ", nrow(combined), ")")
  invisible(combined)
}

#####################################################
# Build sample list from available workflow outputs
#####################################################
files <- list.files(
  path = EUK_DIR,
  pattern = "_collapsed\\.sort\\.filtered\\.agg\\.stat\\.gz$",
  full.names = FALSE
)
if (length(files) == 0) {
  stop(sprintf("No files found in %s matching pattern '*_collapsed.sort.filtered.agg.stat.gz'", EUK_DIR))
}

samples <- sub("_collapsed\\.sort\\.filtered\\.agg\\.stat\\.gz$", "", files)
SAMPLE_LIST <- file.path(OUTPUT_DIR, "sample_list.tsv")
write.table(
  data.frame(library_id = samples),
  SAMPLE_LIST,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
log_msg("Created sample list: ", SAMPLE_LIST, " with ", length(samples), " samples")

#####################################################
# Load taxonomy DB (only builds sqlite once)
#####################################################
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
# Main processing loop
#####################################################
sample_list <- read.table(SAMPLE_LIST, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lib_ids <- unique(sample_list[[1]])
log_msg("Loaded sample list with ", length(lib_ids), " library IDs")

processed_registry <- cache_read(processed_samples_file, fread, sep = "\t", data.table = FALSE)
if (is.null(processed_registry)) {
  processed_ids <- character()
} else {
  processed_ids <- unique(normalize_library_id(as.character(processed_registry[[1]])))
}

successful_samples <- processed_ids
newly_processed <- character()

log_msg("Sample workers: ", SAMPLE_WORKERS, " | CCC nproc per sample: ", CCC_NPROC)

process_one <- function(sample_id) {
  sample_id_short <- normalize_library_id(sample_id)
  result <- process_sample(sample_id)
  tibble(
    sample_id = as.character(sample_id),
    sample_id_short = as.character(sample_id_short),
    ok = isTRUE(result)
  )
}

if (SAMPLE_WORKERS > 1) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = SAMPLE_WORKERS)
  results <- furrr::future_map_dfr(
    lib_ids,
    process_one,
    .options = furrr::furrr_options(seed = TRUE)
  )
} else {
  results <- purrr::map_dfr(lib_ids, process_one)
}

successful_run <- results |>
  filter(ok) |>
  pull(sample_id_short) |>
  unique()

successful_samples <- unique(c(successful_samples, successful_run))
newly_processed <- setdiff(successful_run, processed_ids)

write_processed_samples(successful_samples)
log_msg("Newly processed samples in this run: ", length(newly_processed))

#####################################################
# Rebuild combined outputs from per-sample outputs
#####################################################
combine_rank_outputs(PER_SAMPLE_SPECIES_DIR, merged_species_file)
combine_rank_outputs(PER_SAMPLE_GENUS_DIR, merged_genus_file)
combine_rank_outputs(PER_SAMPLE_FAMILY_DIR, merged_family_file)

log_msg("Finished incremental R postprocessing for unicorn + metaDMG analysis.")
