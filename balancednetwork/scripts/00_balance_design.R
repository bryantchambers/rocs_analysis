#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

source(here("config.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default = NA_character_) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[1]])
}

bin_width_kyr <- as.numeric(arg_val("bin_width_kyr", "10"))
force <- as.integer(arg_val("force", "0")) == 1L
if (!is.finite(bin_width_kyr) || bin_width_kyr <= 0) stop("Invalid bin width.")

OUT_BASE <- here("balancednetwork", "results")
OUT_TAB <- file.path(OUT_BASE, "tables")
OUT_LOG <- file.path(OUT_BASE, "logs")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOG, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(OUT_TAB, "balanced_baseline_samples.tsv")
if (file.exists(manifest_file) && !force) {
  message("[balancednetwork] balance design exists; skipping (use --force=1 to rebuild).")
  quit(save = "no")
}

vst <- readRDS(file.path(RESULTS$stage1, "prokaryotes_vst.rds"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))

train_cores <- PARAMS$stage1_cores
keep_meta <- meta[core %in% train_cores]
keep_meta[, age_kyr := y_bp / 1000]
keep_meta[, in_vst := label %in% rownames(vst)]

core_range <- keep_meta[in_vst == TRUE, .(
  age_min = min(age_kyr, na.rm = TRUE),
  age_max = max(age_kyr, na.rm = TRUE)
), by = core]
shared_age_min <- max(core_range$age_min)
shared_age_max <- min(core_range$age_max)

if (!is.finite(shared_age_min) || !is.finite(shared_age_max) || shared_age_max <= shared_age_min) {
  stop("Shared age window could not be computed for training cores.")
}

design_pool <- keep_meta[in_vst == TRUE & age_kyr >= shared_age_min & age_kyr <= shared_age_max]
design_pool[, age_bin_lo := floor(age_kyr / bin_width_kyr) * bin_width_kyr]
design_pool[, age_bin_hi := age_bin_lo + bin_width_kyr]
design_pool[, age_bin := sprintf("%.0f-%.0f", age_bin_lo, age_bin_hi)]

bin_counts <- design_pool[, .N, by = .(age_bin, age_bin_lo, age_bin_hi, core)]
setorder(bin_counts, age_bin_lo, core)

bin_wide <- dcast(
  bin_counts,
  age_bin + age_bin_lo + age_bin_hi ~ core,
  value.var = "N",
  fill = 0
)

for (cc in train_cores) {
  if (!cc %in% names(bin_wide)) bin_wide[, (cc) := 0L]
}

bin_wide[, quota_per_core := pmin(get(train_cores[1]), get(train_cores[2]), get(train_cores[3]))]
bin_wide[, retained := quota_per_core > 0]
retained_bins <- bin_wide[retained == TRUE]

if (nrow(retained_bins) == 0) stop("No age bins with non-zero quota across all training cores.")

set.seed(PARAMS$seed)
selected <- rbindlist(lapply(seq_len(nrow(retained_bins)), function(i) {
  b <- retained_bins[i]
  rbindlist(lapply(train_cores, function(core_id) {
    pool <- design_pool[core == core_id & age_bin == b$age_bin, .(label, core, age_kyr, age_bin)]
    take <- pool[sample.int(.N, size = b$quota_per_core, replace = FALSE)]
    take
  }))
}))

selected[, source := "balanced_baseline"]
setorder(selected, age_kyr, core, label)

selected_counts <- selected[, .N, by = core][order(core)]
per_core_n <- unique(selected_counts$N)
if (length(per_core_n) != 1L) stop("Balanced baseline selection failed core-count equality.")

all_train <- keep_meta[core %in% train_cores, .(label, core, y_bp, age_kyr, in_vst)]
all_train[, reason := fifelse(
  in_vst == FALSE, "missing_in_vst",
  fifelse(age_kyr < shared_age_min | age_kyr > shared_age_max, "outside_shared_age_window",
          "eligible_not_selected")
)]
all_train <- merge(all_train, selected[, .(label, selected = TRUE)], by = "label", all.x = TRUE)
all_train[is.na(selected), selected := FALSE]
all_train[selected == TRUE, reason := "selected"]
all_train[reason == "eligible_not_selected", reason := "not_selected_due_to_quota"]

summary_dt <- data.table(
  metric = c("bin_width_kyr", "shared_age_min_kyr", "shared_age_max_kyr", "retained_bins", "samples_per_core", "total_training_samples"),
  value = c(
    bin_width_kyr,
    shared_age_min,
    shared_age_max,
    nrow(retained_bins),
    per_core_n[1],
    nrow(selected)
  )
)

fwrite(bin_counts, file.path(OUT_TAB, "balance_bin_counts_long.tsv"), sep = "\t")
fwrite(bin_wide, file.path(OUT_TAB, "balance_bin_availability.tsv"), sep = "\t")
fwrite(retained_bins, file.path(OUT_TAB, "balance_bin_quotas.tsv"), sep = "\t")
fwrite(selected, manifest_file, sep = "\t")
fwrite(all_train, file.path(OUT_TAB, "balance_excluded_samples.tsv"), sep = "\t")
fwrite(summary_dt, file.path(OUT_TAB, "balance_design_summary.tsv"), sep = "\t")
fwrite(core_range, file.path(OUT_TAB, "balance_core_age_ranges.tsv"), sep = "\t")

message(sprintf(
  "[balancednetwork] design complete: bins=%d, per_core=%d, total=%d",
  nrow(retained_bins), per_core_n[1], nrow(selected)
))
