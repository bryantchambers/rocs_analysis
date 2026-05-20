#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# Interactive-friendly: keep an explicit repo root, consistent with this project style.
repo_root <- "/maps/projects/caeg/people/ngm902/apps/repos/rocs"
if (dir.exists(repo_root)) {
  setwd(repo_root)
}

metadata_path <- "data/metadata_v5.tsv"
xrf_path <- "data/combined_xrf_not_normalized.tsv"
out_path <- "results/common/metadata_xrf_matched.tsv"

xrf_elements <- c("ba", "ti", "ca", "zr", "rb", "sr", "al", "fe", "mn", "br", "p", "si")
ratio_specs <- list(
  c("ba", "ti", "ba_ti"),
  c("ca", "ti", "ca_ti"),
  c("zr", "rb", "zr_rb"),
  c("ti", "al", "ti_al"),
  c("mn", "fe", "mn_fe")
)

required_metadata_cols <- c("core", "depth_in_core_cm")
required_xrf_cols <- c("core", "depth_in_core_cm", xrf_elements)

assert_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required file(s): %s", paste(missing, collapse = ", ")))
  }
}

assert_required_cols <- function(dt, required, dt_name) {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required column(s): %s", dt_name, paste(missing, collapse = ", ")))
  }
}

core_group_from_core <- function(x) {
  sub("_R[12]$", "", x)
}

corewise_z <- function(values, cores) {
  z <- rep(NA_real_, length(values))
  for (cc in unique(cores)) {
    idx <- which(cores == cc & is.finite(values))
    if (length(idx) < 2) next
    s <- sd(values[idx])
    if (!is.finite(s) || s == 0) next
    z[idx] <- (values[idx] - mean(values[idx])) / s
  }
  z
}

match_xrf_means <- function(samples_dt, xrf_dt, xrf_elements, half_width_cm = 0.5) {
  out <- copy(samples_dt)

  for (col in xrf_elements) out[, (col) := NA_real_]
  out[, `:=`(
    match_n_points = as.integer(NA),
    match_quality = NA_character_,
    match_min_abs_depth_cm = NA_real_
  )]

  for (i in seq_len(nrow(out))) {
    cg <- out$core_group[i]
    d0 <- out$depth_in_core_cm[i]

    if (is.na(cg) || !nzchar(cg) || !is.finite(d0)) next

    xsub <- xrf_dt[core_group == cg & is.finite(depth_in_core_cm)]
    if (nrow(xsub) == 0) next

    in_window <- xsub[depth_in_core_cm >= (d0 - half_width_cm) & depth_in_core_cm <= (d0 + half_width_cm)]

    if (nrow(in_window) > 0) {
      matched <- in_window
      quality <- "interval"
      min_abs_depth_cm <- 0
    } else {
      xsub[, abs_dist := abs(depth_in_core_cm - d0)]
      matched <- xsub[which.min(abs_dist)]
      quality <- "fallback_nearest"
      min_abs_depth_cm <- matched$abs_dist[1]
      matched[, abs_dist := NULL]
    }

    means <- matched[, lapply(.SD, function(v) {
      m <- mean(v, na.rm = TRUE)
      if (is.nan(m)) NA_real_ else m
    }), .SDcols = xrf_elements]

    out[i, (xrf_elements) := as.list(means[1])]
    out[i, `:=`(
      match_n_points = nrow(matched),
      match_quality = quality,
      match_min_abs_depth_cm = min_abs_depth_cm
    )]
  }

  out
}

add_ratio_log_z <- function(dt, num, den, out_prefix) {
  ratio_col <- paste0(out_prefix, "_ratio")
  log_col <- paste0(out_prefix, "_log")
  z_col <- paste0(out_prefix, "_z")

  dt[, (ratio_col) := fifelse(get(num) > 0 & get(den) > 0, get(num) / get(den), NA_real_)]
  dt[, (log_col) := fifelse(get(ratio_col) > 0, log(get(ratio_col)), NA_real_)]
  dt[, (z_col) := corewise_z(get(log_col), core)]
}

assert_files_exist(c(metadata_path, xrf_path))

metadata <- fread(metadata_path)
xrf <- fread(xrf_path)

assert_required_cols(metadata, required_metadata_cols, "metadata_v5")
assert_required_cols(xrf, required_xrf_cols, "combined_xrf_not_normalized")

metadata[, core_group := core_group_from_core(core)]
xrf[, core_group := core_group_from_core(core)]

xrf_lookup <- xrf[, c("core_group", "depth_in_core_cm", xrf_elements), with = FALSE]

metadata_xrf <- match_xrf_means(
  samples_dt = metadata,
  xrf_dt = xrf_lookup,
  xrf_elements = xrf_elements,
  half_width_cm = 0.5
)

for (spec in ratio_specs) {
  add_ratio_log_z(metadata_xrf, num = spec[1], den = spec[2], out_prefix = spec[3])
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(metadata_xrf, out_path, sep = "\t", na = "NA")

message(sprintf("Saved matched metadata + XRF to: %s", out_path))
