#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))

meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
me_main <- fread(file.path(DIRS$main, "module_eigengenes_main.tsv"))
xrf <- fread(INPUTS$xrf)

me_main[, label := sample]
me_main[, sample := NULL]

dt <- merge(meta_main, me_main, by = "label", all.x = FALSE)
dt <- dt[core %in% PARAMS$training_cores]
dt[, core_group := core_group(core)]

me_cols <- setdiff(grep("^ME", names(dt), value = TRUE), "MEgrey")
me_scaled <- corewise_z_transform(
  mat = as.matrix(dt[, ..me_cols]),
  core_ids = dt$core,
  feature_names = me_cols,
  variant = "corewise_z"
)
dt[, (me_cols) := as.data.table(me_scaled$scaled)]

xrf[, core_group := core_group(core)]
xrf_elements <- c("ba", "ti", "ca", "zr", "rb", "sr", "al", "fe", "mn", "br", "p", "si")
xrf_missing <- setdiff(xrf_elements, names(xrf))
if (length(xrf_missing) > 0) {
  stop(sprintf("XRF input is missing expected element columns: %s", paste(xrf_missing, collapse = ", ")))
}
xrf_lookup <- xrf[, c("core_group", "depth_in_core_cm", xrf_elements), with = FALSE]

match_xrf_means <- function(samples_dt, xrf_dt, xrf_elements, half_width_cm) {
  required_cols <- c("label", "core", "core_group", "y_bp", "depth_in_core_cm")
  missing_cols <- setdiff(required_cols, names(samples_dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("samples_dt is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  valid <- samples_dt[
    !is.na(core_group) & nzchar(core_group) &
      is.finite(depth_in_core_cm)
  ]

  id_cols <- c("label", "core", "core_group", "y_bp", "depth_in_core_cm")
  out_cols <- c(id_cols, xrf_elements, "match_n_points", "match_quality", "match_min_abs_depth_cm")
  if (nrow(valid) == 0) {
    out <- data.table::as.data.table(setNames(replicate(length(out_cols), logical(0), simplify = FALSE), out_cols))
    return(out)
  }

  rows <- vector("list", nrow(valid))
  for (i in seq_len(nrow(valid))) {
    r <- valid[i]
    xsub <- xrf_dt[core_group == r$core_group]
    if (nrow(xsub) == 0) next

    d0 <- r$depth_in_core_cm
    in_win <- xsub[depth_in_core_cm >= (d0 - half_width_cm) & depth_in_core_cm <= (d0 + half_width_cm)]
    if (nrow(in_win) > 0) {
      matched <- in_win
      q <- "interval"
      dist <- 0
    } else {
      xsub[, abs_dist := abs(depth_in_core_cm - d0)]
      matched <- xsub[which.min(abs_dist)]
      q <- "fallback_nearest"
      dist <- matched$abs_dist[1]
      matched[, abs_dist := NULL]
    }

    s <- matched[, lapply(.SD, function(v) {
      m <- mean(v, na.rm = TRUE)
      if (is.nan(m)) NA_real_ else m
    }), .SDcols = xrf_elements]
    s[, `:=`(
      label = r$label,
      core = r$core,
      core_group = r$core_group,
      y_bp = r$y_bp,
      depth_in_core_cm = r$depth_in_core_cm,
      match_n_points = nrow(matched),
      match_quality = q,
      match_min_abs_depth_cm = dist
    )]
    rows[[i]] <- s
  }

  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0) {
    out <- data.table::as.data.table(setNames(replicate(length(out_cols), logical(0), simplify = FALSE), out_cols))
  }
  out
}

ratio_log_z <- function(dt, num, den, out_prefix) {
  raw_col <- paste0(out_prefix, "_ratio")
  log_col <- paste0(out_prefix, "_log")
  z_col <- paste0(out_prefix, "_z")
  dt[, (raw_col) := fifelse(get(num) > 0 & get(den) > 0, get(num) / get(den), NA_real_)]
  dt[, (log_col) := fifelse(get(raw_col) > 0, log(get(raw_col)), NA_real_)]
  z_obj <- corewise_z_transform(
    mat = matrix(dt[[log_col]], ncol = 1),
    core_ids = dt$core,
    feature_names = log_col,
    variant = "corewise_z"
  )
  dt[, (z_col) := as.numeric(z_obj$scaled[, 1])]
}

xrf_match <- match_xrf_means(
  samples_dt = dt,
  xrf_dt = xrf_lookup,
  xrf_elements = xrf_elements,
  half_width_cm = PARAMS$xrf_half_width_cm
)
dtx <- merge(dt, xrf_match, by = c("label", "core", "core_group", "y_bp", "depth_in_core_cm"), all.x = TRUE)

ratio_log_z(dtx, "ba", "ti", "ba_ti")
ratio_log_z(dtx, "ca", "ti", "ca_ti")
ratio_log_z(dtx, "zr", "rb", "zr_rb")
ratio_log_z(dtx, "ti", "al", "ti_al")
ratio_log_z(dtx, "mn", "fe", "mn_fe")

meta_all <- fread(INPUTS$metadata)
meta_all[, core_group := core_group(core)]
meta_all_valid <- meta_all[!is.na(core) & nzchar(core) & is.finite(depth_in_core_cm)]

meta_all_xrf_match <- match_xrf_means(
  samples_dt = meta_all_valid,
  xrf_dt = xrf_lookup,
  xrf_elements = xrf_elements,
  half_width_cm = PARAMS$xrf_half_width_cm
)

meta_all_xrf <- merge(
  meta_all,
  meta_all_xrf_match,
  by = c("label", "core", "core_group", "y_bp", "depth_in_core_cm"),
  all.x = TRUE
)

ratio_log_z(meta_all_xrf, "ba", "ti", "ba_ti")
ratio_log_z(meta_all_xrf, "ca", "ti", "ca_ti")
ratio_log_z(meta_all_xrf, "zr", "rb", "zr_rb")
ratio_log_z(meta_all_xrf, "ti", "al", "ti_al")
ratio_log_z(meta_all_xrf, "mn", "fe", "mn_fe")

capture_fit <- function(expr) {
  warn <- character()
  res <- withCallingHandlers(tryCatch(expr, error = function(e) e), warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(result = res, ok = !inherits(res, "error"), error = if (inherits(res, "error")) conditionMessage(res) else NA_character_, warnings = paste(unique(warn), collapse = " | "))
}

proxy_panel <- PARAMS$core_proxy_panel
rows_res <- list()
diag_rows <- list()

for (me in me_cols) {
  for (px in proxy_panel) {
    dd <- dtx[is.finite(get(me)) & is.finite(get(px))]
    if (nrow(dd) < 30 || uniqueN(dd$core) < 2) next
    dd[, age_kyr := y_bp / 1000]
    dd[, core_f := factor(core)]
    fml <- as.formula(sprintf("%s ~ %s + age_kyr + core_f", me, px))
    lm_res <- capture_fit(lm(fml, data = dd))
    # age_kyr is continuous (non-integer) and irregularly spaced by design;
    # use continuous-time AR(1) rather than corAR1/corARMA (integer index required).
    gls_res <- capture_fit(gls(model = fml, data = dd, correlation = corCAR1(form = ~ age_kyr | core_f), method = "REML"))

    diag_rows[[length(diag_rows) + 1]] <- data.table(module = me, proxy = px, n = nrow(dd), lm_ok = lm_res$ok, gls_ok = gls_res$ok, lm_error = lm_res$error, gls_error = gls_res$error)
    if (!lm_res$ok) next
    cf <- summary(lm_res$result)$coefficients
    if (!(px %in% rownames(cf))) next
    beta <- cf[px, "Estimate"]
    p <- cf[px, "Pr(>|t|)"]
    beta_ar1 <- p_ar1 <- NA_real_
    if (gls_res$ok) {
      c2 <- summary(gls_res$result)$tTable
      if (px %in% rownames(c2)) {
        beta_ar1 <- c2[px, "Value"]
        p_ar1 <- c2[px, "p-value"]
      }
    }
    rows_res[[length(rows_res) + 1]] <- data.table(module = me, proxy = px, n = nrow(dd), beta = beta, p_value = p, beta_ar1 = beta_ar1, p_value_ar1 = p_ar1)
  }
}

mod_xrf <- if (length(rows_res) > 0) rbindlist(rows_res, fill = TRUE) else data.table()
if (nrow(mod_xrf) > 0) {
  mod_xrf[, q_value := p.adjust(p_value, method = "BH")]
  mod_xrf[, q_value_ar1 := p.adjust(p_value_ar1, method = "BH")]
}

# LOCO
loco_rows <- list()
if (nrow(mod_xrf) > 0) {
  for (i in seq_len(nrow(mod_xrf))) {
    me <- mod_xrf$module[i]
    px <- mod_xrf$proxy[i]
    dd0 <- dtx[is.finite(get(me)) & is.finite(get(px))]
    if (nrow(dd0) < 30 || uniqueN(dd0$core) < 2) next
    dd0[, age_kyr := y_bp / 1000]
    dd0[, core_f := factor(core)]
    fml <- as.formula(sprintf("%s ~ %s + age_kyr + core_f", me, px))
    fit_full <- tryCatch(lm(fml, data = dd0), error = function(e) NULL)
    if (is.null(fit_full)) next
    cf_full <- summary(fit_full)$coefficients
    if (!(px %in% rownames(cf_full))) next
    out <- data.table(module = me, proxy = px, beta_full = cf_full[px, "Estimate"], p_full = cf_full[px, "Pr(>|t|)"])
    for (lg in PARAMS$training_cores) {
      ddk <- dd0[core != lg]
      ddk[, core_f := factor(core)]
      if (nrow(ddk) < 30 || uniqueN(ddk$core) < 2) {
        out[, paste0("beta_leave_", lg) := NA_real_]
        out[, paste0("p_leave_", lg) := NA_real_]
        next
      }
      fk <- tryCatch(lm(fml, data = ddk), error = function(e) NULL)
      if (is.null(fk)) {
        out[, paste0("beta_leave_", lg) := NA_real_]
        out[, paste0("p_leave_", lg) := NA_real_]
        next
      }
      cfk <- summary(fk)$coefficients
      out[, paste0("beta_leave_", lg) := if (px %in% rownames(cfk)) cfk[px, "Estimate"] else NA_real_]
      out[, paste0("p_leave_", lg) := if (px %in% rownames(cfk)) cfk[px, "Pr(>|t|)"] else NA_real_]
    }
    loco_rows[[length(loco_rows) + 1]] <- out
  }
}

loco <- if (length(loco_rows) > 0) rbindlist(loco_rows, fill = TRUE) else data.table()
if (nrow(loco) > 0) {
  loco[, q_full := p.adjust(p_full, method = "BH")]
  loco[, sign_preserved_all := mapply(function(bf, b1, b2, b3) {
    b <- c(b1, b2, b3); b <- b[is.finite(b)]
    if (!is.finite(bf) || bf == 0 || length(b) == 0) return(FALSE)
    all(sign(b) == sign(bf))
  }, beta_full, beta_leave_ST8, beta_leave_ST13, beta_leave_GeoB25202_R1)]
  loco[, all_leave_p_lt_010 := mapply(function(p1, p2, p3) {
    p <- c(p1, p2, p3); p <- p[is.finite(p)]
    length(p) > 0 && all(p < PARAMS$q_threshold)
  }, p_leave_ST8, p_leave_ST13, p_leave_GeoB25202_R1)]
  loco[, stability_flag := fifelse(q_full < PARAMS$q_threshold & sign_preserved_all & all_leave_p_lt_010, "stable", "unstable_or_mixed")]
}

fwrite(dtx, file.path(DIRS$main, "module_xrf_dataset_training.tsv"), sep = "\t")
fwrite(meta_all_xrf, file.path(DIRS$main, "metadata_xrf_dataset_all_samples.tsv"), sep = "\t")
fwrite(mod_xrf, file.path(DIRS$main, "module_xrf_results.tsv"), sep = "\t")
fwrite(mod_xrf[q_value < PARAMS$q_threshold][order(q_value, -abs(beta))], file.path(DIRS$main, "module_xrf_significant.tsv"), sep = "\t")
fwrite(if (length(diag_rows) > 0) rbindlist(diag_rows, fill = TRUE) else data.table(), file.path(DIRS$main, "module_xrf_model_diagnostics.tsv"), sep = "\t")
fwrite(loco, file.path(DIRS$main, "loco_module_xrf.tsv"), sep = "\t")
fwrite(loco[stability_flag == "stable"][order(q_full, -abs(beta_full))], file.path(DIRS$main, "loco_stable_associations.tsv"), sep = "\t")

write_run_metadata(
  file.path(DIRS$main, "04_xrf_loco_main_run_metadata.tsv"),
  "04_xrf_loco_main.R",
  extra = list(
    xrf_half_width_cm = PARAMS$xrf_half_width_cm,
    q_threshold = PARAMS$q_threshold,
    training_cores = PARAMS$training_cores,
    harmonization = "corewise_z"
  )
)
write_session_info(file.path(DIRS$main, "04_xrf_loco_main_sessionInfo.txt"))

log_msg("04_xrf_loco_main complete.")
