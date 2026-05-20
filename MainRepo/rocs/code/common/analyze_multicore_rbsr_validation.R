#!/usr/bin/env Rscript

# Multi-core validation of Rb/Sr-like association with library concentration.
# Guardrails:
# - observational, core-specific associations only (not causal proof)
# - downstream sequencing-output metrics are excluded as explanatory variables
# - ST5/ST8 use ST13-like layer aggregation + depth-aware models when feasible
# - GeoB25202 is intentionally basic and focused on rb_sr only

options(stringsAsFactors = FALSE)

suppressWarnings({
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required.")
  }
})

out_dir <- file.path("analysis", "multicore_rbsr_library_validation")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path("data", "metadata_v5.tsv")
xrf_path <- file.path("data", "combined_xrf_geochemistry_curated.csv")
foram_path <- file.path("data", "combined_foraminifera_geochem.tsv")
sst_path <- file.path("data", "combined_sst_proxies_separate_columns.tsv")

required_files <- c(metadata_path, xrf_path, foram_path, sst_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files (run from repo root): ", paste(missing_files, collapse = ", "))
}

params <- list(
  st13_like_cores = c("ST13", "ST5", "ST8"),
  geob_base_core = "GeoB25202",
  geob_replicate_regex = "_R[0-9]+$",
  sparse_match_tolerance_cm = 3,
  depth_spline_df = 4,
  screening_min_n = 10,
  model_min_n = 20,
  interaction_min_n = 18,
  stratum_min_n = 10,
  block_bootstrap_n = 600,
  support_p_threshold = 0.10
)

safe_min_positive <- function(x, default = 1e-6) {
  v <- x[is.finite(x) & x > 0]
  if (length(v) == 0) return(default)
  min(v)
}

safe_spearman <- function(x, y, min_n = 10) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  if (length(x2) < min_n || length(unique(x2)) < 3 || length(unique(y2)) < 3) {
    return(list(n = length(x2), rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x2, y2, method = "spearman", exact = FALSE))
  list(n = length(x2), rho = unname(ct$estimate), p = ct$p.value)
}

depth_residuals <- function(v, depth, df = 4) {
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v) & is.finite(depth)
  if (sum(ok) < 10 || length(unique(depth[ok])) < 6) return(out)
  fit <- tryCatch(lm(v[ok] ~ splines::ns(depth[ok], df = df)), error = function(e) NULL)
  if (is.null(fit)) return(out)
  out[ok] <- residuals(fit)
  out
}

layer_intervals <- function(depths) {
  d <- sort(unique(depths[is.finite(depths)]))
  mids <- (d[-1] + d[-length(d)]) / 2
  data.frame(
    depth_in_core_cm = d,
    lower = c(-Inf, mids),
    upper = c(mids, Inf)
  )
}

aggregate_xrf_to_layers <- function(xrf_df, intervals_df, numeric_cols, fun = c("median", "mean")) {
  fun <- match.arg(fun)
  out <- data.frame(depth_in_core_cm = intervals_df$depth_in_core_cm)
  out$n_xrf_rows <- NA_integer_
  for (nm in numeric_cols) out[[nm]] <- NA_real_

  for (i in seq_len(nrow(intervals_df))) {
    lo <- intervals_df$lower[i]
    hi <- intervals_df$upper[i]
    idx <- which(is.finite(xrf_df$depth_in_core_cm) & xrf_df$depth_in_core_cm > lo & xrf_df$depth_in_core_cm <= hi)
    out$n_xrf_rows[i] <- length(idx)
    if (length(idx) == 0) next
    for (nm in numeric_cols) {
      vals <- xrf_df[[nm]][idx]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) next
      out[[nm]][i] <- if (fun == "median") stats::median(vals) else mean(vals)
    }
  }
  out
}

add_ratio <- function(df, num_col, den_col, out_col) {
  if (!(num_col %in% names(df)) || !(den_col %in% names(df))) {
    df[[out_col]] <- NA_real_
    return(df)
  }
  den <- df[[den_col]]
  eps <- safe_min_positive(den, default = 1e-6) / 2
  den_adj <- ifelse(is.na(den), NA_real_, ifelse(den > 0, den, eps))
  df[[out_col]] <- df[[num_col]] / den_adj
  df
}

nearest_sparse_match <- function(layers_df, sparse_df, tol_cm = 3, source_label = "sparse") {
  out <- data.frame(depth_in_core_cm = layers_df$depth_in_core_cm, stringsAsFactors = FALSE)
  if (nrow(sparse_df) == 0) return(out)

  sparse_depth <- as.numeric(sparse_df$depth_in_core_cm)
  layer_depth <- as.numeric(layers_df$depth_in_core_cm)

  nearest_idx <- rep(NA_integer_, length(layer_depth))
  depth_diff <- rep(NA_real_, length(layer_depth))
  for (i in seq_along(layer_depth)) {
    if (!is.finite(layer_depth[i])) next
    j <- which.min(abs(sparse_depth - layer_depth[i]))
    nearest_idx[i] <- j
    depth_diff[i] <- abs(sparse_depth[j] - layer_depth[i])
  }

  within_tol <- is.finite(depth_diff) & depth_diff <= tol_cm
  num_cols <- names(sparse_df)[vapply(sparse_df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, c("depth_in_core_cm", "y_bp"))
  for (nm in num_cols) {
    vals <- rep(NA_real_, length(layer_depth))
    ok <- which(within_tol & !is.na(nearest_idx))
    if (length(ok) > 0) vals[ok] <- sparse_df[[nm]][nearest_idx[ok]]
    out[[paste0(source_label, "_", nm)]] <- vals
  }
  out[[paste0(source_label, "_matched_within_tolerance")]] <- within_tol
  out[[paste0(source_label, "_depth_diff_cm")]] <- depth_diff
  out
}

fit_depth_model <- function(df, response_col, trend_col, predictors, use_gls_if_available = TRUE) {
  use_cols <- unique(c(response_col, trend_col, predictors))
  use_cols <- use_cols[use_cols %in% names(df)]
  dat <- df[, use_cols, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < params$model_min_n) {
    return(list(model = NULL, model_type = "none", data = dat, coefficients = NULL, error = "too few complete rows"))
  }

  z_names <- character(0)
  for (p in predictors) {
    if (!(p %in% names(dat))) next
    s <- stats::sd(dat[[p]], na.rm = TRUE)
    if (!is.finite(s) || s == 0) next
    zn <- paste0("z_", p)
    dat[[zn]] <- (dat[[p]] - mean(dat[[p]], na.rm = TRUE)) / s
    z_names <- c(z_names, zn)
  }
  if (length(z_names) == 0) {
    return(list(model = NULL, model_type = "none", data = dat, coefficients = NULL, error = "no variable with non-zero variance"))
  }

  rhs <- paste(c(sprintf("splines::ns(%s, df = %d)", trend_col, params$depth_spline_df), z_names), collapse = " + ")
  fml <- as.formula(paste(response_col, "~", rhs))

  model <- NULL
  model_type <- "lm"
  if (use_gls_if_available && requireNamespace("nlme", quietly = TRUE)) {
    model_try <- tryCatch(
      nlme::gls(
        model = fml,
        data = dat,
        method = "REML",
        correlation = nlme::corCAR1(form = stats::as.formula(paste("~", trend_col)))
      ),
      error = function(e) NULL
    )
    if (!is.null(model_try)) {
      model <- model_try
      model_type <- "gls_corCAR1"
    }
  }
  if (is.null(model)) {
    model <- stats::lm(fml, data = dat)
    model_type <- "lm"
  }

  if (inherits(model, "gls")) {
    tt <- summary(model)$tTable
    coef_df <- data.frame(
      term = rownames(tt),
      estimate = tt[, "Value"],
      std_error = tt[, "Std.Error"],
      stat = tt[, "t-value"],
      p_value = tt[, "p-value"],
      stringsAsFactors = FALSE
    )
    ci <- tryCatch(nlme::intervals(model, which = "coef")$coef, error = function(e) NULL)
    coef_df$ci_low <- if (!is.null(ci)) ci[coef_df$term, "lower"] else NA_real_
    coef_df$ci_high <- if (!is.null(ci)) ci[coef_df$term, "upper"] else NA_real_
  } else {
    tt <- summary(model)$coefficients
    coef_df <- data.frame(
      term = rownames(tt),
      estimate = tt[, "Estimate"],
      std_error = tt[, "Std. Error"],
      stat = tt[, "t value"],
      p_value = tt[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    ci <- tryCatch(confint(model), error = function(e) NULL)
    coef_df$ci_low <- if (!is.null(ci)) ci[coef_df$term, 1] else NA_real_
    coef_df$ci_high <- if (!is.null(ci)) ci[coef_df$term, 2] else NA_real_
  }

  list(model = model, model_type = model_type, data = dat, coefficients = coef_df)
}

get_block_bootstrap_coefficients <- function(lm_model, source_data, n_boot = 600, block_len = NULL, seed = 13) {
  if (is.null(lm_model) || !inherits(lm_model, "lm") || is.null(source_data) || nrow(source_data) < params$model_min_n) return(NULL)
  n <- nrow(source_data)
  if (is.null(block_len)) block_len <- max(3, floor(sqrt(n)))
  n_blocks <- ceiling(n / block_len)
  coef_names <- names(stats::coef(lm_model))
  if (length(coef_names) == 0) return(NULL)

  set.seed(seed)
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(coef_names), dimnames = list(NULL, coef_names))
  fml <- stats::formula(lm_model)

  for (b in seq_len(n_boot)) {
    starts <- sample.int(n, size = n_blocks, replace = TRUE)
    idx <- integer(0)
    for (s in starts) idx <- c(idx, ((s - 1 + seq_len(block_len) - 1) %% n) + 1)
    idx <- idx[seq_len(n)]
    fit_b <- tryCatch(stats::lm(fml, data = source_data[idx, , drop = FALSE]), error = function(e) NULL)
    if (is.null(fit_b)) next
    cb <- stats::coef(fit_b)
    boot_mat[b, names(cb)] <- cb
  }

  est <- stats::coef(lm_model)
  out <- data.frame(
    term = coef_names,
    estimate = as.numeric(est[coef_names]),
    std_error = NA_real_,
    stat = NA_real_,
    p_value = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(coef_names)) {
    vv <- boot_mat[, coef_names[i]]
    vv <- vv[is.finite(vv)]
    if (length(vv) < max(100, floor(0.4 * n_boot))) next
    se <- stats::sd(vv)
    ci <- as.numeric(stats::quantile(vv, probs = c(0.025, 0.975), na.rm = TRUE, type = 7))
    pboot <- 2 * min(mean(vv <= 0, na.rm = TRUE), mean(vv >= 0, na.rm = TRUE))
    out$std_error[i] <- se
    out$stat[i] <- ifelse(is.finite(se) && se > 0, out$estimate[i] / se, NA_real_)
    out$p_value[i] <- min(1, pboot)
    out$ci_low[i] <- ci[1]
    out$ci_high[i] <- ci[2]
  }
  attr(out, "block_len") <- block_len
  out
}

write_tsv <- function(df, file) {
  utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

safe_spline_df <- function(x, target_df = 4) {
  ux <- length(unique(x[is.finite(x)]))
  if (!is.finite(ux) || ux < 2) return(NA_integer_)
  as.integer(max(1, min(target_df, ux - 1)))
}

safe_corr_sq <- function(obs, fit) {
  ok <- is.finite(obs) & is.finite(fit)
  if (sum(ok) < 3) return(NA_real_)
  r <- suppressWarnings(stats::cor(obs[ok], fit[ok]))
  if (!is.finite(r)) return(NA_real_)
  as.numeric(r^2)
}

pick_trend_col <- function(df) {
  if ("y_bp" %in% names(df)) {
    ok <- is.finite(df$y_bp)
    if (sum(ok) >= params$model_min_n && length(unique(df$y_bp[ok])) >= 10) return("y_bp")
  }
  if ("depth_in_core_cm" %in% names(df)) {
    ok <- is.finite(df$depth_in_core_cm)
    if (sum(ok) >= params$model_min_n && length(unique(df$depth_in_core_cm[ok])) >= 10) return("depth_in_core_cm")
  }
  NA_character_
}

fit_temporal_vs_plus_rb <- function(dat, trend_col, spline_df = 4) {
  if (nrow(dat) < params$model_min_n) return(NULL)
  df_use <- safe_spline_df(dat[[trend_col]], target_df = spline_df)
  if (!is.finite(df_use)) return(NULL)
  f0 <- as.formula(paste("library_concentration_log10 ~ splines::ns(", trend_col, ", df = ", df_use, ")", sep = ""))
  f1 <- as.formula(paste("library_concentration_log10 ~ splines::ns(", trend_col, ", df = ", df_use, ") + z_rb_sr", sep = ""))
  m0 <- tryCatch(stats::lm(f0, data = dat), error = function(e) NULL)
  m1 <- tryCatch(stats::lm(f1, data = dat), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(NULL)
  a0 <- tryCatch(AIC(m0), error = function(e) NA_real_)
  a1 <- tryCatch(AIC(m1), error = function(e) NA_real_)
  r20 <- safe_corr_sq(dat$library_concentration_log10, stats::fitted(m0))
  r21 <- safe_corr_sq(dat$library_concentration_log10, stats::fitted(m1))
  data.frame(
    n_obs = nrow(dat),
    spline_df = df_use,
    temporal_only_aic = a0,
    temporal_plus_rbsr_aic = a1,
    delta_aic_temporal_minus_plus = a0 - a1,
    temporal_only_corr_sq = r20,
    temporal_plus_rbsr_corr_sq = r21,
    delta_corr_sq_temporal_minus_plus = r21 - r20,
    stringsAsFactors = FALSE
  )
}

fit_rbsr_age_interaction <- function(dat, trend_col, spline_df = 4) {
  if (nrow(dat) < params$interaction_min_n) return(NULL)
  df_use <- safe_spline_df(dat[[trend_col]], target_df = spline_df)
  if (!is.finite(df_use)) return(NULL)
  f_int <- as.formula(paste(
    "library_concentration_log10 ~ splines::ns(", trend_col, ", df = ", df_use, ") + z_rb_sr + z_rb_sr:z_trend",
    sep = ""
  ))
  m_int <- tryCatch(stats::lm(f_int, data = dat), error = function(e) NULL)
  if (is.null(m_int)) return(NULL)
  tt <- summary(m_int)$coefficients
  b_rb <- if ("z_rb_sr" %in% rownames(tt)) tt["z_rb_sr", "Estimate"] else NA_real_
  p_rb <- if ("z_rb_sr" %in% rownames(tt)) tt["z_rb_sr", "Pr(>|t|)"] else NA_real_
  b_int <- if ("z_rb_sr:z_trend" %in% rownames(tt)) tt["z_rb_sr:z_trend", "Estimate"] else NA_real_
  p_int <- if ("z_rb_sr:z_trend" %in% rownames(tt)) tt["z_rb_sr:z_trend", "Pr(>|t|)"] else NA_real_
  zq <- stats::quantile(dat$z_trend, probs = c(0.25, 0.75), na.rm = TRUE, type = 7)
  slope_young <- if (is.finite(b_rb) && is.finite(b_int)) b_rb + b_int * zq[[1]] else NA_real_
  slope_old <- if (is.finite(b_rb) && is.finite(b_int)) b_rb + b_int * zq[[2]] else NA_real_
  data.frame(
    n_obs = nrow(dat),
    spline_df = df_use,
    rbsr_main_beta = b_rb,
    rbsr_main_p = p_rb,
    rbsr_age_interaction_beta = b_int,
    rbsr_age_interaction_p = p_int,
    implied_rbsr_slope_p25_age = slope_young,
    implied_rbsr_slope_p75_age = slope_old,
    stringsAsFactors = FALSE
  )
}

make_age_strata <- function(x, scheme = c("median_half", "tertile")) {
  scheme <- match.arg(scheme)
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (sum(ok) == 0) return(out)
  if (scheme == "median_half") {
    med <- stats::median(x[ok], na.rm = TRUE)
    out[ok & x <= med] <- "younger_half"
    out[ok & x > med] <- "older_half"
  } else {
    q <- stats::quantile(x[ok], probs = c(1/3, 2/3), na.rm = TRUE, type = 7)
    out[ok & x <= q[[1]]] <- "young_tertile"
    out[ok & x > q[[1]] & x <= q[[2]]] <- "mid_tertile"
    out[ok & x > q[[2]]] <- "old_tertile"
  }
  out
}

safe_linear_residuals <- function(v, trend) {
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v) & is.finite(trend)
  if (sum(ok) < 8 || length(unique(trend[ok])) < 4) return(out)
  fit <- tryCatch(stats::lm(v[ok] ~ trend[ok]), error = function(e) NULL)
  if (is.null(fit)) return(out)
  out[ok] <- stats::residuals(fit)
  out
}

pick_best_from_group <- function(group_vars, screen_df, min_n = 20) {
  sub <- screen_df[screen_df$variable %in% group_vars & is.finite(screen_df$rho_detrended) & screen_df$n_detrended >= min_n, , drop = FALSE]
  if (nrow(sub) == 0) {
    sub2 <- screen_df[screen_df$variable %in% group_vars & is.finite(screen_df$rho_raw) & screen_df$n_raw >= min_n, , drop = FALSE]
    if (nrow(sub2) == 0) return(NA_character_)
    return(sub2$variable[which.max(abs(sub2$rho_raw))[1]])
  }
  sub <- sub[order(-abs(sub$rho_detrended), sub$p_detrended), , drop = FALSE]
  sub$variable[1]
}

collapse_numeric <- function(v) {
  vv <- as.numeric(v)
  vv <- vv[is.finite(vv)]
  if (length(vv) == 0) return(NA_real_)
  mean(vv)
}

build_layer_table_st13_like <- function(core_label, metadata, xrf, foram, sst) {
  meta <- metadata[metadata$core == core_label, , drop = FALSE]
  if (nrow(meta) == 0) return(NULL)

  layer_key <- interaction(meta$depth_in_core_cm, meta$y_bp, drop = TRUE)
  idx_split <- split(seq_len(nrow(meta)), layer_key)

  layers <- data.frame(
    core = core_label,
    layer_id = sprintf("%s_layer_%03d", core_label, seq_along(idx_split)),
    depth_in_core_cm = vapply(idx_split, function(ix) collapse_numeric(meta$depth_in_core_cm[ix]), numeric(1)),
    y_bp = vapply(idx_split, function(ix) collapse_numeric(meta$y_bp[ix]), numeric(1)),
    library_concentration = vapply(idx_split, function(ix) collapse_numeric(meta$library_concentration[ix]), numeric(1)),
    temp = vapply(idx_split, function(ix) collapse_numeric(meta$temp[ix]), numeric(1)),
    mis = vapply(idx_split, function(ix) collapse_numeric(meta$mis[ix]), numeric(1)),
    initial = vapply(idx_split, function(ix) collapse_numeric(meta$initial[ix]), numeric(1)),
    derep = vapply(idx_split, function(ix) collapse_numeric(meta$derep[ix]), numeric(1)),
    avg_leng_initial = vapply(idx_split, function(ix) collapse_numeric(meta$avg_leng_initial[ix]), numeric(1)),
    avg_len_derep = vapply(idx_split, function(ix) collapse_numeric(meta$avg_len_derep[ix]), numeric(1)),
    n_rows_collapsed = vapply(idx_split, length, integer(1)),
    stringsAsFactors = FALSE
  )
  layers <- layers[order(layers$depth_in_core_cm), ]

  pseudocount <- safe_min_positive(layers$library_concentration, default = 1e-6) / 2
  layers$library_concentration_log10 <- log10(layers$library_concentration + pseudocount)

  xrf_core <- xrf[xrf$core == core_label, , drop = FALSE]
  if (nrow(xrf_core) == 0) return(list(layer_tbl = layers, pseudocount = pseudocount))

  numeric_xrf <- names(xrf_core)[vapply(xrf_core, is.numeric, logical(1))]
  numeric_xrf <- setdiff(numeric_xrf, c("depth_in_core_cm", "y_bp"))

  ints <- layer_intervals(layers$depth_in_core_cm)
  xrf_med <- aggregate_xrf_to_layers(xrf_core, ints, numeric_xrf, fun = "median")

  ratio_defs <- list(
    c("ba", "ti", "ba_ti"),
    c("p", "ti", "p_ti"),
    c("br", "ti", "br_ti"),
    c("si", "ti", "si_ti"),
    c("ca", "ti", "ca_ti"),
    c("al", "si", "al_si"),
    c("ti", "al", "ti_al"),
    c("fe", "al", "fe_al"),
    c("rb", "sr", "rb_sr"),
    c("mn", "fe", "mn_fe")
  )
  for (r in ratio_defs) xrf_med <- add_ratio(xrf_med, r[1], r[2], r[3])

  layer_tbl <- merge(layers, xrf_med, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
  layer_tbl <- layer_tbl[order(layer_tbl$depth_in_core_cm), ]

  foram_core <- foram[foram$core == core_label, , drop = FALSE]
  sst_core <- sst[sst$core == core_label, , drop = FALSE]
  if (nrow(foram_core) > 0) {
    fm <- nearest_sparse_match(layer_tbl, foram_core, tol_cm = params$sparse_match_tolerance_cm, source_label = "foram")
    layer_tbl <- merge(layer_tbl, fm, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
  }
  if (nrow(sst_core) > 0) {
    sm <- nearest_sparse_match(layer_tbl, sst_core, tol_cm = params$sparse_match_tolerance_cm, source_label = "sst")
    layer_tbl <- merge(layer_tbl, sm, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
  }
  layer_tbl <- layer_tbl[order(layer_tbl$depth_in_core_cm), ]

  list(layer_tbl = layer_tbl, pseudocount = pseudocount)
}

metadata <- read.delim(metadata_path, sep = "\t", check.names = FALSE)
xrf <- read.csv(xrf_path, check.names = FALSE)
foram <- read.delim(foram_path, sep = "\t", check.names = FALSE)
sst <- read.delim(sst_path, sep = "\t", check.names = FALSE)

cross_summary <- data.frame()
rb_effect_table <- data.frame()
per_core_interaction <- data.frame()
per_core_stratified <- data.frame()
temporal_vs_plus <- data.frame()

downstream_excluded <- c("initial", "derep", "avg_leng_initial", "avg_len_derep")

for (core_label in params$st13_like_cores) {
  message("Running ST13-like analysis for ", core_label, " ...")
  built <- build_layer_table_st13_like(core_label, metadata, xrf, foram, sst)
  if (is.null(built)) next
  layer_tbl <- built$layer_tbl

  write_tsv(layer_tbl, file.path(out_dir, sprintf("%s_layer_level_analysis_table.tsv", tolower(core_label))))

  metadata_env <- intersect(c("temp", "mis"), names(layer_tbl))
  xrf_ratio_candidates <- intersect(c("rb_sr", "ba_ti", "p_ti", "br_ti", "si_ti", "ca_ti", "al_si", "ti_al", "fe_al", "mn_fe"), names(layer_tbl))

  numeric_cols <- names(layer_tbl)[vapply(layer_tbl, is.numeric, logical(1))]
  xrf_raw_candidates <- intersect(numeric_cols, names(layer_tbl))
  xrf_raw_candidates <- setdiff(xrf_raw_candidates, c(
    "depth_in_core_cm", "y_bp", "library_concentration", "library_concentration_log10",
    downstream_excluded, "n_rows_collapsed", "n_xrf_rows"
  ))

  foram_candidates <- grep("^foram_", names(layer_tbl), value = TRUE)
  foram_candidates <- setdiff(foram_candidates, c("foram_matched_within_tolerance", "foram_depth_diff_cm"))
  sst_candidates <- grep("^sst_", names(layer_tbl), value = TRUE)
  sst_candidates <- setdiff(sst_candidates, c("sst_matched_within_tolerance", "sst_depth_diff_cm"))

  family_map <- c(
    setNames(rep("metadata_env", length(metadata_env)), metadata_env),
    setNames(rep("xrf_raw", length(xrf_raw_candidates)), xrf_raw_candidates),
    setNames(rep("xrf_ratio", length(xrf_ratio_candidates)), xrf_ratio_candidates),
    setNames(rep("foram_sparse", length(foram_candidates)), foram_candidates),
    setNames(rep("sst_sparse", length(sst_candidates)), sst_candidates)
  )
  all_candidates <- unique(names(family_map))
  all_candidates <- setdiff(all_candidates, downstream_excluded)

  response <- layer_tbl$library_concentration_log10
  depth <- layer_tbl$depth_in_core_cm
  res_response <- depth_residuals(response, depth, df = params$depth_spline_df)

  screen <- data.frame(
    variable = all_candidates,
    family = unname(family_map[all_candidates]),
    n_raw = NA_integer_,
    rho_raw = NA_real_,
    p_raw = NA_real_,
    n_detrended = NA_integer_,
    rho_detrended = NA_real_,
    p_detrended = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(screen))) {
    v <- screen$variable[i]
    x <- as.numeric(layer_tbl[[v]])
    raw <- safe_spearman(x, response, min_n = params$screening_min_n)
    screen$n_raw[i] <- raw$n
    screen$rho_raw[i] <- raw$rho
    screen$p_raw[i] <- raw$p

    x_res <- depth_residuals(x, depth, df = params$depth_spline_df)
    det <- safe_spearman(x_res, res_response, min_n = params$screening_min_n)
    screen$n_detrended[i] <- det$n
    screen$rho_detrended[i] <- det$rho
    screen$p_detrended[i] <- det$p
  }

  screen <- screen[order(-abs(screen$rho_detrended), screen$p_detrended), ]
  write_tsv(screen, file.path(out_dir, sprintf("%s_screening_correlations_raw_and_detrended.tsv", tolower(core_label))))

  prod_var <- pick_best_from_group(intersect(c("ba_ti", "p_ti", "br_ti", "ba", "p", "br"), names(layer_tbl)), screen, min_n = 20)
  terrig_var <- pick_best_from_group(intersect(c("rb_sr", "al_si", "ti_al", "fe_al", "si_ti", "rb", "si", "zr"), names(layer_tbl)), screen, min_n = 20)
  ca_var <- pick_best_from_group(intersect(c("ca_ti", "ca"), names(layer_tbl)), screen, min_n = 20)
  redox_var <- pick_best_from_group(intersect(c("mn_fe", "mn", "fe"), names(layer_tbl)), screen, min_n = 20)

  selected_predictors <- unique(na.omit(c(prod_var, terrig_var, ca_var, redox_var)))
  selection_tbl <- data.frame(
    core = core_label,
    process_group = c("productivity_or_organic", "terrigenous_or_mineralogical", "carbonate", "redox"),
    selected_predictor = c(prod_var, terrig_var, ca_var, redox_var),
    stringsAsFactors = FALSE
  )
  rb_screen_row <- screen[screen$variable == "rb_sr", c("n_raw", "rho_raw", "p_raw", "n_detrended", "rho_detrended", "p_detrended"), drop = FALSE]
  if (nrow(rb_screen_row) == 0) {
    rb_screen_row <- data.frame(n_raw = NA, rho_raw = NA, p_raw = NA, n_detrended = NA, rho_detrended = NA, p_detrended = NA)
  }
  rb_selection <- data.frame(
    core = core_label,
    process_group = "rb_sr_tracking",
    selected_predictor = ifelse("rb_sr" %in% selected_predictors, "rb_sr", "not_selected"),
    rb_sr_n_raw = rb_screen_row$n_raw[1],
    rb_sr_rho_raw = rb_screen_row$rho_raw[1],
    rb_sr_p_raw = rb_screen_row$p_raw[1],
    rb_sr_n_detrended = rb_screen_row$n_detrended[1],
    rb_sr_rho_detrended = rb_screen_row$rho_detrended[1],
    rb_sr_p_detrended = rb_screen_row$p_detrended[1],
    stringsAsFactors = FALSE
  )
  write_tsv(selection_tbl, file.path(out_dir, sprintf("%s_primary_model_selected_predictors.tsv", tolower(core_label))))
  write_tsv(rb_selection, file.path(out_dir, sprintf("%s_rb_sr_selection_tracking.tsv", tolower(core_label))))

  primary_fit <- fit_depth_model(
    df = layer_tbl,
    response_col = "library_concentration_log10",
    trend_col = "depth_in_core_cm",
    predictors = selected_predictors,
    use_gls_if_available = TRUE
  )

  if (!is.null(primary_fit$coefficients)) {
    primary_coef <- primary_fit$coefficients
    primary_coef$core <- core_label
    primary_coef$model_type <- primary_fit$model_type
    write_tsv(primary_coef, file.path(out_dir, sprintf("%s_primary_model_summary_table.tsv", tolower(core_label))))
  } else {
    write_tsv(data.frame(core = core_label, status = "not_fit", reason = primary_fit$error),
              file.path(out_dir, sprintf("%s_primary_model_summary_table.tsv", tolower(core_label))))
  }

  # Explicit rb_sr model for cross-core comparability
  rb_predictors <- c("rb_sr")
  if (!("rb_sr" %in% names(layer_tbl))) {
    rb_predictors <- character(0)
  }

  rb_fit <- NULL
  rb_fit_robust <- NULL
  if (length(rb_predictors) > 0) {
    rb_fit <- fit_depth_model(
      df = layer_tbl,
      response_col = "library_concentration_log10",
      trend_col = "depth_in_core_cm",
      predictors = rb_predictors,
      use_gls_if_available = TRUE
    )

    if (!is.null(rb_fit$model) && rb_fit$model_type == "lm") {
      bb <- get_block_bootstrap_coefficients(
        lm_model = rb_fit$model,
        source_data = rb_fit$data,
        n_boot = params$block_bootstrap_n
      )
      if (!is.null(bb)) {
        bb$core <- core_label
        bb$model_type <- "lm_block_bootstrap"
        rb_fit_robust <- bb
      }
    }
  }

  if (!is.null(rb_fit) && !is.null(rb_fit$coefficients)) {
    rb_coef <- rb_fit$coefficients
    rb_coef$core <- core_label
    rb_coef$model_type <- rb_fit$model_type
    write_tsv(rb_coef, file.path(out_dir, sprintf("%s_rb_sr_depth_adjusted_model_summary.tsv", tolower(core_label))))
    if (!is.null(rb_fit_robust)) {
      write_tsv(rb_fit_robust, file.path(out_dir, sprintf("%s_rb_sr_depth_adjusted_model_block_bootstrap_summary.tsv", tolower(core_label))))
    }
  } else {
    write_tsv(data.frame(core = core_label, status = "not_fit", reason = if (!is.null(rb_fit)) rb_fit$error else "rb_sr missing"),
              file.path(out_dir, sprintf("%s_rb_sr_depth_adjusted_model_summary.tsv", tolower(core_label))))
  }

  # Sparse sensitivity (foram + sst) if feasible
  sst_pref <- intersect(c("sst_sst_uk37_alkenone", "sst_sst_mgca_jonkers_2013", "sst_sst_mgca_kozdon_2009"), names(layer_tbl))
  foram_pref <- intersect(c("foram_g_bulloides_pct", "foram_polar_planktonic_spp_pct", "foram_n_pachyderma_pct"), names(layer_tbl))
  sens_sparse <- NULL
  if (length(sst_pref) > 0 && length(foram_pref) > 0) {
    sens_sparse <- fit_depth_model(
      df = layer_tbl,
      response_col = "library_concentration_log10",
      trend_col = "depth_in_core_cm",
      predictors = unique(c(selected_predictors, sst_pref[1], foram_pref[1])),
      use_gls_if_available = FALSE
    )
  }
  sens_overview <- data.frame(
    core = core_label,
    sensitivity_model = "add_sparse_sst_and_foram",
    status = ifelse(is.null(sens_sparse) || is.null(sens_sparse$model), "not_fit", "fit"),
    n_obs = ifelse(is.null(sens_sparse) || is.null(sens_sparse$model), NA, nrow(sens_sparse$data)),
    note = ifelse(
      is.null(sens_sparse) || is.null(sens_sparse$model),
      "sparse predictors unavailable or too few complete rows",
      paste0("added ", sst_pref[1], " + ", foram_pref[1])
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(sens_overview, file.path(out_dir, sprintf("%s_sensitivity_model_overview.tsv", tolower(core_label))))

  rb_raw <- screen[screen$variable == "rb_sr", , drop = FALSE]
  if (nrow(rb_raw) == 0) rb_raw <- data.frame(rho_raw = NA_real_, p_raw = NA_real_, rho_detrended = NA_real_, p_detrended = NA_real_, n_raw = NA_real_)

  rb_term <- data.frame(estimate = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n_obs = NA_real_, model_type = NA_character_)
  if (!is.null(rb_fit) && !is.null(rb_fit$coefficients)) {
    rr <- rb_fit$coefficients[rb_fit$coefficients$term == "z_rb_sr", , drop = FALSE]
    if (nrow(rr) == 1) {
      rb_term$estimate <- rr$estimate
      rb_term$p_value <- rr$p_value
      rb_term$ci_low <- rr$ci_low
      rb_term$ci_high <- rr$ci_high
      rb_term$n_obs <- nrow(rb_fit$data)
      rb_term$model_type <- rb_fit$model_type
    }
  }

  rb_effect_table <- rbind(rb_effect_table, data.frame(
    core = core_label,
    analysis_level = "st13_like_depth_adjusted_rb_sr",
    n_obs = rb_term$n_obs,
    rb_sr_raw_rho = rb_raw$rho_raw[1],
    rb_sr_raw_p = rb_raw$p_raw[1],
    rb_sr_detrended_rho = rb_raw$rho_detrended[1],
    rb_sr_detrended_p = rb_raw$p_detrended[1],
    rb_sr_model_estimate = rb_term$estimate,
    rb_sr_model_p = rb_term$p_value,
    rb_sr_model_ci_low = rb_term$ci_low,
    rb_sr_model_ci_high = rb_term$ci_high,
    rb_sr_model_type = rb_term$model_type,
    stringsAsFactors = FALSE
  ))

  supported <- ifelse(
    is.finite(rb_term$estimate) & is.finite(rb_term$p_value),
    rb_term$estimate < 0 & rb_term$p_value <= params$support_p_threshold,
    is.finite(rb_raw$rho_detrended[1]) & rb_raw$rho_detrended[1] < 0 & is.finite(rb_raw$p_detrended[1]) & rb_raw$p_detrended[1] <= params$support_p_threshold
  )

  cross_summary <- rbind(cross_summary, data.frame(
    core = core_label,
    analysis_level = "st13_like_multivariable",
    n_obs = ifelse(is.finite(rb_term$n_obs), rb_term$n_obs, rb_raw$n_raw[1]),
    rb_sr_raw_rho = rb_raw$rho_raw[1],
    rb_sr_raw_p = rb_raw$p_raw[1],
    rb_sr_detrended_rho = rb_raw$rho_detrended[1],
    model_estimate = rb_term$estimate,
    rb_sr_direction = ifelse(is.finite(rb_term$estimate), ifelse(rb_term$estimate < 0, "negative", ifelse(rb_term$estimate > 0, "positive", "zero")), ifelse(rb_raw$rho_detrended[1] < 0, "negative", ifelse(rb_raw$rho_detrended[1] > 0, "positive", "zero"))),
    rb_sr_supported_yes_no = ifelse(isTRUE(supported), "yes", "no"),
    notes = "ST13-like depth-aware model; sparse foram/SST treated as sensitivity-only.",
    stringsAsFactors = FALSE
  ))

  # Age-dependent hypothesis tests: interaction + age-stratified + temporal-vs-plus-rb_sr
  trend_col <- pick_trend_col(layer_tbl)
  if (!is.na(trend_col) && all(c("library_concentration_log10", "rb_sr") %in% names(layer_tbl))) {
    dat_h <- layer_tbl[, c("library_concentration_log10", "rb_sr", trend_col), drop = FALSE]
    names(dat_h)[names(dat_h) == trend_col] <- "trend"
    dat_h <- dat_h[complete.cases(dat_h), , drop = FALSE]
    if (nrow(dat_h) >= params$interaction_min_n && stats::sd(dat_h$rb_sr) > 0 && stats::sd(dat_h$trend) > 0) {
      dat_h$z_rb_sr <- as.numeric(scale(dat_h$rb_sr))
      dat_h$z_trend <- as.numeric(scale(dat_h$trend))
      dat_h$y_bp <- if (trend_col == "y_bp") dat_h$trend else NA_real_
      dat_h$depth_in_core_cm <- if (trend_col == "depth_in_core_cm") dat_h$trend else NA_real_

      int_row <- fit_rbsr_age_interaction(
        dat = dat_h,
        trend_col = if (trend_col == "y_bp") "y_bp" else "depth_in_core_cm",
        spline_df = params$depth_spline_df
      )
      if (!is.null(int_row)) {
        int_row$core <- core_label
        int_row$trend_col <- trend_col
        int_row$status <- "fit"
      } else {
        int_row <- data.frame(
          core = core_label,
          trend_col = trend_col,
          n_obs = nrow(dat_h),
          spline_df = NA_real_,
          rbsr_main_beta = NA_real_,
          rbsr_main_p = NA_real_,
          rbsr_age_interaction_beta = NA_real_,
          rbsr_age_interaction_p = NA_real_,
          implied_rbsr_slope_p25_age = NA_real_,
          implied_rbsr_slope_p75_age = NA_real_,
          status = "not_fit",
          stringsAsFactors = FALSE
        )
      }
      per_core_interaction <- rbind(per_core_interaction, int_row)

      tvp <- fit_temporal_vs_plus_rb(
        dat = dat_h,
        trend_col = if (trend_col == "y_bp") "y_bp" else "depth_in_core_cm",
        spline_df = params$depth_spline_df
      )
      if (!is.null(tvp)) {
        tvp$core <- core_label
        tvp$scope <- "whole_core"
        tvp$scheme <- "none"
        tvp$stratum <- "all"
        tvp$trend_col <- trend_col
        tvp$status <- "fit"
        temporal_vs_plus <- rbind(temporal_vs_plus, tvp)
      }

      for (scheme in c("median_half", "tertile")) {
        strata <- make_age_strata(dat_h$trend, scheme = scheme)
        levs <- unique(stats::na.omit(strata))
        for (lv in levs) {
          sub <- dat_h[strata == lv, , drop = FALSE]
          nsub <- nrow(sub)
          raw <- safe_spearman(sub$rb_sr, sub$library_concentration_log10, min_n = params$stratum_min_n)
          y_res <- safe_linear_residuals(sub$library_concentration_log10, sub$trend)
          x_res <- safe_linear_residuals(sub$rb_sr, sub$trend)
          det <- safe_spearman(x_res, y_res, min_n = params$stratum_min_n)

          per_core_stratified <- rbind(per_core_stratified, data.frame(
            core = core_label,
            trend_col = trend_col,
            scheme = scheme,
            stratum = lv,
            n = nsub,
            age_min = min(sub$trend, na.rm = TRUE),
            age_max = max(sub$trend, na.rm = TRUE),
            rb_sr_rho_raw = raw$rho,
            rb_sr_p_raw = raw$p,
            rb_sr_rho_age_adjusted_linear = det$rho,
            rb_sr_p_age_adjusted_linear = det$p,
            stringsAsFactors = FALSE
          ))

          if (nsub >= params$model_min_n) {
            tvp_s <- fit_temporal_vs_plus_rb(
              dat = sub,
              trend_col = if (trend_col == "y_bp") "y_bp" else "depth_in_core_cm",
              spline_df = params$depth_spline_df
            )
            if (!is.null(tvp_s)) {
              tvp_s$core <- core_label
              tvp_s$scope <- "age_stratum"
              tvp_s$scheme <- scheme
              tvp_s$stratum <- lv
              tvp_s$trend_col <- trend_col
              tvp_s$status <- "fit"
              temporal_vs_plus <- rbind(temporal_vs_plus, tvp_s)
            }
          }
        }
      }
    } else {
      per_core_interaction <- rbind(per_core_interaction, data.frame(
        core = core_label,
        trend_col = trend_col,
        n_obs = nrow(dat_h),
        spline_df = NA_real_,
        rbsr_main_beta = NA_real_,
        rbsr_main_p = NA_real_,
        rbsr_age_interaction_beta = NA_real_,
        rbsr_age_interaction_p = NA_real_,
        implied_rbsr_slope_p25_age = NA_real_,
        implied_rbsr_slope_p75_age = NA_real_,
        status = "not_fit_low_n_or_variance",
        stringsAsFactors = FALSE
      ))
    }
  } else {
    per_core_interaction <- rbind(per_core_interaction, data.frame(
      core = core_label,
      trend_col = ifelse(is.na(trend_col), NA_character_, trend_col),
      n_obs = NA_real_,
      spline_df = NA_real_,
      rbsr_main_beta = NA_real_,
      rbsr_main_p = NA_real_,
      rbsr_age_interaction_beta = NA_real_,
      rbsr_age_interaction_p = NA_real_,
      implied_rbsr_slope_p25_age = NA_real_,
      implied_rbsr_slope_p75_age = NA_real_,
      status = "not_fit_missing_trend_or_rbsr",
      stringsAsFactors = FALSE
    ))
  }
}

# GeoB25202 basic rb_sr validation
message("Running basic rb_sr analysis for GeoB25202 ...")
meta_geob <- metadata[grepl(paste0("^", params$geob_base_core, "(_R[0-9]+)?$"), metadata$core), , drop = FALSE]
if (nrow(meta_geob) > 0) {
  meta_geob$core_clean <- sub(params$geob_replicate_regex, "", meta_geob$core)
  key <- interaction(meta_geob$core_clean, meta_geob$depth_in_core_cm, meta_geob$y_bp, drop = TRUE)
  idx_split <- split(seq_len(nrow(meta_geob)), key)

  geob_layers <- data.frame(
    core = params$geob_base_core,
    core_clean = vapply(idx_split, function(ix) unique(meta_geob$core_clean[ix])[1], character(1)),
    depth_in_core_cm = vapply(idx_split, function(ix) collapse_numeric(meta_geob$depth_in_core_cm[ix]), numeric(1)),
    y_bp = vapply(idx_split, function(ix) collapse_numeric(meta_geob$y_bp[ix]), numeric(1)),
    library_concentration = vapply(idx_split, function(ix) collapse_numeric(meta_geob$library_concentration[ix]), numeric(1)),
    n_technical_replicates = vapply(idx_split, length, integer(1)),
    stringsAsFactors = FALSE
  )
  geob_layers <- geob_layers[order(geob_layers$depth_in_core_cm), ]
  geob_pc <- safe_min_positive(geob_layers$library_concentration, default = 1e-6) / 2
  geob_layers$library_concentration_log10 <- log10(geob_layers$library_concentration + geob_pc)

  xrf_geob <- xrf[xrf$core == params$geob_base_core, , drop = FALSE]
  if (nrow(xrf_geob) > 0) {
    ints <- layer_intervals(geob_layers$depth_in_core_cm)
    xrf_num <- names(xrf_geob)[vapply(xrf_geob, is.numeric, logical(1))]
    xrf_num <- setdiff(xrf_num, c("depth_in_core_cm", "y_bp"))
    xrf_ag <- aggregate_xrf_to_layers(xrf_geob, ints, xrf_num, fun = "median")
    if (!("rb_sr" %in% names(xrf_ag))) xrf_ag <- add_ratio(xrf_ag, "rb", "sr", "rb_sr")
    geob_layers <- merge(geob_layers, xrf_ag[, c("depth_in_core_cm", "rb_sr", "n_xrf_rows"), drop = FALSE], by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
    geob_layers <- geob_layers[order(geob_layers$depth_in_core_cm), ]
  } else {
    geob_layers$rb_sr <- NA_real_
  }

  write_tsv(geob_layers, file.path(out_dir, "geob25202_layer_level_analysis_table.tsv"))

  raw <- safe_spearman(geob_layers$rb_sr, geob_layers$library_concentration_log10, min_n = params$screening_min_n)
  rb_res <- depth_residuals(geob_layers$rb_sr, geob_layers$depth_in_core_cm, df = params$depth_spline_df)
  y_res <- depth_residuals(geob_layers$library_concentration_log10, geob_layers$depth_in_core_cm, df = params$depth_spline_df)
  det <- safe_spearman(rb_res, y_res, min_n = params$screening_min_n)

  geob_model <- NULL
  geob_model_row <- data.frame(model_estimate = NA_real_, model_p = NA_real_, model_ci_low = NA_real_, model_ci_high = NA_real_, n_model = NA_real_)
  if (sum(is.finite(geob_layers$rb_sr) & is.finite(geob_layers$library_concentration_log10) & is.finite(geob_layers$depth_in_core_cm)) >= params$model_min_n) {
    dat <- geob_layers[, c("library_concentration_log10", "depth_in_core_cm", "rb_sr"), drop = FALSE]
    dat <- dat[complete.cases(dat), , drop = FALSE]
    dat$z_rb_sr <- as.numeric(scale(dat$rb_sr))
    if (all(is.finite(dat$z_rb_sr))) {
      geob_model <- lm(library_concentration_log10 ~ splines::ns(depth_in_core_cm, df = params$depth_spline_df) + z_rb_sr, data = dat)
      tt <- summary(geob_model)$coefficients
      if ("z_rb_sr" %in% rownames(tt)) {
        ci <- tryCatch(confint(geob_model), error = function(e) NULL)
        geob_model_row <- data.frame(
          model_estimate = tt["z_rb_sr", "Estimate"],
          model_p = tt["z_rb_sr", "Pr(>|t|)"],
          model_ci_low = if (!is.null(ci)) ci["z_rb_sr", 1] else NA_real_,
          model_ci_high = if (!is.null(ci)) ci["z_rb_sr", 2] else NA_real_,
          n_model = nrow(dat)
        )
      }
    }
  }

  geob_basic <- data.frame(
    core = params$geob_base_core,
    n_layers = nrow(geob_layers),
    n_raw = raw$n,
    rb_sr_raw_rho = raw$rho,
    rb_sr_raw_p = raw$p,
    n_detrended = det$n,
    rb_sr_detrended_rho = det$rho,
    rb_sr_detrended_p = det$p,
    rb_sr_depth_adjusted_model_estimate = geob_model_row$model_estimate,
    rb_sr_depth_adjusted_model_p = geob_model_row$model_p,
    rb_sr_depth_adjusted_model_ci_low = geob_model_row$model_ci_low,
    rb_sr_depth_adjusted_model_ci_high = geob_model_row$model_ci_high,
    rb_sr_depth_adjusted_model_n = geob_model_row$n_model,
    replicate_collapse_rule = "core_clean + depth_in_core_cm + y_bp mean collapse",
    stringsAsFactors = FALSE
  )
  write_tsv(geob_basic, file.path(out_dir, "geob25202_basic_rb_sr_validation.tsv"))

  rb_effect_table <- rbind(rb_effect_table, data.frame(
    core = params$geob_base_core,
    analysis_level = "basic_depth_adjusted_rb_sr",
    n_obs = geob_model_row$n_model,
    rb_sr_raw_rho = raw$rho,
    rb_sr_raw_p = raw$p,
    rb_sr_detrended_rho = det$rho,
    rb_sr_detrended_p = det$p,
    rb_sr_model_estimate = geob_model_row$model_estimate,
    rb_sr_model_p = geob_model_row$model_p,
    rb_sr_model_ci_low = geob_model_row$model_ci_low,
    rb_sr_model_ci_high = geob_model_row$model_ci_high,
    rb_sr_model_type = "lm_depth_spline",
    stringsAsFactors = FALSE
  ))

  geob_supported <- ifelse(
    is.finite(geob_model_row$model_estimate) & is.finite(geob_model_row$model_p),
    geob_model_row$model_estimate < 0 & geob_model_row$model_p <= params$support_p_threshold,
    is.finite(det$rho) & is.finite(det$p) & det$rho < 0 & det$p <= params$support_p_threshold
  )

  cross_summary <- rbind(cross_summary, data.frame(
    core = params$geob_base_core,
    analysis_level = "basic_rb_sr",
    n_obs = ifelse(is.finite(geob_model_row$n_model), geob_model_row$n_model, raw$n),
    rb_sr_raw_rho = raw$rho,
    rb_sr_raw_p = raw$p,
    rb_sr_detrended_rho = det$rho,
    model_estimate = geob_model_row$model_estimate,
    rb_sr_direction = ifelse(is.finite(geob_model_row$model_estimate), ifelse(geob_model_row$model_estimate < 0, "negative", ifelse(geob_model_row$model_estimate > 0, "positive", "zero")), ifelse(det$rho < 0, "negative", ifelse(det$rho > 0, "positive", "zero"))),
    rb_sr_supported_yes_no = ifelse(isTRUE(geob_supported), "yes", "no"),
    notes = "Basic GeoB25202-only rb_sr validation after conservative technical replicate collapse.",
    stringsAsFactors = FALSE
  ))

  # Age-dependent hypothesis tests for GeoB25202 basic table
  trend_col <- pick_trend_col(geob_layers)
  if (!is.na(trend_col) && all(c("library_concentration_log10", "rb_sr") %in% names(geob_layers))) {
    dat_h <- geob_layers[, c("library_concentration_log10", "rb_sr", trend_col), drop = FALSE]
    names(dat_h)[names(dat_h) == trend_col] <- "trend"
    dat_h <- dat_h[complete.cases(dat_h), , drop = FALSE]
    if (nrow(dat_h) >= params$interaction_min_n && stats::sd(dat_h$rb_sr) > 0 && stats::sd(dat_h$trend) > 0) {
      dat_h$z_rb_sr <- as.numeric(scale(dat_h$rb_sr))
      dat_h$z_trend <- as.numeric(scale(dat_h$trend))
      dat_h$y_bp <- if (trend_col == "y_bp") dat_h$trend else NA_real_
      dat_h$depth_in_core_cm <- if (trend_col == "depth_in_core_cm") dat_h$trend else NA_real_

      int_row <- fit_rbsr_age_interaction(
        dat = dat_h,
        trend_col = if (trend_col == "y_bp") "y_bp" else "depth_in_core_cm",
        spline_df = params$depth_spline_df
      )
      if (!is.null(int_row)) {
        int_row$core <- params$geob_base_core
        int_row$trend_col <- trend_col
        int_row$status <- "fit"
      } else {
        int_row <- data.frame(
          core = params$geob_base_core,
          trend_col = trend_col,
          n_obs = nrow(dat_h),
          spline_df = NA_real_,
          rbsr_main_beta = NA_real_,
          rbsr_main_p = NA_real_,
          rbsr_age_interaction_beta = NA_real_,
          rbsr_age_interaction_p = NA_real_,
          implied_rbsr_slope_p25_age = NA_real_,
          implied_rbsr_slope_p75_age = NA_real_,
          status = "not_fit",
          stringsAsFactors = FALSE
        )
      }
      per_core_interaction <- rbind(per_core_interaction, int_row)

      tvp <- fit_temporal_vs_plus_rb(
        dat = dat_h,
        trend_col = if (trend_col == "y_bp") "y_bp" else "depth_in_core_cm",
        spline_df = params$depth_spline_df
      )
      if (!is.null(tvp)) {
        tvp$core <- params$geob_base_core
        tvp$scope <- "whole_core"
        tvp$scheme <- "none"
        tvp$stratum <- "all"
        tvp$trend_col <- trend_col
        tvp$status <- "fit"
        temporal_vs_plus <- rbind(temporal_vs_plus, tvp)
      }

      for (scheme in c("median_half", "tertile")) {
        strata <- make_age_strata(dat_h$trend, scheme = scheme)
        levs <- unique(stats::na.omit(strata))
        for (lv in levs) {
          sub <- dat_h[strata == lv, , drop = FALSE]
          nsub <- nrow(sub)
          raw <- safe_spearman(sub$rb_sr, sub$library_concentration_log10, min_n = params$stratum_min_n)
          y_res <- safe_linear_residuals(sub$library_concentration_log10, sub$trend)
          x_res <- safe_linear_residuals(sub$rb_sr, sub$trend)
          det <- safe_spearman(x_res, y_res, min_n = params$stratum_min_n)

          per_core_stratified <- rbind(per_core_stratified, data.frame(
            core = params$geob_base_core,
            trend_col = trend_col,
            scheme = scheme,
            stratum = lv,
            n = nsub,
            age_min = min(sub$trend, na.rm = TRUE),
            age_max = max(sub$trend, na.rm = TRUE),
            rb_sr_rho_raw = raw$rho,
            rb_sr_p_raw = raw$p,
            rb_sr_rho_age_adjusted_linear = det$rho,
            rb_sr_p_age_adjusted_linear = det$p,
            stringsAsFactors = FALSE
          ))

          if (nsub >= params$model_min_n) {
            tvp_s <- fit_temporal_vs_plus_rb(
              dat = sub,
              trend_col = if (trend_col == "y_bp") "y_bp" else "depth_in_core_cm",
              spline_df = params$depth_spline_df
            )
            if (!is.null(tvp_s)) {
              tvp_s$core <- params$geob_base_core
              tvp_s$scope <- "age_stratum"
              tvp_s$scheme <- scheme
              tvp_s$stratum <- lv
              tvp_s$trend_col <- trend_col
              tvp_s$status <- "fit"
              temporal_vs_plus <- rbind(temporal_vs_plus, tvp_s)
            }
          }
        }
      }
    } else {
      per_core_interaction <- rbind(per_core_interaction, data.frame(
        core = params$geob_base_core,
        trend_col = trend_col,
        n_obs = nrow(dat_h),
        spline_df = NA_real_,
        rbsr_main_beta = NA_real_,
        rbsr_main_p = NA_real_,
        rbsr_age_interaction_beta = NA_real_,
        rbsr_age_interaction_p = NA_real_,
        implied_rbsr_slope_p25_age = NA_real_,
        implied_rbsr_slope_p75_age = NA_real_,
        status = "not_fit_low_n_or_variance",
        stringsAsFactors = FALSE
      ))
    }
  } else {
    per_core_interaction <- rbind(per_core_interaction, data.frame(
      core = params$geob_base_core,
      trend_col = ifelse(is.na(trend_col), NA_character_, trend_col),
      n_obs = NA_real_,
      spline_df = NA_real_,
      rbsr_main_beta = NA_real_,
      rbsr_main_p = NA_real_,
      rbsr_age_interaction_beta = NA_real_,
      rbsr_age_interaction_p = NA_real_,
      implied_rbsr_slope_p25_age = NA_real_,
      implied_rbsr_slope_p75_age = NA_real_,
      status = "not_fit_missing_trend_or_rbsr",
      stringsAsFactors = FALSE
    ))
  }
}

write_tsv(rb_effect_table, file.path(out_dir, "multicore_rbsr_effect_details.tsv"))
write_tsv(cross_summary, file.path(out_dir, "multicore_rbsr_summary.tsv"))

if (nrow(per_core_interaction) == 0) {
  per_core_interaction <- data.frame(
    core = character(), trend_col = character(), n_obs = numeric(), spline_df = numeric(),
    rbsr_main_beta = numeric(), rbsr_main_p = numeric(), rbsr_age_interaction_beta = numeric(),
    rbsr_age_interaction_p = numeric(), implied_rbsr_slope_p25_age = numeric(),
    implied_rbsr_slope_p75_age = numeric(), status = character(), stringsAsFactors = FALSE
  )
}
if (nrow(per_core_stratified) == 0) {
  per_core_stratified <- data.frame(
    core = character(), trend_col = character(), scheme = character(), stratum = character(),
    n = numeric(), age_min = numeric(), age_max = numeric(), rb_sr_rho_raw = numeric(),
    rb_sr_p_raw = numeric(), rb_sr_rho_age_adjusted_linear = numeric(),
    rb_sr_p_age_adjusted_linear = numeric(), stringsAsFactors = FALSE
  )
}
if (nrow(temporal_vs_plus) == 0) {
  temporal_vs_plus <- data.frame(
    n_obs = numeric(), spline_df = numeric(), temporal_only_aic = numeric(), temporal_plus_rbsr_aic = numeric(),
    delta_aic_temporal_minus_plus = numeric(), temporal_only_corr_sq = numeric(), temporal_plus_rbsr_corr_sq = numeric(),
    delta_corr_sq_temporal_minus_plus = numeric(), core = character(), scope = character(), scheme = character(),
    stratum = character(), trend_col = character(), status = character(), stringsAsFactors = FALSE
  )
}

write_tsv(per_core_interaction, file.path(out_dir, "per_core_rbsr_age_interaction_summary.tsv"))
write_tsv(per_core_stratified, file.path(out_dir, "per_core_age_stratified_rbsr_summary.tsv"))
write_tsv(temporal_vs_plus, file.path(out_dir, "per_core_temporal_vs_temporal_plus_rbsr_model_comparison.tsv"))

# Build compact hypothesis summary
older_younger <- per_core_stratified[
  per_core_stratified$scheme == "median_half" & per_core_stratified$stratum %in% c("younger_half", "older_half"),
  , drop = FALSE
]
dy <- data.frame()
if (nrow(older_younger) > 0) {
  for (cc in unique(older_younger$core)) {
    sub <- older_younger[older_younger$core == cc, , drop = FALSE]
    yo <- sub[sub$stratum == "younger_half", , drop = FALSE]
    od <- sub[sub$stratum == "older_half", , drop = FALSE]
    if (nrow(yo) == 1 && nrow(od) == 1) {
      dy <- rbind(dy, data.frame(
        core = cc,
        younger_rho = yo$rb_sr_rho_age_adjusted_linear,
        older_rho = od$rb_sr_rho_age_adjusted_linear,
        older_minus_younger = od$rb_sr_rho_age_adjusted_linear - yo$rb_sr_rho_age_adjusted_linear,
        younger_n = yo$n,
        older_n = od$n,
        stringsAsFactors = FALSE
      ))
    }
  }
}

core_interpret <- data.frame()
for (cc in unique(c(per_core_interaction$core, dy$core))) {
  ir <- per_core_interaction[per_core_interaction$core == cc & per_core_interaction$status == "fit", , drop = FALSE]
  dr <- dy[dy$core == cc, , drop = FALSE]
  support_flag <- NA_character_
  note <- ""
  if (nrow(ir) == 1 && is.finite(ir$rbsr_age_interaction_beta) && is.finite(ir$rbsr_age_interaction_p)) {
    if (ir$rbsr_age_interaction_beta < 0 && ir$rbsr_age_interaction_p <= 0.10) {
      support_flag <- "supports"
      note <- "interaction term indicates stronger negative rb_sr slope with older age"
    }
  }
  if (is.na(support_flag) && nrow(dr) == 1 && is.finite(dr$older_minus_younger)) {
    if (dr$older_minus_younger < 0) {
      support_flag <- "supports"
      note <- ifelse(note == "", "older half is more negative than younger half", paste(note, "; older half more negative"))
    } else {
      support_flag <- "not_supporting"
      note <- ifelse(note == "", "older half is not more negative than younger half", note)
    }
  }
  if (is.na(support_flag)) {
    support_flag <- "inconclusive"
    note <- ifelse(note == "", "insufficient information for age-dependent direction", note)
  }
  core_interpret <- rbind(core_interpret, data.frame(
    core = cc,
    support_flag = support_flag,
    note = note,
    stringsAsFactors = FALSE
  ))
}

overall_call <- "mixed"
if (nrow(core_interpret) > 0) {
  n_support <- sum(core_interpret$support_flag == "supports", na.rm = TRUE)
  n_not <- sum(core_interpret$support_flag == "not_supporting", na.rm = TRUE)
  if (n_support > 0 && n_not == 0) overall_call <- "supported"
  if (n_support == 0 && n_not > 0) overall_call <- "weakened"
}

hyp_summary <- data.frame(
  summary_scope = "multicore_age_dependent_rbsr_hypothesis",
  overall_call = overall_call,
  n_cores_with_interaction_fit = sum(per_core_interaction$status == "fit", na.rm = TRUE),
  n_cores_supporting_pattern = sum(core_interpret$support_flag == "supports", na.rm = TRUE),
  n_cores_not_supporting_pattern = sum(core_interpret$support_flag == "not_supporting", na.rm = TRUE),
  n_cores_inconclusive = sum(core_interpret$support_flag == "inconclusive", na.rm = TRUE),
  interpretation_guardrail = "Across-core differences confound age span/maximum age with core identity; pattern support is suggestive not causal proof.",
  stringsAsFactors = FALSE
)
write_tsv(hyp_summary, file.path(out_dir, "age_dependent_rbsr_hypothesis_summary.tsv"))
write_tsv(core_interpret, file.path(out_dir, "age_dependent_rbsr_hypothesis_per_core_interpretation.tsv"))

# Figure: age-stratified effect overview
png(file.path(plot_dir, "age_stratified_rbsr_effects_by_core.png"), width = 1700, height = 1000, res = 150)
op2 <- par(mfrow = c(1, 2), mar = c(8, 5, 3, 1))

pp <- per_core_stratified[per_core_stratified$scheme == "median_half", , drop = FALSE]
if (nrow(pp) > 0) {
  pp <- pp[order(pp$core, pp$stratum), ]
  lbl <- paste(pp$core, pp$stratum, sep = "\n")
  vals <- pp$rb_sr_rho_age_adjusted_linear
  cols <- ifelse(pp$stratum == "older_half", "#1f77b4", "#ff7f0e")
  bp <- barplot(vals, names.arg = lbl, las = 2, col = cols,
                ylab = "Age-adjusted Spearman rho (rb_sr vs log10 library)",
                main = "Median split: younger vs older")
  abline(h = 0, lty = 2)
  text(bp, vals, labels = paste0("n=", pp$n), pos = ifelse(vals >= 0, 3, 1), cex = 0.75)
} else {
  plot.new(); text(0.5, 0.5, "No median-half stratified results")
}

tt <- temporal_vs_plus[temporal_vs_plus$scope == "age_stratum" & temporal_vs_plus$scheme == "median_half", , drop = FALSE]
if (nrow(tt) > 0) {
  tt <- tt[order(tt$core, tt$stratum), ]
  lbl2 <- paste(tt$core, tt$stratum, sep = "\n")
  vals2 <- tt$delta_aic_temporal_minus_plus
  cols2 <- ifelse(tt$stratum == "older_half", "#1f77b4", "#ff7f0e")
  bp2 <- barplot(vals2, names.arg = lbl2, las = 2, col = cols2,
                 ylab = "ΔAIC (temporal-only minus temporal+rb_sr)",
                 main = "Model gain from adding rb_sr by stratum")
  abline(h = 0, lty = 2)
  text(bp2, vals2, labels = paste0("n=", tt$n_obs), pos = ifelse(vals2 >= 0, 3, 1), cex = 0.75)
} else {
  plot.new(); text(0.5, 0.5, "No stratified model-comparison results")
}
par(op2)
dev.off()

# Cross-core comparison figure
png(file.path(plot_dir, "multicore_rbsr_comparison.png"), width = 1600, height = 900, res = 150)
op <- par(mfrow = c(1, 2), mar = c(8, 5, 3, 1))

if (nrow(cross_summary) > 0) {
  cores <- cross_summary$core
  vals <- cross_summary$rb_sr_raw_rho
  cols <- ifelse(vals < 0, "#1f77b4", "#d62728")
  barplot(vals, names.arg = cores, las = 2, col = cols,
          ylab = "Raw Spearman rho (log10 library vs rb_sr)",
          main = "Raw rb_sr association")
  abline(h = 0, lty = 2)
} else {
  plot.new(); text(0.5, 0.5, "No cross-core summary available")
}

if (nrow(rb_effect_table) > 0) {
  ord <- seq_len(nrow(rb_effect_table))
  x <- rb_effect_table$rb_sr_model_estimate
  y <- ord
  xlim <- range(c(rb_effect_table$rb_sr_model_ci_low, rb_effect_table$rb_sr_model_ci_high, rb_effect_table$rb_sr_detrended_rho, 0), na.rm = TRUE)
  if (!all(is.finite(xlim))) xlim <- c(-1, 1)
  plot(x, y, pch = 16,
       col = ifelse(x < 0, "#1f77b4", "#d62728"),
       xlim = xlim, yaxt = "n",
       xlab = "Depth-adjusted rb_sr metric",
       ylab = "",
       main = "Depth-adjusted rb_sr effect")
  axis(2, at = y, labels = rb_effect_table$core, las = 1)
  if (all(c("rb_sr_model_ci_low", "rb_sr_model_ci_high") %in% names(rb_effect_table))) {
    segments(rb_effect_table$rb_sr_model_ci_low, y, rb_effect_table$rb_sr_model_ci_high, y, col = "#444444")
  }
  # fallback points for cores where model estimate is missing
  miss <- which(!is.finite(x) & is.finite(rb_effect_table$rb_sr_detrended_rho))
  if (length(miss) > 0) {
    points(rb_effect_table$rb_sr_detrended_rho[miss], y[miss], pch = 17, col = "#2ca02c")
  }
  abline(v = 0, lty = 2)
  legend("bottomright", legend = c("Model estimate", "Detrended rho fallback"), pch = c(16, 17), col = c("#444444", "#2ca02c"), bty = "n")
} else {
  plot.new(); text(0.5, 0.5, "No rb_sr effect table available")
}
par(op)
dev.off()

writeLines(c(
  "This multicore validation is observational and core-specific; associations are not causal proof.",
  "ST5/ST8 are analyzed with a ST13-like framework (layer aggregation, depth-aware screening/modeling).",
  "GeoB25202 is intentionally basic and rb_sr-focused after conservative technical replicate collapse.",
  "Downstream sequencing-output metrics (initial, derep, avg_leng_initial, avg_len_derep) were excluded from explanatory screening and model selection.",
  "A shared negative sign across cores supports generality of a terrigenous/mineralogical association, not a proven mechanism.",
  "Lack of support in any core can reflect missingness, lower power, or different sedimentary regimes."
), con = file.path(out_dir, "analysis_caveats.txt"))

param_tbl <- data.frame(
  parameter = names(params),
  value = vapply(params, function(x) paste(x, collapse = ","), character(1)),
  stringsAsFactors = FALSE
)
write_tsv(param_tbl, file.path(out_dir, "analysis_parameters.tsv"))

message("Multi-core rb_sr validation complete.")
message("Outputs written to: ", out_dir)
