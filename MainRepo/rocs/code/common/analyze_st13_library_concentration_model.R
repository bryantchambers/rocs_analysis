#!/usr/bin/env Rscript

# Conservative ST13-only analysis of library concentration variation.
# Guardrails:
# - single-core observational associations only (not causal proof)
# - depth/age structure is a major confounder and is explicitly modeled
# - sparse foram/SST predictors are used for sensitivity only

options(stringsAsFactors = FALSE)

suppressWarnings({
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required (base-recommended package not found).")
  }
})

# -------------------------------
# Paths and parameters
# -------------------------------
out_dir <- file.path("analysis", "st13_library_concentration_model")
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
  st13_core_label = "ST13",
  sparse_match_tolerance_cm = 3,
  response_pseudocount_rule = "min_positive_library_concentration / 2",
  xrf_layer_aggregation_primary = "median",
  xrf_layer_aggregation_sensitivity = "mean",
  depth_spline_df = 4,
  valley_rule = "residual <= 10th percentile from depth-only model",
  screening_min_n = 10,
  covariation_min_nonmissing = 30
)

# -------------------------------
# Helpers
# -------------------------------
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
  dsub <- depth[ok]
  vsub <- v[ok]
  fit <- tryCatch(
    lm(vsub ~ splines::ns(dsub, df = df)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(out)
  out[ok] <- residuals(fit)
  out
}

layer_intervals <- function(depths) {
  d <- sort(unique(depths))
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
      if (length(vals) == 0) {
        out[[nm]][i] <- NA_real_
      } else if (fun == "median") {
        out[[nm]][i] <- stats::median(vals)
      } else {
        out[[nm]][i] <- mean(vals)
      }
    }
  }
  out
}

add_ratio <- function(df, num_col, den_col, out_col, ratio_rules) {
  if (!(num_col %in% names(df)) || !(den_col %in% names(df))) {
    df[[out_col]] <- NA_real_
    ratio_rules[[length(ratio_rules) + 1]] <- data.frame(
      ratio = out_col,
      numerator = num_col,
      denominator = den_col,
      denominator_eps = NA_real_,
      note = "missing numerator or denominator column"
    )
    return(list(df = df, ratio_rules = ratio_rules))
  }
  den <- df[[den_col]]
  eps <- safe_min_positive(den, default = 1e-6) / 2
  den_adj <- ifelse(is.na(den), NA_real_, ifelse(den > 0, den, eps))
  df[[out_col]] <- df[[num_col]] / den_adj
  ratio_rules[[length(ratio_rules) + 1]] <- data.frame(
    ratio = out_col,
    numerator = num_col,
    denominator = den_col,
    denominator_eps = eps,
    note = "denominator <= 0 replaced with eps"
  )
  list(df = df, ratio_rules = ratio_rules)
}

nearest_sparse_match <- function(layers_df, sparse_df, tol_cm = 3, source_label = "sparse") {
  out <- data.frame(
    depth_in_core_cm = layers_df$depth_in_core_cm,
    stringsAsFactors = FALSE
  )

  if (nrow(sparse_df) == 0) {
    out[[paste0(source_label, "_matched_depth_cm")]] <- NA_real_
    out[[paste0(source_label, "_depth_diff_cm")]] <- NA_real_
    out[[paste0(source_label, "_matched_within_tolerance")]] <- FALSE
    return(out)
  }

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

  matched_depth <- ifelse(is.na(nearest_idx), NA_real_, sparse_depth[nearest_idx])
  within_tol <- is.finite(depth_diff) & depth_diff <= tol_cm

  out[[paste0(source_label, "_matched_depth_cm")]] <- matched_depth
  out[[paste0(source_label, "_depth_diff_cm")]] <- depth_diff
  out[[paste0(source_label, "_matched_within_tolerance")]] <- within_tol

  numeric_cols <- names(sparse_df)[vapply(sparse_df, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("depth_in_core_cm", "y_bp"))
  for (nm in numeric_cols) {
    vals <- rep(NA_real_, length(layer_depth))
    ok <- which(within_tol & !is.na(nearest_idx))
    if (length(ok) > 0) vals[ok] <- sparse_df[[nm]][nearest_idx[ok]]
    out[[paste0(source_label, "_", nm)]] <- vals
  }

  out
}

fit_primary_model <- function(df, response_col, depth_col, predictors, trend_col = "depth_in_core_cm", use_gls_if_available = TRUE) {
  use_cols <- unique(c(response_col, depth_col, trend_col, predictors))
  use_cols <- use_cols[use_cols %in% names(df)]
  dat <- df[, use_cols, drop = FALSE]
  complete_rows <- complete.cases(dat)
  dat <- dat[complete_rows, , drop = FALSE]

  if (nrow(dat) < 20) {
    return(list(model = NULL, model_type = "none", data = dat, coefficients = NULL, fit = NULL, error = "too few complete rows"))
  }

  # Standardize selected predictors
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
    return(list(model = NULL, model_type = "none", data = dat, coefficients = NULL, fit = NULL, error = "no variable with non-zero variance"))
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
        correlation = nlme::corCAR1(form = stats::as.formula(paste("~", depth_col)))
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

  # Coefficients
  coef_df <- NULL
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
    ci_low <- rep(NA_real_, nrow(coef_df))
    ci_high <- rep(NA_real_, nrow(coef_df))
    cint <- tryCatch(nlme::intervals(model, which = "coef")$coef, error = function(e) NULL)
    if (!is.null(cint)) {
      m <- match(coef_df$term, rownames(cint))
      ci_low <- cint[m, "lower"]
      ci_high <- cint[m, "upper"]
    }
    coef_df$ci_low <- ci_low
    coef_df$ci_high <- ci_high
  } else {
    sm <- summary(model)
    tt <- sm$coefficients
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

  fitted_vals <- as.numeric(stats::fitted(model))
  obs_vals <- dat[[response_col]]
  fit_df <- data.frame(
    metric = c("n_obs", "AIC", "BIC", "corr_sq_observed_vs_fitted"),
    value = c(
      nrow(dat),
      AIC(model),
      BIC(model),
      suppressWarnings(cor(obs_vals, fitted_vals, use = "complete.obs")^2)
    ),
    model_type = model_type,
    stringsAsFactors = FALSE
  )

  list(model = model, model_type = model_type, data = dat, coefficients = coef_df, fit = fit_df, predictors_z = z_names)
}

write_tsv_base <- function(df, file) {
  utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_matrix_tsv <- function(mat, file) {
  out <- cbind(variable = rownames(mat), as.data.frame(mat, check.names = FALSE))
  write_tsv_base(out, file)
}

extract_metric_value <- function(fit_df, metric_name) {
  if (is.null(fit_df) || nrow(fit_df) == 0 || !("metric" %in% names(fit_df))) return(NA_real_)
  ii <- which(fit_df$metric == metric_name)
  if (length(ii) == 0) return(NA_real_)
  as.numeric(fit_df$value[ii[1]])
}

compute_residual_acf_diagnostics <- function(model, model_id, max_lag = 3) {
  if (is.null(model)) {
    return(data.frame(
      model_id = model_id,
      residual_type = NA_character_,
      n_obs = NA_integer_,
      lag1_acf = NA_real_,
      lag2_acf = NA_real_,
      lag3_acf = NA_real_,
      ljung_box_lag = NA_integer_,
      ljung_box_p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  r <- tryCatch({
    if (inherits(model, "gls")) {
      as.numeric(residuals(model, type = "normalized"))
    } else {
      as.numeric(residuals(model))
    }
  }, error = function(e) as.numeric(residuals(model)))

  r <- r[is.finite(r)]
  if (length(r) < 5) {
    return(data.frame(
      model_id = model_id,
      residual_type = ifelse(inherits(model, "gls"), "normalized", "response"),
      n_obs = length(r),
      lag1_acf = NA_real_,
      lag2_acf = NA_real_,
      lag3_acf = NA_real_,
      ljung_box_lag = NA_integer_,
      ljung_box_p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  ac <- stats::acf(r, lag.max = max_lag, plot = FALSE, na.action = na.pass)$acf
  ac <- as.numeric(ac)
  lag_vals <- rep(NA_real_, 3)
  for (k in seq_len(3)) {
    idx <- k + 1
    if (idx <= length(ac)) lag_vals[k] <- ac[idx]
  }

  lb_lag <- min(10, max(3, floor(length(r) / 5)))
  lb <- tryCatch(stats::Box.test(r, lag = lb_lag, type = "Ljung-Box"), error = function(e) NULL)

  data.frame(
    model_id = model_id,
    residual_type = ifelse(inherits(model, "gls"), "normalized", "response"),
    n_obs = length(r),
    lag1_acf = lag_vals[1],
    lag2_acf = lag_vals[2],
    lag3_acf = lag_vals[3],
    ljung_box_lag = lb_lag,
    ljung_box_p_value = if (!is.null(lb)) lb$p.value else NA_real_,
    stringsAsFactors = FALSE
  )
}

get_hac_coefficients <- function(lm_model) {
  if (is.null(lm_model) || !inherits(lm_model, "lm")) return(NULL)
  if (!requireNamespace("sandwich", quietly = TRUE) || !requireNamespace("lmtest", quietly = TRUE)) return(NULL)

  n <- stats::nobs(lm_model)
  nw_lag <- max(1, floor(4 * (n / 100)^(2 / 9)))
  vc <- tryCatch(sandwich::NeweyWest(lm_model, lag = nw_lag, prewhite = FALSE, adjust = TRUE), error = function(e) NULL)
  if (is.null(vc)) return(NULL)

  ct <- tryCatch(lmtest::coeftest(lm_model, vcov. = vc), error = function(e) NULL)
  if (is.null(ct)) return(NULL)

  ci <- tryCatch(lmtest::coefci(lm_model, vcov. = vc), error = function(e) NULL)
  if (is.null(ci)) {
    ci <- matrix(NA_real_, nrow = nrow(ct), ncol = 2, dimnames = list(rownames(ct), c("2.5 %", "97.5 %")))
  }

  out <- data.frame(
    term = rownames(ct),
    estimate = ct[, 1],
    std_error = ct[, 2],
    stat = ct[, 3],
    p_value = ct[, 4],
    ci_low = ci[rownames(ct), 1],
    ci_high = ci[rownames(ct), 2],
    stringsAsFactors = FALSE
  )
  attr(out, "nw_lag") <- nw_lag
  out
}

get_block_bootstrap_coefficients <- function(lm_model, source_data = NULL, n_boot = 600, block_len = NULL, seed = 13) {
  if (is.null(lm_model) || !inherits(lm_model, "lm")) return(NULL)
  dat <- source_data
  if (is.null(dat)) {
    dat <- tryCatch(stats::model.frame(lm_model), error = function(e) NULL)
  }
  if (is.null(dat) || nrow(dat) < 20) return(NULL)

  n <- nrow(dat)
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
    for (s in starts) {
      idx <- c(idx, ((s - 1 + seq_len(block_len) - 1) %% n) + 1)
    }
    idx <- idx[seq_len(n)]

    fit_b <- tryCatch(stats::lm(fml, data = dat[idx, , drop = FALSE]), error = function(e) NULL)
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

  attr(out, "n_boot") <- n_boot
  attr(out, "block_len") <- block_len
  attr(out, "successful_replicates") <- sum(rowSums(is.finite(boot_mat)) > 0)
  out
}

# -------------------------------
# Read data
# -------------------------------
metadata <- read.delim(metadata_path, sep = "\t", check.names = FALSE)
xrf <- read.csv(xrf_path, check.names = FALSE)
foram <- read.delim(foram_path, sep = "\t", check.names = FALSE)
sst <- read.delim(sst_path, sep = "\t", check.names = FALSE)

# -------------------------------
# ST13 layer table (analysis unit)
# -------------------------------
meta_st13 <- metadata[metadata$core == params$st13_core_label, , drop = FALSE]
if (nrow(meta_st13) == 0) {
  stop("No ST13 rows found in metadata_v5.tsv")
}

# one observation per ST13 library layer (depth + age)
layer_key <- interaction(meta_st13$depth_in_core_cm, meta_st13$y_bp, drop = TRUE)
idx_split <- split(seq_len(nrow(meta_st13)), layer_key)

collapse_numeric <- function(v) {
  vv <- as.numeric(v)
  vv <- vv[is.finite(vv)]
  if (length(vv) == 0) return(NA_real_)
  mean(vv)
}

layers <- data.frame(
  layer_id = sprintf("ST13_layer_%03d", seq_along(idx_split)),
  depth_in_core_cm = vapply(idx_split, function(ix) collapse_numeric(meta_st13$depth_in_core_cm[ix]), numeric(1)),
  y_bp = vapply(idx_split, function(ix) collapse_numeric(meta_st13$y_bp[ix]), numeric(1)),
  library_concentration = vapply(idx_split, function(ix) collapse_numeric(meta_st13$library_concentration[ix]), numeric(1)),
  temp = vapply(idx_split, function(ix) collapse_numeric(meta_st13$temp[ix]), numeric(1)),
  mis = vapply(idx_split, function(ix) collapse_numeric(meta_st13$mis[ix]), numeric(1)),
  initial = vapply(idx_split, function(ix) collapse_numeric(meta_st13$initial[ix]), numeric(1)),
  derep = vapply(idx_split, function(ix) collapse_numeric(meta_st13$derep[ix]), numeric(1)),
  avg_leng_initial = vapply(idx_split, function(ix) collapse_numeric(meta_st13$avg_leng_initial[ix]), numeric(1)),
  avg_len_derep = vapply(idx_split, function(ix) collapse_numeric(meta_st13$avg_len_derep[ix]), numeric(1)),
  n_rows_collapsed = vapply(idx_split, length, integer(1)),
  stringsAsFactors = FALSE
)

layers <- layers[order(layers$depth_in_core_cm), ]
rownames(layers) <- NULL

min_pos <- safe_min_positive(layers$library_concentration, default = 1e-6)
pseudocount <- min_pos / 2
layers$library_concentration_log10 <- log10(layers$library_concentration + pseudocount)

# -------------------------------
# XRF aggregation by layer interval
# -------------------------------
xrf_st13 <- xrf[xrf$core == params$st13_core_label, , drop = FALSE]
if (nrow(xrf_st13) == 0) {
  stop("No ST13 rows found in combined_xrf_geochemistry_curated.csv")
}

numeric_xrf <- names(xrf_st13)[vapply(xrf_st13, is.numeric, logical(1))]
numeric_xrf <- setdiff(numeric_xrf, c("depth_in_core_cm", "y_bp"))

intervals <- layer_intervals(layers$depth_in_core_cm)
xrf_med <- aggregate_xrf_to_layers(xrf_st13, intervals, numeric_cols = numeric_xrf, fun = "median")
xrf_mean <- aggregate_xrf_to_layers(xrf_st13, intervals, numeric_cols = numeric_xrf, fun = "mean")

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

ratio_rules <- list()
for (r in ratio_defs) {
  rr <- add_ratio(xrf_med, r[1], r[2], r[3], ratio_rules)
  xrf_med <- rr$df
  ratio_rules <- rr$ratio_rules
}

ratio_rules_mean <- list()
for (r in ratio_defs) {
  rr <- add_ratio(xrf_mean, r[1], r[2], paste0(r[3], "_meanagg"), ratio_rules_mean)
  xrf_mean <- rr$df
  ratio_rules_mean <- rr$ratio_rules
}

ratio_rule_table <- do.call(rbind, c(ratio_rules, ratio_rules_mean))

names(xrf_med)[names(xrf_med) %in% c("n_xrf_rows")] <- "n_xrf_rows_median"
names(xrf_mean)[names(xrf_mean) %in% c("n_xrf_rows")] <- "n_xrf_rows_mean"

for (nm in setdiff(names(xrf_mean), "depth_in_core_cm")) {
  if (!(nm %in% names(xrf_med))) next
  names(xrf_mean)[names(xrf_mean) == nm] <- paste0(nm, "_meanagg")
}

layer_tbl <- merge(layers, xrf_med, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
layer_tbl <- merge(layer_tbl, xrf_mean, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
layer_tbl <- layer_tbl[order(layer_tbl$depth_in_core_cm), ]
rownames(layer_tbl) <- NULL

# -------------------------------
# Sparse foram/SST nearest-depth matching
# -------------------------------
foram_st13 <- foram[foram$core == params$st13_core_label, , drop = FALSE]
sst_st13 <- sst[sst$core == params$st13_core_label, , drop = FALSE]

foram_match <- nearest_sparse_match(layer_tbl, foram_st13, tol_cm = params$sparse_match_tolerance_cm, source_label = "foram")
sst_match <- nearest_sparse_match(layer_tbl, sst_st13, tol_cm = params$sparse_match_tolerance_cm, source_label = "sst")

layer_tbl <- merge(layer_tbl, foram_match, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
layer_tbl <- merge(layer_tbl, sst_match, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
layer_tbl <- layer_tbl[order(layer_tbl$depth_in_core_cm), ]

# -------------------------------
# Candidate variable sets
# -------------------------------
metadata_env <- intersect(c("temp", "mis"), names(layer_tbl))
metadata_tech_downstream <- intersect(c("initial", "derep", "avg_leng_initial", "avg_len_derep"), names(layer_tbl))

xrf_raw_candidates <- intersect(numeric_xrf, names(layer_tbl))
xrf_ratio_candidates <- intersect(c("ba_ti", "p_ti", "br_ti", "si_ti", "ca_ti", "al_si", "ti_al", "fe_al", "rb_sr", "mn_fe"), names(layer_tbl))

foram_candidates <- grep("^foram_", names(layer_tbl), value = TRUE)
foram_candidates <- setdiff(foram_candidates, c("foram_matched_depth_cm", "foram_depth_diff_cm", "foram_matched_within_tolerance"))
sst_candidates <- grep("^sst_", names(layer_tbl), value = TRUE)
sst_candidates <- setdiff(sst_candidates, c("sst_matched_depth_cm", "sst_depth_diff_cm", "sst_matched_within_tolerance"))

family_map <- c(
  setNames(rep("metadata_env", length(metadata_env)), metadata_env),
  setNames(rep("xrf_raw", length(xrf_raw_candidates)), xrf_raw_candidates),
  setNames(rep("xrf_ratio", length(xrf_ratio_candidates)), xrf_ratio_candidates),
  setNames(rep("foram_sparse", length(foram_candidates)), foram_candidates),
  setNames(rep("sst_sparse", length(sst_candidates)), sst_candidates)
)

all_explanatory_candidates <- names(family_map)

downstream_exclusion_reasons <- c(
  initial = "Downstream sequencing-output metric; excluded from explanatory modeling.",
  derep = "Downstream sequencing-output metric; excluded from explanatory modeling.",
  avg_leng_initial = "Downstream sequencing-output metric; excluded from explanatory modeling.",
  avg_len_derep = "Downstream sequencing-output metric; excluded from explanatory modeling."
)

# -------------------------------
# Coverage table
# -------------------------------
coverage <- data.frame(
  variable = c("library_concentration", "library_concentration_log10", all_explanatory_candidates),
  family = c("response", "response", unname(family_map)),
  n_non_missing = NA_integer_,
  pct_non_missing = NA_real_,
  n_unique = NA_integer_,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(coverage))) {
  v <- coverage$variable[i]
  x <- layer_tbl[[v]]
  ok <- is.finite(x)
  coverage$n_non_missing[i] <- sum(ok)
  coverage$pct_non_missing[i] <- round(100 * mean(ok), 2)
  coverage$n_unique[i] <- length(unique(x[ok]))
}

downstream_qc <- data.frame(
  variable = metadata_tech_downstream,
  classification = "downstream_sequencing_metric",
  exclusion_reason = unname(downstream_exclusion_reasons[metadata_tech_downstream]),
  n_non_missing = NA_integer_,
  pct_non_missing = NA_real_,
  n_unique = NA_integer_,
  stringsAsFactors = FALSE
)
if (nrow(downstream_qc) > 0) {
  for (i in seq_len(nrow(downstream_qc))) {
    x <- layer_tbl[[downstream_qc$variable[i]]]
    ok <- is.finite(x)
    downstream_qc$n_non_missing[i] <- sum(ok)
    downstream_qc$pct_non_missing[i] <- round(100 * mean(ok), 2)
    downstream_qc$n_unique[i] <- length(unique(x[ok]))
  }
}

# -------------------------------
# Screening correlations (raw + detrended)
# -------------------------------
response <- layer_tbl$library_concentration_log10
depth <- layer_tbl$depth_in_core_cm
res_response <- depth_residuals(response, depth, df = params$depth_spline_df)

screen <- data.frame(
  variable = all_explanatory_candidates,
  family = unname(family_map[all_explanatory_candidates]),
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

screen$q_raw_bh_family <- NA_real_
screen$q_detrended_bh_family <- NA_real_
for (fam in unique(screen$family)) {
  i_f <- which(screen$family == fam & is.finite(screen$p_raw))
  if (length(i_f) > 0) screen$q_raw_bh_family[i_f] <- p.adjust(screen$p_raw[i_f], method = "BH")
  j_f <- which(screen$family == fam & is.finite(screen$p_detrended))
  if (length(j_f) > 0) screen$q_detrended_bh_family[j_f] <- p.adjust(screen$p_detrended[j_f], method = "BH")
}

screen <- screen[order(-abs(screen$rho_detrended), screen$p_detrended), ]

# -------------------------------
# Co-variation among explanatory variables
# -------------------------------
covar_vars <- coverage$variable[
  coverage$variable %in% all_explanatory_candidates &
    coverage$n_non_missing >= params$covariation_min_nonmissing &
    coverage$n_unique >= 4
]

covar_mat <- matrix(NA_real_, nrow = length(covar_vars), ncol = length(covar_vars), dimnames = list(covar_vars, covar_vars))
for (i in seq_along(covar_vars)) {
  for (j in seq_along(covar_vars)) {
    xi <- as.numeric(layer_tbl[[covar_vars[i]]])
    xj <- as.numeric(layer_tbl[[covar_vars[j]]])
    covar_mat[i, j] <- suppressWarnings(cor(xi, xj, method = "spearman", use = "pairwise.complete.obs"))
  }
}

covar_long <- data.frame(
  var1 = character(0), var2 = character(0), spearman_rho = numeric(0), abs_rho = numeric(0),
  stringsAsFactors = FALSE
)
if (length(covar_vars) > 1) {
  for (i in seq_len(length(covar_vars) - 1)) {
    for (j in (i + 1):length(covar_vars)) {
      rho <- covar_mat[i, j]
      covar_long <- rbind(covar_long, data.frame(
        var1 = covar_vars[i],
        var2 = covar_vars[j],
        spearman_rho = rho,
        abs_rho = abs(rho),
        stringsAsFactors = FALSE
      ))
    }
  }
  covar_long <- covar_long[order(-covar_long$abs_rho), ]
}

# -------------------------------
# Choose interpretable predictor set for primary model
# -------------------------------
pick_best_from_group <- function(group_vars, screen_df, min_n = 30) {
  sub <- screen_df[screen_df$variable %in% group_vars & is.finite(screen_df$rho_detrended) & screen_df$n_detrended >= min_n, , drop = FALSE]
  if (nrow(sub) == 0) {
    sub2 <- screen_df[screen_df$variable %in% group_vars & is.finite(screen_df$rho_raw) & screen_df$n_raw >= min_n, , drop = FALSE]
    if (nrow(sub2) == 0) return(NA_character_)
    return(sub2$variable[which.max(abs(sub2$rho_raw))[1]])
  }
  sub$ord <- abs(sub$rho_detrended)
  sub <- sub[order(-sub$ord), ]
  sub$variable[1]
}

prod_var <- pick_best_from_group(intersect(c("ba_ti", "p_ti", "br_ti"), names(layer_tbl)), screen, min_n = 40)
terrig_var <- pick_best_from_group(intersect(c("al_si", "ti_al", "fe_al", "rb_sr", "si_ti"), names(layer_tbl)), screen, min_n = 40)
ca_var <- if ("ca_ti" %in% names(layer_tbl)) "ca_ti" else NA_character_
mnfe_var <- if ("mn_fe" %in% names(layer_tbl)) "mn_fe" else NA_character_

selected_predictors <- unique(na.omit(c(prod_var, terrig_var, ca_var, mnfe_var)))

selected_tbl <- data.frame(
  process_group = c("productivity_or_organic", "terrigenous_or_mineralogical", "carbonate", "redox"),
  selected_predictor = c(prod_var, terrig_var, ca_var, mnfe_var),
  stringsAsFactors = FALSE
)

# -------------------------------
# Primary model (depth trend + autocorrelation if available)
# -------------------------------
primary_fit <- fit_primary_model(
  df = layer_tbl,
  response_col = "library_concentration_log10",
  depth_col = "depth_in_core_cm",
  trend_col = "depth_in_core_cm",
  predictors = selected_predictors,
  use_gls_if_available = TRUE
)

primary_fitted_table <- data.frame()
if (!is.null(primary_fit$model)) {
  dat <- primary_fit$data
  primary_fitted_table <- data.frame(
    depth_in_core_cm = dat$depth_in_core_cm,
    y_bp = if ("y_bp" %in% names(dat)) dat$y_bp else NA_real_,
    observed_log10_library_concentration = dat$library_concentration_log10,
    fitted_log10_library_concentration = as.numeric(stats::fitted(primary_fit$model)),
    residual = as.numeric(stats::residuals(primary_fit$model)),
    stringsAsFactors = FALSE
  )
}

# -------------------------------
# Autocorrelation-robust inference for primary predictor structure
# -------------------------------
autocorr_summary <- data.frame()
autocorr_fit_metrics <- data.frame()
autocorr_diagnostics <- data.frame()

# Method 1 (preferred): GLS with continuous AR(1) correlation in depth order.
auto_gls <- fit_primary_model(
  df = layer_tbl,
  response_col = "library_concentration_log10",
  depth_col = "depth_in_core_cm",
  trend_col = "depth_in_core_cm",
  predictors = selected_predictors,
  use_gls_if_available = TRUE
)

if (!is.null(auto_gls$model)) {
  auto_gls_id <- ifelse(auto_gls$model_type == "gls_corCAR1", "gls_corCAR1_primary", "lm_primary_fallback")
  auto_gls_note <- ifelse(
    auto_gls$model_type == "gls_corCAR1",
    "Autocorrelation handled via nlme::gls with corCAR1(~ depth_in_core_cm).",
    "GLS unavailable or failed; no explicit residual autocorrelation term in this fit."
  )
  gls_coef <- auto_gls$coefficients
  gls_coef$model_id <- auto_gls_id
  gls_coef$autocorrelation_handling <- auto_gls_note
  autocorr_summary <- rbind(autocorr_summary, gls_coef[, c("model_id", "autocorrelation_handling", "term", "estimate", "std_error", "stat", "p_value", "ci_low", "ci_high")])

  autocorr_fit_metrics <- rbind(autocorr_fit_metrics, data.frame(
    model_id = auto_gls_id,
    status = "fit",
    model_type = auto_gls$model_type,
    autocorrelation_handling = auto_gls_note,
    n_obs = extract_metric_value(auto_gls$fit, "n_obs"),
    AIC = extract_metric_value(auto_gls$fit, "AIC"),
    BIC = extract_metric_value(auto_gls$fit, "BIC"),
    corr_sq_observed_vs_fitted = extract_metric_value(auto_gls$fit, "corr_sq_observed_vs_fitted"),
    note = "Primary predictor structure",
    stringsAsFactors = FALSE
  ))
  autocorr_diagnostics <- rbind(autocorr_diagnostics, compute_residual_acf_diagnostics(auto_gls$model, auto_gls_id))
} else {
  autocorr_fit_metrics <- rbind(autocorr_fit_metrics, data.frame(
    model_id = "gls_corCAR1_primary",
    status = "not_fit",
    model_type = "none",
    autocorrelation_handling = "GLS corCAR1 not fit",
    n_obs = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    corr_sq_observed_vs_fitted = NA_real_,
    note = auto_gls$error,
    stringsAsFactors = FALSE
  ))
}

# Method 2 (defensible fallback/comparison): LM with HAC/Newey-West SE if available.
auto_lm <- fit_primary_model(
  df = layer_tbl,
  response_col = "library_concentration_log10",
  depth_col = "depth_in_core_cm",
  trend_col = "depth_in_core_cm",
  predictors = selected_predictors,
  use_gls_if_available = FALSE
)

if (!is.null(auto_lm$model)) {
  hac_coef <- get_hac_coefficients(auto_lm$model)
  if (!is.null(hac_coef)) {
    nw_lag <- attr(hac_coef, "nw_lag")
    hac_coef$model_id <- "lm_newey_west_primary"
    hac_coef$autocorrelation_handling <- paste0("LM with Newey-West HAC SE (lag=", nw_lag, ")")
    autocorr_summary <- rbind(autocorr_summary, hac_coef[, c("model_id", "autocorrelation_handling", "term", "estimate", "std_error", "stat", "p_value", "ci_low", "ci_high")])
    autocorr_fit_metrics <- rbind(autocorr_fit_metrics, data.frame(
      model_id = "lm_newey_west_primary",
      status = "fit",
      model_type = "lm_hac",
      autocorrelation_handling = paste0("HAC robust uncertainty (Newey-West lag=", nw_lag, ")"),
      n_obs = extract_metric_value(auto_lm$fit, "n_obs"),
      AIC = extract_metric_value(auto_lm$fit, "AIC"),
      BIC = extract_metric_value(auto_lm$fit, "BIC"),
      corr_sq_observed_vs_fitted = extract_metric_value(auto_lm$fit, "corr_sq_observed_vs_fitted"),
      note = "Primary predictor structure; HAC uncertainty for serial dependence",
      stringsAsFactors = FALSE
    ))
  } else {
    bb_coef <- get_block_bootstrap_coefficients(auto_lm$model, source_data = auto_lm$data, n_boot = 600)
    if (!is.null(bb_coef)) {
      bb_nboot <- attr(bb_coef, "n_boot")
      bb_block <- attr(bb_coef, "block_len")
      bb_success <- attr(bb_coef, "successful_replicates")
      bb_coef$model_id <- "lm_block_bootstrap_primary"
      bb_coef$autocorrelation_handling <- paste0("LM with moving block bootstrap by depth order (n_boot=", bb_nboot, ", block_len=", bb_block, ")")
      autocorr_summary <- rbind(autocorr_summary, bb_coef[, c("model_id", "autocorrelation_handling", "term", "estimate", "std_error", "stat", "p_value", "ci_low", "ci_high")])
      autocorr_fit_metrics <- rbind(autocorr_fit_metrics, data.frame(
        model_id = "lm_block_bootstrap_primary",
        status = "fit",
        model_type = "lm_block_bootstrap",
        autocorrelation_handling = paste0("Block bootstrap uncertainty (block_len=", bb_block, ")"),
        n_obs = extract_metric_value(auto_lm$fit, "n_obs"),
        AIC = extract_metric_value(auto_lm$fit, "AIC"),
        BIC = extract_metric_value(auto_lm$fit, "BIC"),
        corr_sq_observed_vs_fitted = extract_metric_value(auto_lm$fit, "corr_sq_observed_vs_fitted"),
        note = paste0("Primary predictor structure; bootstrap successes=", bb_success, "/", bb_nboot),
        stringsAsFactors = FALSE
      ))
    } else {
      lm_coef <- auto_lm$coefficients
      lm_coef$model_id <- "lm_conventional_primary"
      lm_coef$autocorrelation_handling <- "LM fit; HAC and block-bootstrap robust uncertainty unavailable (conventional SE shown)"
      autocorr_summary <- rbind(autocorr_summary, lm_coef[, c("model_id", "autocorrelation_handling", "term", "estimate", "std_error", "stat", "p_value", "ci_low", "ci_high")])
      autocorr_fit_metrics <- rbind(autocorr_fit_metrics, data.frame(
        model_id = "lm_conventional_primary",
        status = "fit",
        model_type = "lm",
        autocorrelation_handling = "No explicit autocorrelation correction in SE",
        n_obs = extract_metric_value(auto_lm$fit, "n_obs"),
        AIC = extract_metric_value(auto_lm$fit, "AIC"),
        BIC = extract_metric_value(auto_lm$fit, "BIC"),
        corr_sq_observed_vs_fitted = extract_metric_value(auto_lm$fit, "corr_sq_observed_vs_fitted"),
        note = "Fallback only; interpret uncertainty cautiously if residual autocorrelation remains",
        stringsAsFactors = FALSE
      ))
    }
  }
  autocorr_diagnostics <- rbind(autocorr_diagnostics, compute_residual_acf_diagnostics(auto_lm$model, "lm_primary"))
} else {
  autocorr_fit_metrics <- rbind(autocorr_fit_metrics, data.frame(
    model_id = "lm_primary",
    status = "not_fit",
    model_type = "none",
    autocorrelation_handling = "LM fallback not fit",
    n_obs = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    corr_sq_observed_vs_fitted = NA_real_,
    note = auto_lm$error,
    stringsAsFactors = FALSE
  ))
}

# -------------------------------
# Geochemical-axis stability across collinear formulations
# -------------------------------
prod_for_stability <- if ("ba_ti" %in% names(layer_tbl)) "ba_ti" else prod_var
core_covariates <- unique(na.omit(c(prod_for_stability, "ca_ti", "mn_fe")))

raw_axis_candidate <- NA_character_
for (rv in c("rb", "si", "zr")) {
  if (rv %in% names(layer_tbl)) {
    ok <- is.finite(layer_tbl[[rv]])
    if (sum(ok) >= 30 && length(unique(layer_tbl[[rv]][ok])) >= 4) {
      raw_axis_candidate <- rv
      break
    }
  }
}

axis_candidates <- unique(na.omit(c("rb_sr", "al_si", "si_ti", "fe_al", raw_axis_candidate)))
geochem_stability <- data.frame(
  mineral_representation = character(0),
  representation_family = character(0),
  status = character(0),
  model_type = character(0),
  autocorrelation_handling = character(0),
  n_obs = numeric(0),
  mineral_term = character(0),
  estimate = numeric(0),
  estimate_sign = character(0),
  std_error = numeric(0),
  p_value = numeric(0),
  ci_low = numeric(0),
  ci_high = numeric(0),
  AIC = numeric(0),
  BIC = numeric(0),
  corr_sq_observed_vs_fitted = numeric(0),
  note = character(0),
  stringsAsFactors = FALSE
)

for (mv in axis_candidates) {
  model_predictors <- unique(c(mv, core_covariates))
  model_predictors <- model_predictors[model_predictors %in% names(layer_tbl)]
  fit_mv <- fit_primary_model(
    df = layer_tbl,
    response_col = "library_concentration_log10",
    depth_col = "depth_in_core_cm",
    trend_col = "depth_in_core_cm",
    predictors = model_predictors,
    use_gls_if_available = TRUE
  )

  if (is.null(fit_mv$model)) {
    geochem_stability <- rbind(geochem_stability, data.frame(
      mineral_representation = mv,
      representation_family = ifelse(mv %in% c("rb", "si", "zr"), "raw_element", "ratio"),
      status = "not_fit",
      model_type = "none",
      autocorrelation_handling = "not_fit",
      n_obs = NA_real_,
      mineral_term = paste0("z_", mv),
      estimate = NA_real_,
      estimate_sign = NA_character_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      corr_sq_observed_vs_fitted = NA_real_,
      note = fit_mv$error,
      stringsAsFactors = FALSE
    ))
    next
  }

  term_name <- paste0("z_", mv)
  cf <- fit_mv$coefficients
  row_mv <- cf[cf$term == term_name, , drop = FALSE]
  if (nrow(row_mv) == 0) {
    geochem_stability <- rbind(geochem_stability, data.frame(
      mineral_representation = mv,
      representation_family = ifelse(mv %in% c("rb", "si", "zr"), "raw_element", "ratio"),
      status = "fit_but_term_not_estimable",
      model_type = fit_mv$model_type,
      autocorrelation_handling = ifelse(fit_mv$model_type == "gls_corCAR1", "gls_corCAR1", "lm_no_correlation_term"),
      n_obs = extract_metric_value(fit_mv$fit, "n_obs"),
      mineral_term = term_name,
      estimate = NA_real_,
      estimate_sign = NA_character_,
      std_error = NA_real_,
      p_value = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      AIC = extract_metric_value(fit_mv$fit, "AIC"),
      BIC = extract_metric_value(fit_mv$fit, "BIC"),
      corr_sq_observed_vs_fitted = extract_metric_value(fit_mv$fit, "corr_sq_observed_vs_fitted"),
      note = "Mineralogical term dropped or not estimable (possible collinearity/zero variance).",
      stringsAsFactors = FALSE
    ))
    next
  }

  geochem_stability <- rbind(geochem_stability, data.frame(
    mineral_representation = mv,
    representation_family = ifelse(mv %in% c("rb", "si", "zr"), "raw_element", "ratio"),
    status = "fit",
    model_type = fit_mv$model_type,
    autocorrelation_handling = ifelse(fit_mv$model_type == "gls_corCAR1", "gls_corCAR1", "lm_no_correlation_term"),
    n_obs = extract_metric_value(fit_mv$fit, "n_obs"),
    mineral_term = term_name,
    estimate = row_mv$estimate[1],
    estimate_sign = ifelse(is.na(row_mv$estimate[1]), NA_character_, ifelse(row_mv$estimate[1] > 0, "positive", ifelse(row_mv$estimate[1] < 0, "negative", "zero"))),
    std_error = row_mv$std_error[1],
    p_value = row_mv$p_value[1],
    ci_low = row_mv$ci_low[1],
    ci_high = row_mv$ci_high[1],
    AIC = extract_metric_value(fit_mv$fit, "AIC"),
    BIC = extract_metric_value(fit_mv$fit, "BIC"),
    corr_sq_observed_vs_fitted = extract_metric_value(fit_mv$fit, "corr_sq_observed_vs_fitted"),
    note = paste0("Depth spline + ", paste(setdiff(model_predictors, mv), collapse = " + "), "; mineralogical term represented by ", mv, "."),
    stringsAsFactors = FALSE
  ))
}

png(file.path(plot_dir, "geochemical_axis_stability_coefficients.png"), width = 1400, height = 900, res = 150)
stable_plot_df <- geochem_stability[geochem_stability$status == "fit" & is.finite(geochem_stability$estimate), , drop = FALSE]
if (nrow(stable_plot_df) > 0) {
  ord <- order(stable_plot_df$estimate)
  d <- stable_plot_df[ord, , drop = FALSE]
  x <- d$estimate
  y <- seq_len(nrow(d))
  xlim <- range(c(d$ci_low, d$ci_high, 0), na.rm = TRUE)
  plot(x, y,
       xlim = xlim,
       yaxt = "n",
       pch = 16,
       col = "#1f77b4",
       xlab = "Standardized coefficient (mineralogical term)",
       ylab = "",
       main = "Geochemical-axis stability across collinear formulations")
  axis(2, at = y, labels = d$mineral_representation, las = 1)
  segments(d$ci_low, y, d$ci_high, y, col = "#1f77b4", lwd = 2)
  abline(v = 0, lty = 2, col = "#d62728")
} else {
  plot.new(); text(0.5, 0.5, "No geochemical stability models successfully estimated")
}
dev.off()

# -------------------------------
# Sensitivity models
# -------------------------------
sensitivity_overview <- data.frame(
  sensitivity_model = character(0),
  status = character(0),
  n_obs = numeric(0),
  model_type = character(0),
  note = character(0),
  stringsAsFactors = FALSE
)

# Sensitivity 1: age trend instead of depth trend
sens_age <- fit_primary_model(
  df = layer_tbl,
  response_col = "library_concentration_log10",
  depth_col = "depth_in_core_cm",
  trend_col = "y_bp",
  predictors = selected_predictors,
  use_gls_if_available = TRUE
)
sensitivity_overview <- rbind(sensitivity_overview, data.frame(
  sensitivity_model = "age_trend_instead_of_depth",
  status = ifelse(is.null(sens_age$model), "not_fit", "fit"),
  n_obs = ifelse(is.null(sens_age$model), NA, nrow(sens_age$data)),
  model_type = ifelse(is.null(sens_age$model), "none", sens_age$model_type),
  note = ifelse(is.null(sens_age$model), sens_age$error, "y_bp spline trend"),
  stringsAsFactors = FALSE
))

# Sensitivity 2: mean XRF aggregation instead of median
mean_predictors <- paste0(selected_predictors, "_meanagg")
mean_predictors <- mean_predictors[mean_predictors %in% names(layer_tbl)]
sens_mean <- fit_primary_model(
  df = layer_tbl,
  response_col = "library_concentration_log10",
  depth_col = "depth_in_core_cm",
  trend_col = "depth_in_core_cm",
  predictors = mean_predictors,
  use_gls_if_available = TRUE
)
sensitivity_overview <- rbind(sensitivity_overview, data.frame(
  sensitivity_model = "mean_xrf_aggregation",
  status = ifelse(is.null(sens_mean$model), "not_fit", "fit"),
  n_obs = ifelse(is.null(sens_mean$model), NA, nrow(sens_mean$data)),
  model_type = ifelse(is.null(sens_mean$model), "none", sens_mean$model_type),
  note = ifelse(is.null(sens_mean$model), sens_mean$error, "XRF layer summaries use means"),
  stringsAsFactors = FALSE
))

# Sensitivity 3: add one SST + one foram predictor, if feasible
sst_pref <- intersect(c("sst_sst_uk37_alkenone", "sst_sst_mgca_jonkers_2013", "sst_sst_mgca_kozdon_2009"), names(layer_tbl))
foram_pref <- intersect(c("foram_g_bulloides_pct", "foram_polar_planktonic_spp_pct", "foram_n_pachyderma_pct"), names(layer_tbl))
sens_sparse <- NULL
if (length(sst_pref) > 0 && length(foram_pref) > 0) {
  sparse_predictors <- unique(c(selected_predictors, sst_pref[1], foram_pref[1]))
  sens_sparse <- fit_primary_model(
    df = layer_tbl,
    response_col = "library_concentration_log10",
    depth_col = "depth_in_core_cm",
    trend_col = "depth_in_core_cm",
    predictors = sparse_predictors,
    use_gls_if_available = FALSE
  )
  sparse_n <- ifelse(is.null(sens_sparse$model), 0, nrow(sens_sparse$data))
  # conservative threshold for sparse sensitivity fit
  if (!is.null(sens_sparse$model) && sparse_n < 20) {
    sens_sparse <- NULL
    sensitivity_overview <- rbind(sensitivity_overview, data.frame(
      sensitivity_model = "add_sparse_sst_and_foram",
      status = "not_fit",
      n_obs = sparse_n,
      model_type = "none",
      note = "too few complete rows (<20) with sparse predictors",
      stringsAsFactors = FALSE
    ))
  } else {
    sensitivity_overview <- rbind(sensitivity_overview, data.frame(
      sensitivity_model = "add_sparse_sst_and_foram",
      status = ifelse(is.null(sens_sparse$model), "not_fit", "fit"),
      n_obs = ifelse(is.null(sens_sparse$model), NA, nrow(sens_sparse$data)),
      model_type = ifelse(is.null(sens_sparse$model), "none", sens_sparse$model_type),
      note = ifelse(is.null(sens_sparse$model), sens_sparse$error, paste0("added ", sst_pref[1], " + ", foram_pref[1])),
      stringsAsFactors = FALSE
    ))
  }
} else {
  sensitivity_overview <- rbind(sensitivity_overview, data.frame(
    sensitivity_model = "add_sparse_sst_and_foram",
    status = "not_fit",
    n_obs = NA,
    model_type = "none",
    note = "required sparse predictors not available",
    stringsAsFactors = FALSE
  ))
}

# -------------------------------
# Valley analysis (secondary)
# -------------------------------
base_depth_dat <- layer_tbl[, c("layer_id", "depth_in_core_cm", "y_bp", "library_concentration", "library_concentration_log10"), drop = FALSE]
base_depth_dat <- base_depth_dat[complete.cases(base_depth_dat[, c("depth_in_core_cm", "library_concentration_log10")]), ]

valley_tbl <- data.frame()
valley_compare <- data.frame()
if (nrow(base_depth_dat) >= 20) {
  baseline <- lm(library_concentration_log10 ~ splines::ns(depth_in_core_cm, df = params$depth_spline_df), data = base_depth_dat)
  base_depth_dat$baseline_fitted <- as.numeric(fitted(baseline))
  base_depth_dat$baseline_residual <- as.numeric(residuals(baseline))
  thr <- as.numeric(quantile(base_depth_dat$baseline_residual, probs = 0.10, na.rm = TRUE))
  base_depth_dat$valley_candidate <- base_depth_dat$baseline_residual <= thr

  # Attach selected predictors for context
  valley_tbl <- merge(base_depth_dat, layer_tbl[, c("layer_id", selected_predictors), drop = FALSE], by = "layer_id", all.x = TRUE)
  valley_tbl <- valley_tbl[order(valley_tbl$depth_in_core_cm), ]

  # simple descriptive comparison
  vars_for_compare <- unique(c(selected_predictors, "library_concentration", "library_concentration_log10"))
  comp <- data.frame()
  for (v in vars_for_compare) {
    if (!(v %in% names(valley_tbl))) next
    vv <- valley_tbl[[v]]
    a <- vv[valley_tbl$valley_candidate & is.finite(vv)]
    b <- vv[!valley_tbl$valley_candidate & is.finite(vv)]
    comp <- rbind(comp, data.frame(
      variable = v,
      valley_n = length(a),
      non_valley_n = length(b),
      valley_median = ifelse(length(a) > 0, median(a), NA_real_),
      non_valley_median = ifelse(length(b) > 0, median(b), NA_real_),
      valley_iqr = ifelse(length(a) > 1, IQR(a), NA_real_),
      non_valley_iqr = ifelse(length(b) > 1, IQR(b), NA_real_),
      stringsAsFactors = FALSE
    ))
  }
  valley_compare <- comp
}

# -------------------------------
# Plot outputs (base R only)
# -------------------------------
png(file.path(plot_dir, "st13_library_concentration_vs_depth.png"), width = 1200, height = 900, res = 150)
plot(layer_tbl$depth_in_core_cm, layer_tbl$library_concentration_log10,
     pch = 16, col = "#1f77b4", cex = 0.8,
     xlab = "Depth in core (cm)", ylab = "log10(library concentration + pseudocount)",
     main = "ST13 library concentration vs depth")
if (nrow(base_depth_dat) >= 20) {
  ord <- order(base_depth_dat$depth_in_core_cm)
  lines(base_depth_dat$depth_in_core_cm[ord], base_depth_dat$baseline_fitted[ord], col = "#d62728", lwd = 2)
}
legend("topright", legend = c("Observed", "Depth-only trend"), col = c("#1f77b4", "#d62728"), pch = c(16, NA), lty = c(NA, 1), bty = "n")
dev.off()

png(file.path(plot_dir, "st13_library_concentration_vs_age.png"), width = 1200, height = 900, res = 150)
plot(layer_tbl$y_bp, layer_tbl$library_concentration_log10,
     pch = 16, col = "#2ca02c", cex = 0.8,
     xlab = "Age (y BP)", ylab = "log10(library concentration + pseudocount)",
     main = "ST13 library concentration vs age")
dev.off()

top_vars <- screen$variable[is.finite(screen$rho_detrended)]
top_vars <- top_vars[!grepl("^(foram_|sst_)", top_vars)]
top_vars <- unique(top_vars)[seq_len(min(4, length(unique(top_vars))))]

png(file.path(plot_dir, "top_explanatory_variables_vs_library_concentration.png"), width = 1600, height = 1200, res = 150)
if (length(top_vars) > 0) {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (v in top_vars) {
    x <- layer_tbl[[v]]
    y <- layer_tbl$library_concentration_log10
    plot(x, y, pch = 16, col = "#9467bd", cex = 0.75,
         xlab = v, ylab = "log10(library concentration + pseudocount)", main = v)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 10) {
      abline(lm(y[ok] ~ x[ok]), col = "#d62728", lwd = 1.5)
    }
  }
} else {
  plot.new(); text(0.5, 0.5, "No top variables available for plotting")
}
dev.off()

png(file.path(plot_dir, "predictor_covariation_heatmap.png"), width = 1300, height = 1300, res = 150)
if (length(covar_vars) >= 2) {
  cm <- covar_mat
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  hc <- hclust(as.dist(1 - abs(cm)))
  ord <- hc$order
  cm2 <- cm[ord, ord]
  op <- par(mar = c(8, 8, 3, 2))
  image(1:nrow(cm2), 1:ncol(cm2), t(cm2[nrow(cm2):1, ]), col = colorRampPalette(c("#313695", "white", "#a50026"))(200),
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "Predictor co-variation (Spearman)")
  axis(1, at = 1:ncol(cm2), labels = colnames(cm2), las = 2, cex.axis = 0.5)
  axis(2, at = 1:nrow(cm2), labels = rev(rownames(cm2)), las = 2, cex.axis = 0.5)
  par(op)
} else {
  plot.new(); text(0.5, 0.5, "Not enough variables for heatmap")
}
dev.off()

png(file.path(plot_dir, "primary_model_observed_vs_fitted.png"), width = 1200, height = 900, res = 150)
if (nrow(primary_fitted_table) > 0) {
  plot(primary_fitted_table$fitted_log10_library_concentration,
       primary_fitted_table$observed_log10_library_concentration,
       pch = 16, col = "#1f77b4", cex = 0.8,
       xlab = "Fitted log10(library concentration)",
       ylab = "Observed log10(library concentration)",
       main = "Primary model: observed vs fitted")
  abline(a = 0, b = 1, col = "#d62728", lwd = 2)
} else {
  plot.new(); text(0.5, 0.5, "Primary model not available")
}
dev.off()

png(file.path(plot_dir, "valley_diagnostic_depth_plot.png"), width = 1200, height = 900, res = 150)
if (nrow(valley_tbl) > 0) {
  cols <- ifelse(valley_tbl$valley_candidate, "#d62728", "#1f77b4")
  plot(valley_tbl$depth_in_core_cm, valley_tbl$library_concentration_log10,
       pch = 16, col = cols, cex = 0.9,
       xlab = "Depth in core (cm)", ylab = "log10(library concentration + pseudocount)",
       main = "Depth trend and valley candidates")
  ord <- order(valley_tbl$depth_in_core_cm)
  lines(valley_tbl$depth_in_core_cm[ord], valley_tbl$baseline_fitted[ord], col = "black", lwd = 2)
  legend("topright", legend = c("Non-valley", "Valley candidate", "Depth-only fitted"),
         col = c("#1f77b4", "#d62728", "black"), pch = c(16, 16, NA), lty = c(NA, NA, 1), bty = "n")
} else {
  plot.new(); text(0.5, 0.5, "Valley model not available")
}
dev.off()

# -------------------------------
# Write outputs
# -------------------------------
write_tsv_base(
  data.frame(parameter = names(params), value = unlist(params), stringsAsFactors = FALSE),
  file.path(out_dir, "analysis_parameters.tsv")
)

write_tsv_base(layer_tbl, file.path(out_dir, "st13_layer_level_analysis_table.tsv"))
write_tsv_base(coverage, file.path(out_dir, "variable_coverage_table.tsv"))
write_tsv_base(downstream_qc, file.path(out_dir, "downstream_sequencing_metrics_excluded_from_explanatory_modeling.tsv"))
write_tsv_base(screen, file.path(out_dir, "screening_correlations_raw_and_detrended.tsv"))
write_matrix_tsv(covar_mat, file.path(out_dir, "predictor_covariation_matrix.tsv"))
write_tsv_base(covar_long, file.path(out_dir, "predictor_covariation_pairs.tsv"))
write_tsv_base(selected_tbl, file.path(out_dir, "primary_model_selected_predictors.tsv"))
write_tsv_base(ratio_rule_table, file.path(out_dir, "xrf_ratio_small_value_rules.tsv"))
write_tsv_base(autocorr_summary, file.path(out_dir, "autocorrelation_robust_model_summary.tsv"))
write_tsv_base(autocorr_fit_metrics, file.path(out_dir, "autocorrelation_model_fit_metrics.tsv"))
write_tsv_base(autocorr_diagnostics, file.path(out_dir, "residual_autocorrelation_diagnostics.tsv"))
write_tsv_base(geochem_stability, file.path(out_dir, "geochemical_axis_stability_summary.tsv"))

if (!is.null(primary_fit$coefficients)) {
  write_tsv_base(primary_fit$coefficients, file.path(out_dir, "primary_model_summary_table.tsv"))
  write_tsv_base(primary_fit$fit, file.path(out_dir, "primary_model_fit_metrics.tsv"))
  write_tsv_base(primary_fitted_table, file.path(out_dir, "primary_model_fitted_vs_observed.tsv"))
}

if (!is.null(sens_age$coefficients)) {
  write_tsv_base(sens_age$coefficients, file.path(out_dir, "sensitivity_model_age_trend_summary.tsv"))
  write_tsv_base(sens_age$fit, file.path(out_dir, "sensitivity_model_age_trend_fit_metrics.tsv"))
}
if (!is.null(sens_mean$coefficients)) {
  write_tsv_base(sens_mean$coefficients, file.path(out_dir, "sensitivity_model_mean_aggregation_summary.tsv"))
  write_tsv_base(sens_mean$fit, file.path(out_dir, "sensitivity_model_mean_aggregation_fit_metrics.tsv"))
}
if (!is.null(sens_sparse) && !is.null(sens_sparse$coefficients)) {
  write_tsv_base(sens_sparse$coefficients, file.path(out_dir, "sensitivity_model_sparse_predictors_summary.tsv"))
  write_tsv_base(sens_sparse$fit, file.path(out_dir, "sensitivity_model_sparse_predictors_fit_metrics.tsv"))
}
write_tsv_base(sensitivity_overview, file.path(out_dir, "sensitivity_model_overview.tsv"))

if (nrow(valley_tbl) > 0) {
  write_tsv_base(valley_tbl[valley_tbl$valley_candidate, , drop = FALSE], file.path(out_dir, "valley_layers_table.tsv"))
  write_tsv_base(valley_tbl, file.path(out_dir, "valley_all_layers_with_residuals.tsv"))
  write_tsv_base(valley_compare, file.path(out_dir, "valley_vs_nonvalley_descriptive_comparison.tsv"))
}

caveats <- c(
  "This analysis is ST13 single-core and observational; associations are not causal proof.",
  "Depth/age structure is a major confounder and was explicitly modeled via spline terms.",
  "Autocorrelation-robust inference was attempted using nlme::gls with corCAR1(~ depth_in_core_cm); fallback uncertainty uses LM+Newey-West HAC if available, otherwise moving block bootstrap by depth order.",
  "Interpretive guardrail: a robust rb_sr coefficient does not by itself prove an Rb/Sr-specific mechanism if other correlated terrigenous/mineralogical representations show similar direction and magnitude.",
  "Interpretive guardrail: if rb_sr is stable while alternatives diverge materially, that pattern is more consistent with a relatively specific Rb/Sr-linked signal (still observational, not causal proof).",
  "Downstream sequencing-output metrics (initial, derep, avg_leng_initial, avg_len_derep) were excluded from explanatory-variable screening, co-variation, and model selection.",
  "Primary model uses XRF predictors aggregated to one value per ST13 library layer.",
  "Sparse foram/SST predictors were treated as sensitivity-only due to limited depth coverage.",
  "Valley analysis is secondary/supporting and based on residual definitions, not event attribution."
)
writeLines(caveats, con = file.path(out_dir, "analysis_caveats.txt"))

message("ST13 library concentration analysis complete.")
message("Outputs written to: ", out_dir)
