#!/usr/bin/env Rscript

# ST13 state-dependent terrigenous/mineralogical association analysis.
# Guardrails:
# - observational associations only (not causal proof)
# - interaction significance is interpreted as effect modification in this model framework
# - additive non-support of modifiers does not preclude interaction testing

options(stringsAsFactors = FALSE)

suppressWarnings({
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required.")
  }
})

# -------------------------------
# Paths and parameters
# -------------------------------
out_dir <- file.path("analysis", "st13_state_dependent_terrigenous_model")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

axis_table_path <- file.path("analysis", "st13_proxy_axes_model", "st13_axis_scores_table.tsv")

required_files <- c(axis_table_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files (run from repo root): ", paste(missing_files, collapse = ", "))
}

params <- list(
  depth_spline_df = 4,
  block_bootstrap_n = 600,
  block_bootstrap_seed = 13,
  stratified_min_n = 20,
  split_rule = "median_split_with_low_<=_median",
  interaction_plot_modifier_quantiles = "25th_vs_75th",
  best_model_rule = "lowest_AIC_among_stable_interaction_models"
)

# -------------------------------
# Helpers
# -------------------------------
write_tsv <- function(df, file) {
  utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

safe_spearman <- function(x, y, min_n = 10) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  if (length(x2) < min_n || length(unique(x2)) < 3 || length(unique(y2)) < 3) {
    return(list(n = length(x2), rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(stats::cor.test(x2, y2, method = "spearman", exact = FALSE))
  list(n = length(x2), rho = unname(ct$estimate), p = ct$p.value)
}

depth_residuals <- function(v, depth, df = 4) {
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v) & is.finite(depth)
  if (sum(ok) < 10 || length(unique(depth[ok])) < 6) return(out)
  fit <- tryCatch(stats::lm(v[ok] ~ splines::ns(depth[ok], df = df)), error = function(e) NULL)
  if (is.null(fit)) return(out)
  out[ok] <- stats::residuals(fit)
  out
}

extract_model_metrics <- function(mod, model_id) {
  y <- stats::model.response(stats::model.frame(mod))
  f <- stats::fitted(mod)
  r2 <- suppressWarnings(stats::cor(y, f, use = "complete.obs")^2)
  data.frame(
    model_id = model_id,
    n_obs = stats::nobs(mod),
    aic = stats::AIC(mod),
    bic = stats::BIC(mod),
    adj_r2 = summary(mod)$adj.r.squared,
    corr_sq_observed_vs_fitted = r2,
    stringsAsFactors = FALSE
  )
}

model_coef_table <- function(mod, model_id) {
  sm <- summary(mod)$coefficients
  ci <- tryCatch(stats::confint(mod), error = function(e) NULL)
  out <- data.frame(
    model_id = model_id,
    term = rownames(sm),
    estimate = sm[, 1],
    std_error = sm[, 2],
    stat = sm[, 3],
    p_value = sm[, 4],
    ci_low = if (!is.null(ci)) ci[rownames(sm), 1] else NA_real_,
    ci_high = if (!is.null(ci)) ci[rownames(sm), 2] else NA_real_,
    stringsAsFactors = FALSE
  )
  out
}

block_bootstrap_coefficients <- function(lm_model, model_data, n_boot = 600, block_len = NULL, seed = 13) {
  if (!inherits(lm_model, "lm")) return(NULL)
  dat <- model_data
  if (is.null(dat) || nrow(dat) < 20) return(NULL)

  n <- nrow(dat)
  if (is.null(block_len)) block_len <- max(3, floor(sqrt(n)))
  n_blocks <- ceiling(n / block_len)

  coef_names <- names(stats::coef(lm_model))
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(coef_names), dimnames = list(NULL, coef_names))
  fml <- stats::formula(lm_model)

  set.seed(seed)
  for (b in seq_len(n_boot)) {
    starts <- sample.int(n, size = n_blocks, replace = TRUE)
    idx <- integer(0)
    for (s in starts) idx <- c(idx, ((s - 1 + seq_len(block_len) - 1) %% n) + 1)
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
    bootstrap_se = NA_real_,
    bootstrap_ci_low = NA_real_,
    bootstrap_ci_high = NA_real_,
    bootstrap_p_two_sided = NA_real_,
    n_boot = n_boot,
    block_len = block_len,
    successful_replicates = sum(rowSums(is.finite(boot_mat)) > 0),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(coef_names)) {
    vv <- boot_mat[, coef_names[i]]
    vv <- vv[is.finite(vv)]
    if (length(vv) < max(100, floor(0.4 * n_boot))) next
    out$bootstrap_se[i] <- stats::sd(vv)
    ci <- as.numeric(stats::quantile(vv, probs = c(0.025, 0.975), na.rm = TRUE, type = 7))
    out$bootstrap_ci_low[i] <- ci[1]
    out$bootstrap_ci_high[i] <- ci[2]
    pboot <- 2 * min(mean(vv <= 0), mean(vv >= 0))
    out$bootstrap_p_two_sided[i] <- min(1, pboot)
  }
  out
}

compare_nested <- function(m0, m1, model_without, model_with, modifier_type) {
  cmp <- tryCatch(stats::anova(m0, m1), error = function(e) NULL)
  out <- data.frame(
    modifier_type = modifier_type,
    model_without_interaction = model_without,
    model_with_interaction = model_with,
    n_obs = stats::nobs(m1),
    aic_without = stats::AIC(m0),
    aic_with = stats::AIC(m1),
    delta_aic_with_minus_without = stats::AIC(m1) - stats::AIC(m0),
    df_change = NA_real_,
    f_stat = NA_real_,
    p_nested = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(cmp) && nrow(cmp) >= 2) {
    out$df_change <- cmp$Df[2]
    out$f_stat <- cmp$F[2]
    out$p_nested <- cmp$`Pr(>F)`[2]
  }
  out
}

stratify_axis_median <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  ifelse(x <= med, "low", "high")
}

stratified_effects <- function(df, modifier_col, modifier_label, min_n = 20, depth_df = 4) {
  stopifnot(modifier_col %in% names(df))
  strata <- stratify_axis_median(df[[modifier_col]])
  med <- stats::median(df[[modifier_col]], na.rm = TRUE)
  out <- data.frame(
    modifier_type = modifier_label,
    modifier_axis = modifier_col,
    split_rule = params$split_rule,
    median_cutpoint = med,
    stratum = c("low", "high"),
    n = NA_integer_,
    rho_raw = NA_real_,
    p_raw = NA_real_,
    rho_depth_adjusted = NA_real_,
    p_depth_adjusted = NA_real_,
    terrigenous_depth_adjusted_beta = NA_real_,
    terrigenous_depth_adjusted_se = NA_real_,
    terrigenous_depth_adjusted_p = NA_real_,
    terrigenous_depth_adjusted_ci_low = NA_real_,
    terrigenous_depth_adjusted_ci_high = NA_real_,
    stringsAsFactors = FALSE
  )

  for (s in c("low", "high")) {
    i <- which(out$stratum == s)
    sub <- df[strata == s, , drop = FALSE]
    ok <- stats::complete.cases(sub[, c("library_concentration_log10", "axis_terrigenous_mineralogical", "depth_in_core_cm")])
    sub <- sub[ok, , drop = FALSE]

    out$n[i] <- nrow(sub)
    if (nrow(sub) < 10) next

    raw <- safe_spearman(sub$axis_terrigenous_mineralogical, sub$library_concentration_log10, min_n = 10)
    out$rho_raw[i] <- raw$rho
    out$p_raw[i] <- raw$p

    y_res <- depth_residuals(sub$library_concentration_log10, sub$depth_in_core_cm, df = depth_df)
    x_res <- depth_residuals(sub$axis_terrigenous_mineralogical, sub$depth_in_core_cm, df = depth_df)
    adj <- safe_spearman(x_res, y_res, min_n = 10)
    out$rho_depth_adjusted[i] <- adj$rho
    out$p_depth_adjusted[i] <- adj$p

    if (nrow(sub) < min_n || stats::sd(sub$axis_terrigenous_mineralogical, na.rm = TRUE) == 0) next
    sub$z_axis_terrigenous_mineralogical <- as.numeric(scale(sub$axis_terrigenous_mineralogical))
    fit <- tryCatch(
      stats::lm(library_concentration_log10 ~ splines::ns(depth_in_core_cm, df = depth_df) + z_axis_terrigenous_mineralogical, data = sub),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    co <- summary(fit)$coefficients
    trm <- "z_axis_terrigenous_mineralogical"
    if (!(trm %in% rownames(co))) next
    ci <- tryCatch(stats::confint(fit), error = function(e) NULL)
    out$terrigenous_depth_adjusted_beta[i] <- co[trm, 1]
    out$terrigenous_depth_adjusted_se[i] <- co[trm, 2]
    out$terrigenous_depth_adjusted_p[i] <- co[trm, 4]
    out$terrigenous_depth_adjusted_ci_low[i] <- if (!is.null(ci)) ci[trm, 1] else NA_real_
    out$terrigenous_depth_adjusted_ci_high[i] <- if (!is.null(ci)) ci[trm, 2] else NA_real_
  }

  out
}

plot_interaction_effect <- function(model, model_df, modifier_var, plot_file, plot_title, modifier_label) {
  q <- stats::quantile(model_df[[modifier_var]], probs = c(0.25, 0.75), na.rm = TRUE)
  terrig_seq <- seq(min(model_df$z_axis_terrigenous_mineralogical, na.rm = TRUE),
                    max(model_df$z_axis_terrigenous_mineralogical, na.rm = TRUE),
                    length.out = 120)
  depth_fix <- stats::median(model_df$depth_in_core_cm, na.rm = TRUE)

  nd_low <- data.frame(
    depth_in_core_cm = depth_fix,
    z_axis_terrigenous_mineralogical = terrig_seq,
    z_axis_productivity_carbonate = if ("z_axis_productivity_carbonate" %in% names(model_df)) stats::median(model_df$z_axis_productivity_carbonate, na.rm = TRUE) else 0,
    z_axis_redox = if ("z_axis_redox" %in% names(model_df)) stats::median(model_df$z_axis_redox, na.rm = TRUE) else 0
  )
  nd_high <- nd_low
  nd_low[[modifier_var]] <- q[1]
  nd_high[[modifier_var]] <- q[2]

  p_low <- stats::predict(model, newdata = nd_low, se.fit = TRUE)
  p_high <- stats::predict(model, newdata = nd_high, se.fit = TRUE)

  png(plot_file, width = 1200, height = 900, res = 150)
  y_all <- c(p_low$fit - 1.96 * p_low$se.fit, p_low$fit + 1.96 * p_low$se.fit,
             p_high$fit - 1.96 * p_high$se.fit, p_high$fit + 1.96 * p_high$se.fit)
  plot(
    terrig_seq, p_low$fit,
    type = "l", lwd = 2, col = "#1f77b4",
    ylim = range(y_all, finite = TRUE),
    xlab = "Terrigenous/mineralogical axis (standardized)",
    ylab = "Predicted log10(library concentration)",
    main = plot_title
  )
  lines(terrig_seq, p_high$fit, lwd = 2, col = "#d62728")
  lines(terrig_seq, p_low$fit - 1.96 * p_low$se.fit, lty = 2, col = "#1f77b4")
  lines(terrig_seq, p_low$fit + 1.96 * p_low$se.fit, lty = 2, col = "#1f77b4")
  lines(terrig_seq, p_high$fit - 1.96 * p_high$se.fit, lty = 2, col = "#d62728")
  lines(terrig_seq, p_high$fit + 1.96 * p_high$se.fit, lty = 2, col = "#d62728")
  legend(
    "topright",
    legend = c(
      paste0(modifier_label, " low (25th pct)"),
      paste0(modifier_label, " high (75th pct)")
    ),
    col = c("#1f77b4", "#d62728"),
    lty = 1,
    lwd = 2,
    bty = "n"
  )
  abline(v = 0, lty = 3, col = "#888888")
  dev.off()
}

plot_stratified_summary <- function(strat_df, plot_file) {
  d <- strat_df
  d$label <- paste0(d$modifier_type, "_", d$stratum)

  png(plot_file, width = 1200, height = 900, res = 150)
  op <- par(mar = c(5, 10, 4, 2))
  y <- seq_len(nrow(d))
  x <- d$terrigenous_depth_adjusted_beta
  xl <- d$terrigenous_depth_adjusted_ci_low
  xh <- d$terrigenous_depth_adjusted_ci_high

  xlim <- range(c(xl, xh, 0), na.rm = TRUE)
  plot(
    x, y,
    xlim = xlim,
    yaxt = "n",
    pch = 16,
    col = "#1f77b4",
    xlab = "Depth-adjusted terrigenous coefficient (within stratum)",
    ylab = "",
    main = "Stratified terrigenous effects by sediment state"
  )
  axis(2, at = y, labels = paste0(d$label, " (n=", d$n, ")"), las = 1)
  ok <- is.finite(xl) & is.finite(xh)
  if (any(ok)) segments(xl[ok], y[ok], xh[ok], y[ok], lwd = 2, col = "#1f77b4")
  abline(v = 0, lty = 2, col = "#d62728")
  par(op)
  dev.off()
}

plot_observed_fitted <- function(model, model_df, plot_file, plot_title) {
  obs <- model_df$library_concentration_log10
  fit <- stats::fitted(model)
  ord <- order(model_df$depth_in_core_cm)

  png(plot_file, width = 1200, height = 900, res = 150)
  op <- par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
  plot(
    obs, fit,
    pch = 16, col = "#1f77b4",
    xlab = "Observed log10(library concentration)",
    ylab = "Fitted",
    main = paste0(plot_title, ": observed vs fitted")
  )
  abline(0, 1, lty = 2, col = "#d62728")

  plot(
    model_df$depth_in_core_cm[ord],
    obs[ord],
    type = "p", pch = 16, col = "#2ca02c",
    xlab = "Depth in core (cm)",
    ylab = "log10(library concentration)",
    main = paste0(plot_title, ": observed and fitted by depth")
  )
  lines(model_df$depth_in_core_cm[ord], fit[ord], col = "#1f77b4", lwd = 2)
  legend("topright", legend = c("Observed", "Fitted"), col = c("#2ca02c", "#1f77b4"), pch = c(16, NA), lty = c(NA, 1), bty = "n")
  par(op)
  dev.off()
}

# -------------------------------
# Read and validate axis table
# -------------------------------
axis_tbl <- read.delim(axis_table_path, sep = "\t", check.names = FALSE)

required_cols <- c(
  "depth_in_core_cm",
  "library_concentration",
  "library_concentration_log10",
  "axis_terrigenous_mineralogical",
  "axis_productivity_carbonate",
  "axis_redox"
)
missing_cols <- required_cols[!(required_cols %in% names(axis_tbl))]
if (length(missing_cols) > 0) {
  stop("Axis table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Complete-case table for modeling with all three axes available
model_df <- axis_tbl[, required_cols, drop = FALSE]
model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
if (nrow(model_df) < 30) {
  stop("Too few complete rows for interaction modeling.")
}

# Standardize interaction-participating axes
model_df$z_axis_terrigenous_mineralogical <- as.numeric(scale(model_df$axis_terrigenous_mineralogical))
model_df$z_axis_productivity_carbonate <- as.numeric(scale(model_df$axis_productivity_carbonate))
model_df$z_axis_redox <- as.numeric(scale(model_df$axis_redox))

if (any(!is.finite(model_df$z_axis_terrigenous_mineralogical)) ||
    any(!is.finite(model_df$z_axis_productivity_carbonate)) ||
    any(!is.finite(model_df$z_axis_redox))) {
  stop("Standardization produced non-finite values; check axis variance and missingness.")
}

# -------------------------------
# Interaction models and comparisons
# -------------------------------
f_depth <- library_concentration_log10 ~ splines::ns(depth_in_core_cm, df = params$depth_spline_df)

f_prod_add <- update(f_depth, . ~ . + z_axis_terrigenous_mineralogical + z_axis_productivity_carbonate)
f_prod_int <- update(f_prod_add, . ~ . + z_axis_terrigenous_mineralogical:z_axis_productivity_carbonate)

f_redox_add <- update(f_depth, . ~ . + z_axis_terrigenous_mineralogical + z_axis_redox)
f_redox_int <- update(f_redox_add, . ~ . + z_axis_terrigenous_mineralogical:z_axis_redox)

f_combined_int <- update(
  f_depth,
  . ~ . + z_axis_terrigenous_mineralogical + z_axis_productivity_carbonate + z_axis_redox +
    z_axis_terrigenous_mineralogical:z_axis_productivity_carbonate +
    z_axis_terrigenous_mineralogical:z_axis_redox
)

m_prod_add <- stats::lm(f_prod_add, data = model_df)
m_prod_int <- stats::lm(f_prod_int, data = model_df)
m_redox_add <- stats::lm(f_redox_add, data = model_df)
m_redox_int <- stats::lm(f_redox_int, data = model_df)

m_combined_int <- tryCatch(stats::lm(f_combined_int, data = model_df), error = function(e) NULL)
combined_stable <- FALSE
combined_note <- "combined model fit attempt"
if (!is.null(m_combined_int)) {
  co <- summary(m_combined_int)$coefficients
  kappa_x <- tryCatch(kappa(model.matrix(m_combined_int), exact = FALSE), error = function(e) Inf)
  if (all(is.finite(co[, 2])) && kappa_x < 1000) {
    combined_stable <- TRUE
    combined_note <- paste0("combined model retained (kappa=", signif(kappa_x, 4), ")")
  } else {
    combined_note <- paste0("combined model not considered stable/interpretable (kappa=", signif(kappa_x, 4), ")")
  }
}

model_metrics <- rbind(
  extract_model_metrics(m_prod_add, "prod_additive"),
  extract_model_metrics(m_prod_int, "prod_interaction"),
  extract_model_metrics(m_redox_add, "redox_additive"),
  extract_model_metrics(m_redox_int, "redox_interaction")
)
if (combined_stable) model_metrics <- rbind(model_metrics, extract_model_metrics(m_combined_int, "combined_two_interactions"))

coef_all <- rbind(
  model_coef_table(m_prod_add, "prod_additive"),
  model_coef_table(m_prod_int, "prod_interaction"),
  model_coef_table(m_redox_add, "redox_additive"),
  model_coef_table(m_redox_int, "redox_interaction")
)
if (combined_stable) coef_all <- rbind(coef_all, model_coef_table(m_combined_int, "combined_two_interactions"))

comparison_tbl <- rbind(
  compare_nested(m_prod_add, m_prod_int, "prod_additive", "prod_interaction", "productivity_carbonate"),
  compare_nested(m_redox_add, m_redox_int, "redox_additive", "redox_interaction", "redox")
)

interaction_terms <- c(
  prod_interaction = "z_axis_terrigenous_mineralogical:z_axis_productivity_carbonate",
  redox_interaction = "z_axis_terrigenous_mineralogical:z_axis_redox",
  combined_prod_interaction = "z_axis_terrigenous_mineralogical:z_axis_productivity_carbonate",
  combined_redox_interaction = "z_axis_terrigenous_mineralogical:z_axis_redox"
)

interaction_summary <- coef_all[coef_all$term %in% unique(interaction_terms), , drop = FALSE]
interaction_summary$modifier_type <- ifelse(
  grepl("productivity", interaction_summary$term, fixed = TRUE),
  "productivity_carbonate", "redox"
)
interaction_summary$interaction_supported_p_lt_0_05 <- interaction_summary$p_value < 0.05

# -------------------------------
# Stratified summaries
# -------------------------------
strat_prod <- stratified_effects(
  model_df,
  modifier_col = "axis_productivity_carbonate",
  modifier_label = "productivity_carbonate",
  min_n = params$stratified_min_n,
  depth_df = params$depth_spline_df
)
strat_redox <- stratified_effects(
  model_df,
  modifier_col = "axis_redox",
  modifier_label = "redox",
  min_n = params$stratified_min_n,
  depth_df = params$depth_spline_df
)
strat_tbl <- rbind(strat_prod, strat_redox)

# -------------------------------
# Block bootstrap for interaction coefficients
# -------------------------------
boot_rows <- list()

boot_prod <- block_bootstrap_coefficients(
  lm_model = m_prod_int,
  model_data = model_df,
  n_boot = params$block_bootstrap_n,
  seed = params$block_bootstrap_seed
)
if (!is.null(boot_prod)) {
  boot_prod$model_id <- "prod_interaction"
  boot_rows[[length(boot_rows) + 1]] <- boot_prod
}

boot_redox <- block_bootstrap_coefficients(
  lm_model = m_redox_int,
  model_data = model_df,
  n_boot = params$block_bootstrap_n,
  seed = params$block_bootstrap_seed
)
if (!is.null(boot_redox)) {
  boot_redox$model_id <- "redox_interaction"
  boot_rows[[length(boot_rows) + 1]] <- boot_redox
}

if (combined_stable) {
  boot_combined <- block_bootstrap_coefficients(
    lm_model = m_combined_int,
    model_data = model_df,
    n_boot = params$block_bootstrap_n,
    seed = params$block_bootstrap_seed
  )
  if (!is.null(boot_combined)) {
    boot_combined$model_id <- "combined_two_interactions"
    boot_rows[[length(boot_rows) + 1]] <- boot_combined
  }
}

if (length(boot_rows) > 0) {
  boot_tbl <- do.call(rbind, boot_rows)
  boot_tbl$interaction_term <- boot_tbl$term %in% unique(interaction_terms)
} else {
  boot_tbl <- data.frame(
    model_id = character(0),
    term = character(0),
    estimate = numeric(0),
    bootstrap_se = numeric(0),
    bootstrap_ci_low = numeric(0),
    bootstrap_ci_high = numeric(0),
    bootstrap_p_two_sided = numeric(0),
    n_boot = integer(0),
    block_len = integer(0),
    successful_replicates = integer(0),
    interaction_term = logical(0),
    stringsAsFactors = FALSE
  )
}

# -------------------------------
# Select best interaction model for diagnostics
# -------------------------------
candidate_models <- list(
  prod_interaction = m_prod_int,
  redox_interaction = m_redox_int
)
if (combined_stable) candidate_models$combined_two_interactions <- m_combined_int

candidate_aic <- vapply(candidate_models, stats::AIC, numeric(1))
best_model_id <- names(which.min(candidate_aic))[1]
best_model <- candidate_models[[best_model_id]]

# -------------------------------
# Plots
# -------------------------------
plot_interaction_effect(
  model = m_prod_int,
  model_df = model_df,
  modifier_var = "z_axis_productivity_carbonate",
  plot_file = file.path(plot_dir, "interaction_effect_terrigenous_x_productivity_carbonate.png"),
  plot_title = "Terrigenous × productivity/carbonate interaction model",
  modifier_label = "Productivity/carbonate state"
)

plot_interaction_effect(
  model = m_redox_int,
  model_df = model_df,
  modifier_var = "z_axis_redox",
  plot_file = file.path(plot_dir, "interaction_effect_terrigenous_x_redox.png"),
  plot_title = "Terrigenous × redox interaction model",
  modifier_label = "Redox state"
)

plot_stratified_summary(
  strat_df = strat_tbl,
  plot_file = file.path(plot_dir, "stratified_terrigenous_effect_summary.png")
)

plot_observed_fitted(
  model = best_model,
  model_df = model_df,
  plot_file = file.path(plot_dir, "best_interaction_model_observed_fitted_diagnostic.png"),
  plot_title = best_model_id
)

# -------------------------------
# Caveats
# -------------------------------
prod_cmp <- comparison_tbl[comparison_tbl$modifier_type == "productivity_carbonate", , drop = FALSE]
redox_cmp <- comparison_tbl[comparison_tbl$modifier_type == "redox", , drop = FALSE]

prod_support <- nrow(prod_cmp) == 1 && is.finite(prod_cmp$p_nested) && prod_cmp$p_nested < 0.05
redox_support <- nrow(redox_cmp) == 1 && is.finite(redox_cmp$p_nested) && redox_cmp$p_nested < 0.05

caveats <- c(
  "Interaction terms test whether the terrigenous/mineralogical association differs by sediment state within this observational depth-adjusted framework.",
  "A significant interaction indicates state-dependent slope differences; it does not imply the modifier alone drives DNA preservation.",
  "No interaction support indicates the terrigenous association appears relatively stable across that state dimension in this framework.",
  "Prior additive non-support of productivity/carbonate does not preclude interaction effects; this analysis explicitly tests effect modification.",
  ifelse(prod_support,
         "Productivity/carbonate interaction received statistical support in nested comparison.",
         "Productivity/carbonate interaction did not receive clear support in nested comparison."),
  ifelse(redox_support,
         "Redox interaction received statistical support in nested comparison.",
         "Redox interaction did not receive clear support in nested comparison."),
  combined_note,
  "Median split stratification is an interpretability aid and may reduce power relative to continuous interaction models."
)

# -------------------------------
# Write outputs
# -------------------------------
write_tsv(data.frame(parameter = names(params), value = unlist(params), stringsAsFactors = FALSE), file.path(out_dir, "analysis_parameters.tsv"))
write_tsv(interaction_summary, file.path(out_dir, "interaction_model_summary.tsv"))
write_tsv(comparison_tbl, file.path(out_dir, "interaction_model_comparison.tsv"))
write_tsv(strat_tbl, file.path(out_dir, "state_stratified_terrigenous_effects.tsv"))
write_tsv(boot_tbl, file.path(out_dir, "interaction_bootstrap_summary.tsv"))
write_tsv(model_metrics, file.path(out_dir, "interaction_model_fit_metrics.tsv"))
write_tsv(coef_all, file.path(out_dir, "interaction_model_coefficients_all.tsv"))
writeLines(caveats, con = file.path(out_dir, "analysis_caveats.txt"))

message("ST13 state-dependent terrigenous analysis complete.")
message("Outputs written to: ", out_dir)
message("Best interaction model by AIC: ", best_model_id)
