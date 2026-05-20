#!/usr/bin/env Rscript

# ST13 grouped-proxy / latent-axis analysis for library concentration.
# Guardrails:
# - observational associations only (not causal proof)
# - no downstream sequencing-output metrics used as explanatory variables
# - axis construction is conservative and explicitly documented

options(stringsAsFactors = FALSE)

suppressWarnings({
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required.")
  }
})

# -------------------------------
# Paths and parameters
# -------------------------------
out_dir <- file.path("analysis", "st13_proxy_axes_model")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path("data", "metadata_v5.tsv")
xrf_path <- file.path("data", "combined_xrf_geochemistry_curated.csv")

required_files <- c(metadata_path, xrf_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files (run from repo root): ", paste(missing_files, collapse = ", "))
}

params <- list(
  st13_core_label = "ST13",
  response_pseudocount_rule = "min_positive_library_concentration / 2",
  xrf_layer_aggregation = "median",
  depth_spline_df = 4,
  screening_min_n = 10,
  axis_min_complete_rows = 20,
  axis_min_unique_values = 4,
  block_bootstrap_n = 600,
  block_bootstrap_seed = 13
)

# -------------------------------
# Helpers
# -------------------------------
write_tsv <- function(df, file) {
  utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

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
  d <- sort(unique(depths))
  mids <- (d[-1] + d[-length(d)]) / 2
  data.frame(depth_in_core_cm = d, lower = c(-Inf, mids), upper = c(mids, Inf))
}

aggregate_xrf_to_layers <- function(xrf_df, intervals_df, numeric_cols) {
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
      if (length(vals) > 0) out[[nm]][i] <- median(vals)
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

build_axis <- function(df, axis_name, candidates, orientation_priority, min_complete = 20, min_unique = 4) {
  candidates <- candidates[candidates %in% names(df)]
  membership <- data.frame(
    axis = axis_name,
    variable = candidates,
    included = FALSE,
    exclusion_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  if (length(candidates) == 0) {
    return(list(
      score = rep(NA_real_, nrow(df)),
      method = "none",
      orientation_variable = NA_character_,
      orientation_rho_before = NA_real_,
      membership = membership,
      loadings = data.frame(),
      notes = "No candidate variables available in table."
    ))
  }

  usable <- character(0)
  for (v in candidates) {
    x <- df[[v]]
    n_nonmiss <- sum(is.finite(x))
    n_unique <- length(unique(x[is.finite(x)]))
    if (n_nonmiss < min_complete) {
      membership$exclusion_reason[membership$variable == v] <- paste0("non-missing < ", min_complete)
      next
    }
    if (n_unique < min_unique) {
      membership$exclusion_reason[membership$variable == v] <- paste0("unique values < ", min_unique)
      next
    }
    usable <- c(usable, v)
    membership$included[membership$variable == v] <- TRUE
    membership$exclusion_reason[membership$variable == v] <- "included"
  }

  if (length(usable) == 0) {
    return(list(
      score = rep(NA_real_, nrow(df)),
      method = "none",
      orientation_variable = NA_character_,
      orientation_rho_before = NA_real_,
      membership = membership,
      loadings = data.frame(),
      notes = "No candidate variable passed minimum coverage/variability filters."
    ))
  }

  zmat <- as.data.frame(lapply(df[usable], function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }))
  names(zmat) <- usable

  method <- "z_mean"
  axis_score <- rep(NA_real_, nrow(df))
  load_tbl <- data.frame(
    axis = axis_name,
    variable = usable,
    method = NA_character_,
    raw_loading_or_weight = NA_real_,
    oriented_loading_or_weight = NA_real_,
    stringsAsFactors = FALSE
  )

  complete_rows <- complete.cases(zmat)
  n_complete <- sum(complete_rows)

  if (length(usable) >= 2 && n_complete >= max(min_complete, length(usable) + 3)) {
    pca <- tryCatch(stats::prcomp(zmat[complete_rows, , drop = FALSE], center = FALSE, scale. = FALSE), error = function(e) NULL)
    if (!is.null(pca) && is.matrix(pca$rotation) && ncol(pca$rotation) >= 1) {
      method <- "pca_pc1"
      axis_score[complete_rows] <- pca$x[, 1]
      load_tbl$method <- method
      raw_load <- pca$rotation[, 1]
      load_tbl$raw_loading_or_weight <- as.numeric(raw_load[load_tbl$variable])
      load_tbl$oriented_loading_or_weight <- load_tbl$raw_loading_or_weight
    }
  }

  if (method == "z_mean") {
    axis_score <- rowMeans(zmat, na.rm = TRUE)
    axis_score[!is.finite(axis_score)] <- NA_real_
    load_tbl$method <- method
    load_tbl$raw_loading_or_weight <- 1
    load_tbl$oriented_loading_or_weight <- 1
  }

  orient_var <- NA_character_
  orient_rho <- NA_real_
  for (v in orientation_priority) {
    if (v %in% usable) {
      orient_var <- v
      break
    }
  }
  if (is.na(orient_var) && length(usable) > 0) orient_var <- usable[1]

  if (!is.na(orient_var)) {
    cor0 <- suppressWarnings(cor(axis_score, df[[orient_var]], method = "spearman", use = "pairwise.complete.obs"))
    orient_rho <- as.numeric(cor0)
    if (is.finite(orient_rho) && orient_rho < 0) {
      axis_score <- -1 * axis_score
      load_tbl$oriented_loading_or_weight <- -1 * load_tbl$oriented_loading_or_weight
    }
  }

  note <- paste0("Axis built with method=", method, "; usable vars=", paste(usable, collapse = ", "),
                 "; orientation variable=", orient_var)

  list(
    score = axis_score,
    method = method,
    orientation_variable = orient_var,
    orientation_rho_before = orient_rho,
    membership = membership,
    loadings = load_tbl,
    notes = note
  )
}

extract_model_metrics <- function(mod, model_id, extra = list()) {
  y <- model.response(model.frame(mod))
  f <- fitted(mod)
  r2 <- suppressWarnings(cor(y, f, use = "complete.obs")^2)
  out <- data.frame(
    model_id = model_id,
    n_obs = stats::nobs(mod),
    aic = AIC(mod),
    bic = BIC(mod),
    adj_r2 = summary(mod)$adj.r.squared,
    corr_sq_observed_vs_fitted = r2,
    stringsAsFactors = FALSE
  )
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}

model_coef_table <- function(mod, model_id) {
  sm <- summary(mod)$coefficients
  ci <- tryCatch(confint(mod), error = function(e) NULL)
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

  coef_names <- names(coef(lm_model))
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(coef_names), dimnames = list(NULL, coef_names))
  fml <- formula(lm_model)

  set.seed(seed)
  for (b in seq_len(n_boot)) {
    starts <- sample.int(n, size = n_blocks, replace = TRUE)
    idx <- integer(0)
    for (s in starts) idx <- c(idx, ((s - 1 + seq_len(block_len) - 1) %% n) + 1)
    idx <- idx[seq_len(n)]
    fit_b <- tryCatch(lm(fml, data = dat[idx, , drop = FALSE]), error = function(e) NULL)
    if (is.null(fit_b)) next
    cb <- coef(fit_b)
    boot_mat[b, names(cb)] <- cb
  }

  est <- coef(lm_model)
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
    out$bootstrap_se[i] <- sd(vv)
    ci <- as.numeric(quantile(vv, probs = c(0.025, 0.975), na.rm = TRUE, type = 7))
    out$bootstrap_ci_low[i] <- ci[1]
    out$bootstrap_ci_high[i] <- ci[2]
    pboot <- 2 * min(mean(vv <= 0), mean(vv >= 0))
    out$bootstrap_p_two_sided[i] <- min(1, pboot)
  }
  out
}

# -------------------------------
# Read data and build ST13 layer-level table
# -------------------------------
metadata <- read.delim(metadata_path, sep = "\t", check.names = FALSE)
xrf <- read.csv(xrf_path, check.names = FALSE)

meta_st13 <- metadata[metadata$core == params$st13_core_label, , drop = FALSE]
if (nrow(meta_st13) == 0) stop("No ST13 rows found in metadata_v5.tsv")

collapse_numeric <- function(v) {
  vv <- as.numeric(v)
  vv <- vv[is.finite(vv)]
  if (length(vv) == 0) return(NA_real_)
  mean(vv)
}

layer_key <- interaction(meta_st13$depth_in_core_cm, meta_st13$y_bp, drop = TRUE)
idx_split <- split(seq_len(nrow(meta_st13)), layer_key)

layers <- data.frame(
  layer_id = sprintf("ST13_layer_%03d", seq_along(idx_split)),
  depth_in_core_cm = vapply(idx_split, function(ix) collapse_numeric(meta_st13$depth_in_core_cm[ix]), numeric(1)),
  y_bp = vapply(idx_split, function(ix) collapse_numeric(meta_st13$y_bp[ix]), numeric(1)),
  library_concentration = vapply(idx_split, function(ix) collapse_numeric(meta_st13$library_concentration[ix]), numeric(1)),
  stringsAsFactors = FALSE
)
layers <- layers[order(layers$depth_in_core_cm), ]
rownames(layers) <- NULL

min_pos <- safe_min_positive(layers$library_concentration, default = 1e-6)
pseudocount <- min_pos / 2
layers$library_concentration_log10 <- log10(layers$library_concentration + pseudocount)

xrf_st13 <- xrf[xrf$core == params$st13_core_label, , drop = FALSE]
if (nrow(xrf_st13) == 0) stop("No ST13 rows found in combined_xrf_geochemistry_curated.csv")

numeric_xrf <- names(xrf_st13)[vapply(xrf_st13, is.numeric, logical(1))]
numeric_xrf <- setdiff(numeric_xrf, c("depth_in_core_cm", "y_bp"))

intervals <- layer_intervals(layers$depth_in_core_cm)
xrf_med <- aggregate_xrf_to_layers(xrf_st13, intervals, numeric_cols = numeric_xrf)

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
ratio_rule_table <- do.call(rbind, ratio_rules)

layer_tbl <- merge(layers, xrf_med, by = "depth_in_core_cm", all.x = TRUE, sort = FALSE)
layer_tbl <- layer_tbl[order(layer_tbl$depth_in_core_cm), ]
rownames(layer_tbl) <- NULL

# -------------------------------
# Axis candidate groups
# -------------------------------
axis_groups <- list(
  productivity_carbonate = c("ba_ti", "p_ti", "br_ti", "si_ti", "ca_ti", "ca"),
  terrigenous_mineralogical = c("rb_sr", "al_si", "ti_al", "si_ti", "rb", "zr"),
  redox = c("mn_fe", "fe_al", "mn", "fe")
)

orientation_priority <- list(
  productivity_carbonate = c("ca_ti", "ba_ti", "p_ti", "br_ti", "si_ti", "ca"),
  terrigenous_mineralogical = c("rb_sr", "al_si", "ti_al", "rb", "zr", "si_ti"),
  redox = c("mn_fe", "fe_al", "mn", "fe")
)

axis_results <- list()
membership_all <- data.frame()
loadings_all <- data.frame()
axis_notes <- c()

for (ax in names(axis_groups)) {
  res <- build_axis(
    df = layer_tbl,
    axis_name = ax,
    candidates = axis_groups[[ax]],
    orientation_priority = orientation_priority[[ax]],
    min_complete = params$axis_min_complete_rows,
    min_unique = params$axis_min_unique_values
  )
  axis_results[[ax]] <- res
  layer_tbl[[paste0("axis_", ax)]] <- res$score
  membership_all <- rbind(membership_all, res$membership)
  if (nrow(res$loadings) > 0) loadings_all <- rbind(loadings_all, res$loadings)
  axis_notes <- c(axis_notes, paste0(ax, ": ", res$notes))
}

axis_method_summary <- data.frame(
  axis = names(axis_results),
  method = vapply(axis_results, function(x) x$method, character(1)),
  orientation_variable = vapply(axis_results, function(x) x$orientation_variable, character(1)),
  orientation_rho_before = vapply(axis_results, function(x) x$orientation_rho_before, numeric(1)),
  stringsAsFactors = FALSE
)

# -------------------------------
# Axis screening: raw + depth-detrended Spearman
# -------------------------------
response <- layer_tbl$library_concentration_log10
depth <- layer_tbl$depth_in_core_cm
res_response <- depth_residuals(response, depth, df = params$depth_spline_df)

axis_vars <- paste0("axis_", names(axis_groups))
screen <- data.frame(
  axis = names(axis_groups),
  axis_variable = axis_vars,
  n_raw = NA_integer_,
  rho_raw = NA_real_,
  p_raw = NA_real_,
  n_depth_detrended = NA_integer_,
  rho_depth_detrended = NA_real_,
  p_depth_detrended = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(screen))) {
  v <- screen$axis_variable[i]
  x <- layer_tbl[[v]]
  raw <- safe_spearman(x, response, min_n = params$screening_min_n)
  screen$n_raw[i] <- raw$n
  screen$rho_raw[i] <- raw$rho
  screen$p_raw[i] <- raw$p

  x_res <- depth_residuals(x, depth, df = params$depth_spline_df)
  det <- safe_spearman(x_res, res_response, min_n = params$screening_min_n)
  screen$n_depth_detrended[i] <- det$n
  screen$rho_depth_detrended[i] <- det$rho
  screen$p_depth_detrended[i] <- det$p
}

# -------------------------------
# Primary grouped-axis model and nested comparison
# -------------------------------
model_df <- layer_tbl[, c(
  "depth_in_core_cm",
  "y_bp",
  "library_concentration",
  "library_concentration_log10",
  "axis_terrigenous_mineralogical",
  "axis_redox",
  "axis_productivity_carbonate"
), drop = FALSE]

model_df <- model_df[complete.cases(model_df[, c(
  "depth_in_core_cm",
  "library_concentration_log10",
  "axis_terrigenous_mineralogical",
  "axis_redox",
  "axis_productivity_carbonate"
)]), , drop = FALSE]

if (nrow(model_df) < 25) {
  stop("Too few complete rows for grouped-axis nested modeling after axis construction.")
}

z_axis <- c("axis_terrigenous_mineralogical", "axis_redox", "axis_productivity_carbonate")
for (v in z_axis) {
  s <- sd(model_df[[v]], na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop("Axis has zero/invalid variance in model_df: ", v)
  model_df[[paste0("z_", v)]] <- (model_df[[v]] - mean(model_df[[v]], na.rm = TRUE)) / s
}

f_depth <- library_concentration_log10 ~ splines::ns(depth_in_core_cm, df = params$depth_spline_df)
f_t <- update(f_depth, . ~ . + z_axis_terrigenous_mineralogical)
f_tr <- update(f_t, . ~ . + z_axis_redox)
f_trp <- update(f_tr, . ~ . + z_axis_productivity_carbonate)

m_depth <- lm(f_depth, data = model_df)
m_t <- lm(f_t, data = model_df)
m_tr <- lm(f_tr, data = model_df)
m_trp <- lm(f_trp, data = model_df)

metrics <- rbind(
  extract_model_metrics(m_depth, "depth_only", list(step = 1, includes_productivity_axis = FALSE)),
  extract_model_metrics(m_t, "depth_plus_terrigenous", list(step = 2, includes_productivity_axis = FALSE)),
  extract_model_metrics(m_tr, "depth_plus_terrigenous_plus_redox", list(step = 3, includes_productivity_axis = FALSE)),
  extract_model_metrics(m_trp, "depth_plus_terrigenous_plus_redox_plus_productivity", list(step = 4, includes_productivity_axis = TRUE))
)

anova_cmp <- anova(m_depth, m_t, m_tr, m_trp)
nested_cmp <- data.frame(
  model_id = c(
    "depth_only",
    "depth_plus_terrigenous",
    "depth_plus_terrigenous_plus_redox",
    "depth_plus_terrigenous_plus_redox_plus_productivity"
  )[seq_len(nrow(anova_cmp))],
  residual_df = anova_cmp$Res.Df,
  rss = anova_cmp$RSS,
  df_change = c(NA, diff(anova_cmp$Res.Df) * -1),
  ss_change = c(NA, diff(anova_cmp$RSS) * -1),
  f_stat = anova_cmp$F,
  p_value = anova_cmp$`Pr(>F)`,
  stringsAsFactors = FALSE
)

final_coef <- model_coef_table(m_trp, model_id = "depth_plus_terrigenous_plus_redox_plus_productivity")

boot_coef <- block_bootstrap_coefficients(
  lm_model = m_trp,
  model_data = model_df,
  n_boot = params$block_bootstrap_n,
  seed = params$block_bootstrap_seed
)

if (is.null(boot_coef)) {
  boot_coef <- data.frame(
    term = names(coef(m_trp)),
    estimate = as.numeric(coef(m_trp)),
    bootstrap_se = NA_real_,
    bootstrap_ci_low = NA_real_,
    bootstrap_ci_high = NA_real_,
    bootstrap_p_two_sided = NA_real_,
    n_boot = params$block_bootstrap_n,
    block_len = NA_real_,
    successful_replicates = NA_real_,
    stringsAsFactors = FALSE
  )
}

grouped_summary <- merge(final_coef, boot_coef, by = c("term", "estimate"), all.x = TRUE, sort = FALSE)

# -------------------------------
# Tables for output
# -------------------------------
axis_scores_out <- model_df[, c(
  "depth_in_core_cm", "y_bp", "library_concentration", "library_concentration_log10",
  "axis_productivity_carbonate", "axis_terrigenous_mineralogical", "axis_redox",
  "z_axis_productivity_carbonate", "z_axis_terrigenous_mineralogical", "z_axis_redox"
), drop = FALSE]

axis_membership_out <- membership_all
axis_membership_out$group_candidates_defined <- vapply(
  axis_membership_out$axis,
  function(ax) paste(axis_groups[[ax]], collapse = ";"),
  character(1)
)

axis_loadings_out <- merge(loadings_all, axis_method_summary, by = "axis", all.x = TRUE, sort = FALSE)

# -------------------------------
# Plots
# -------------------------------
png(file.path(plot_dir, "axis_loadings_contribution_plot.png"), width = 1800, height = 900, res = 150)
if (nrow(loadings_all) > 0) {
  op <- par(mfrow = c(1, 3), mar = c(7, 4, 4, 1))
  for (ax in names(axis_groups)) {
    d <- loadings_all[loadings_all$axis == ax, , drop = FALSE]
    if (nrow(d) == 0) {
      plot.new(); title(main = ax); text(0.5, 0.5, "No loadings/weights")
      next
    }
    ord <- order(d$oriented_loading_or_weight)
    dd <- d[ord, , drop = FALSE]
    barplot(
      height = dd$oriented_loading_or_weight,
      names.arg = dd$variable,
      las = 2,
      col = "#1f77b4",
      ylab = "Oriented loading/weight",
      main = paste0(ax, " (", unique(dd$method), ")")
    )
    abline(h = 0, lty = 2)
  }
  par(op)
} else {
  plot.new(); text(0.5, 0.5, "No axis loadings available")
}
dev.off()

png(file.path(plot_dir, "library_concentration_vs_axes.png"), width = 1800, height = 650, res = 150)
op <- par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
for (ax in names(axis_groups)) {
  x <- layer_tbl[[paste0("axis_", ax)]]
  y <- layer_tbl$library_concentration_log10
  plot(x, y,
       pch = 16, col = "#2ca02c", cex = 0.8,
       xlab = paste0("axis_", ax),
       ylab = "log10(library concentration + pseudocount)",
       main = paste0("Library vs ", ax, " axis"))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) >= 10) abline(lm(y[ok] ~ x[ok]), col = "#d62728", lwd = 1.5)
}
par(op)
dev.off()

png(file.path(plot_dir, "grouped_axis_coefficient_summary.png"), width = 1200, height = 900, res = 150)
coef_plot <- grouped_summary[grouped_summary$term %in% c("z_axis_terrigenous_mineralogical", "z_axis_redox", "z_axis_productivity_carbonate"), , drop = FALSE]
if (nrow(coef_plot) > 0) {
  ord <- order(coef_plot$estimate)
  d <- coef_plot[ord, , drop = FALSE]
  y <- seq_len(nrow(d))
  xlim <- range(c(d$ci_low, d$ci_high, d$bootstrap_ci_low, d$bootstrap_ci_high, 0), na.rm = TRUE)
  plot(d$estimate, y,
       xlim = xlim,
       yaxt = "n",
       pch = 16,
       col = "#1f77b4",
       xlab = "Coefficient estimate (standardized axes)",
       ylab = "",
       main = "Grouped-axis model coefficients")
  axis(2, at = y, labels = d$term, las = 1)
  segments(d$ci_low, y, d$ci_high, y, lwd = 2, col = "#1f77b4")
  if (all(is.finite(d$bootstrap_ci_low) | is.finite(d$bootstrap_ci_high))) {
    segments(d$bootstrap_ci_low, y + 0.12, d$bootstrap_ci_high, y + 0.12, lwd = 2, col = "#ff7f0e")
    legend("bottomright", legend = c("Model CI", "Block-bootstrap CI"), col = c("#1f77b4", "#ff7f0e"), lty = 1, bty = "n")
  }
  abline(v = 0, lty = 2, col = "#d62728")
} else {
  plot.new(); text(0.5, 0.5, "No grouped-axis coefficients to plot")
}
dev.off()

png(file.path(plot_dir, "nested_model_gain_plot.png"), width = 1200, height = 800, res = 150)
if (nrow(metrics) > 0) {
  plot(metrics$step, metrics$corr_sq_observed_vs_fitted,
       type = "b", pch = 16, lwd = 2, col = "#1f77b4",
       xaxt = "n",
       xlab = "Nested model step",
       ylab = "corr^2(observed, fitted)",
       main = "Nested model gain (productivity axis added last)")
  axis(1, at = metrics$step, labels = c("Depth", "+Terrig", "+Redox", "+Productivity"))
  if (nrow(nested_cmp) >= 4 && is.finite(nested_cmp$p_value[4])) {
    txt <- paste0("Final step p=", signif(nested_cmp$p_value[4], 3))
    text(x = metrics$step[4], y = metrics$corr_sq_observed_vs_fitted[4], labels = txt, pos = 3, cex = 0.9)
  }
} else {
  plot.new(); text(0.5, 0.5, "No nested metrics available")
}
dev.off()

# -------------------------------
# Caveats and assumptions
# -------------------------------
final_step <- nested_cmp[nrow(nested_cmp), , drop = FALSE]
prod_added_supported <- if (nrow(final_step) == 1 && is.finite(final_step$p_value)) final_step$p_value < 0.05 else NA

caveats <- c(
  "Axes are data-reduced summaries of correlated proxies; they are not direct mechanisms.",
  "Depth structure was controlled with a spline term (depth_in_core_cm) in all nested models.",
  "Candidate overlap remains possible between productivity/carbonate and terrigenous/mineralogical groups (notably si_ti).",
  "Axis orientation was set to keep direction interpretable relative to predefined anchor variables; signs should be interpreted with loadings table.",
  "Downstream sequencing-output variables were excluded from explanatory modeling.",
  ifelse(isTRUE(prod_added_supported),
         "In this framework, the productivity/carbonate axis added additional explanatory information beyond depth+terrigenous+redox (associational, not causal).",
         "In this framework, the productivity/carbonate axis did not show independent additional explanatory support beyond depth+terrigenous+redox.")
)

# -------------------------------
# Write outputs
# -------------------------------
write_tsv(data.frame(parameter = names(params), value = unlist(params), stringsAsFactors = FALSE), file.path(out_dir, "analysis_parameters.tsv"))
write_tsv(ratio_rule_table, file.path(out_dir, "xrf_ratio_small_value_rules.tsv"))

write_tsv(axis_membership_out, file.path(out_dir, "axis_variable_membership.tsv"))
write_tsv(axis_loadings_out, file.path(out_dir, "axis_loadings_or_weights.tsv"))
write_tsv(axis_scores_out, file.path(out_dir, "st13_axis_scores_table.tsv"))
write_tsv(screen, file.path(out_dir, "axis_screening_results.tsv"))
write_tsv(grouped_summary, file.path(out_dir, "grouped_axis_model_summary.tsv"))
write_tsv(metrics, file.path(out_dir, "grouped_axis_model_fit_metrics.tsv"))
write_tsv(nested_cmp, file.path(out_dir, "nested_axis_model_comparison.tsv"))
write_tsv(boot_coef, file.path(out_dir, "axis_bootstrap_summary.tsv"))
write_tsv(axis_method_summary, file.path(out_dir, "axis_method_summary.tsv"))
writeLines(c(caveats, "", "Axis construction notes:", axis_notes), con = file.path(out_dir, "analysis_caveats.txt"))

message("ST13 grouped proxy-axis analysis complete.")
message("Outputs written to: ", out_dir)
