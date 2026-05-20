#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(depmixS4)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)

clr <- readRDS(file.path(DIRS$main, "clr_matrix_train_centered.rds"))
meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
meta_ext <- fread(file.path(DIRS$extended, "sample_metadata_extended.tsv"))
mods <- fread(file.path(DIRS$main, "module_assignments.tsv"))
MEs_train <- fread(file.path(DIRS$main, "module_eigengenes_training.tsv"))
main_states <- fread(file.path(DIRS$main, "hmm_states_main.tsv"))
art <- readRDS(file.path(DIRS$main, "main_hmm_artifacts.rds"))
main_scaling_params <- fread(file.path(DIRS$main, "corewise_z_scaling_parameters_main.tsv"))

module_colors <- setNames(mods$module, mods$taxon)
train_ids <- intersect(MEs_train$sample, rownames(clr))
if (length(train_ids) == 0) stop("No overlapping training IDs between module eigengenes and CLR matrix.")
MEs_train <- MEs_train[match(train_ids, sample)]
train_expr <- clr[train_ids, names(module_colors), drop = FALSE]
MEs_train_mat <- as.data.frame(MEs_train[, setdiff(names(MEs_train), "sample"), with = FALSE])

project_to_training_basis <- function(expr_train, expr_new, module_colors, me_train_df) {
  modules <- sort(unique(module_colors))
  out <- vector("list", length(modules))
  names(out) <- paste0("ME", modules)

  for (i in seq_along(modules)) {
    mod <- modules[i]
    me_name <- paste0("ME", mod)
    taxa <- names(module_colors)[module_colors == mod]
    xtr <- as.matrix(expr_train[, taxa, drop = FALSE])
    xnw <- as.matrix(expr_new[, taxa, drop = FALSE])
    mu <- colMeans(xtr, na.rm = TRUE)
    sdv <- apply(xtr, 2, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    ztr <- sweep(sweep(xtr, 2, mu, "-"), 2, sdv, "/")
    znw <- sweep(sweep(xnw, 2, mu, "-"), 2, sdv, "/")
    if (ncol(ztr) == 1) {
      tr_pc1 <- as.numeric(ztr[, 1])
      nw_pc1 <- as.numeric(znw[, 1])
    } else {
      p <- prcomp(ztr, center = FALSE, scale. = FALSE)
      l1 <- p$rotation[, 1]
      tr_pc1 <- as.numeric(ztr %*% l1)
      nw_pc1 <- as.numeric(znw %*% l1)
    }

    sign_cor <- suppressWarnings(cor(tr_pc1, me_train_df[[me_name]], use = "pairwise.complete.obs"))
    if (is.finite(sign_cor) && sign_cor < 0) nw_pc1 <- -nw_pc1

    out[[i]] <- nw_pc1
  }
  as.data.table(out)
}

ext_ids <- intersect(meta_ext$label, rownames(clr))
expr_ext <- clr[ext_ids, names(module_colors), drop = FALSE]
me_ext <- project_to_training_basis(train_expr, expr_ext, module_colors, MEs_train_mat)
me_ext[, sample := rownames(expr_ext)]
setcolorder(me_ext, c("sample", setdiff(names(me_ext), "sample")))

dt <- merge(me_ext, meta_ext[, .(label, core, y_bp, mis, sst)], by.x = "sample", by.y = "label")
dt[, age_kyr := y_bp / 1000]
setnames(dt, "mis", "d18O")

me_cols <- setdiff(grep("^ME", names(dt), value = TRUE), "MEgrey")

feature_cols <- if (!is.null(art$feature_cols)) art$feature_cols else if (!is.null(art$resid_cols)) art$resid_cols else me_cols
missing_feature_cols <- setdiff(feature_cols, names(dt))
if (length(missing_feature_cols) > 0) {
  stop(
    "Extended projection is missing required feature columns from main HMM artifacts: ",
    paste(missing_feature_cols, collapse = ", ")
  )
}

build_extended_scaling_from_imported_main <- function(d, feature_cols, main_params, main_max_age_kyr) {
  x <- as.matrix(d[, ..feature_cols])
  core_ids <- as.character(d$core)
  main_params_use <- copy(main_params)

  required_cols <- c("core", "feature", "mu", "sd")
  missing_required <- setdiff(required_cols, names(main_params_use))
  if (length(missing_required) > 0) {
    stop("Main scaling parameter table is missing required columns: ", paste(missing_required, collapse = ", "))
  }

  main_params_use <- main_params_use[feature %in% feature_cols, .(core, feature, mu, sd)]

  if (nrow(main_params_use) == 0) {
    warning("No matching rows found in imported main scaling parameters for requested features; using conservative fallback.")
  }

  ref_mask <- d$age_kyr <= main_max_age_kyr

  if (!any(ref_mask)) {
    stop("No main-window rows available for conservative fallback scaling in extended projection.")
  }

  ref_global <- x[ref_mask, , drop = FALSE]
  global_mu <- colMeans(ref_global, na.rm = TRUE)
  global_sd <- apply(ref_global, 2, sd, na.rm = TRUE)
  global_sd[!is.finite(global_sd) | global_sd == 0] <- 1

  out <- x
  param_rows <- list()
  cores <- unique(core_ids)

  for (co in cores) {
    idx_core <- which(core_ids == co)
    idx_ref_core <- idx_core[d$age_kyr[idx_core] <= main_max_age_kyr]

    mu <- numeric(length(feature_cols))
    sdv <- numeric(length(feature_cols))
    source <- character(length(feature_cols))
    n_ref <- integer(length(feature_cols))
    fallback_used <- logical(length(feature_cols))

    for (j in seq_along(feature_cols)) {
      feat <- feature_cols[j]

      imported_row <- main_params_use[core == co & feature == feat]
      if (nrow(imported_row) > 0 && is.finite(imported_row$mu[1]) && is.finite(imported_row$sd[1]) && imported_row$sd[1] > 0) {
        mu[j] <- as.numeric(imported_row$mu[1])
        sdv[j] <- as.numeric(imported_row$sd[1])
        source[j] <- "imported_main_params"
        n_ref[j] <- NA_integer_
        fallback_used[j] <- FALSE
      } else if (length(idx_ref_core) > 0) {
        # Conservative fallback 1: core-specific rows from main window only.
        x_ref <- x[idx_ref_core, j, drop = FALSE]
        mu[j] <- mean(x_ref[, 1], na.rm = TRUE)
        sdv[j] <- stats::sd(x_ref[, 1], na.rm = TRUE)
        source[j] <- "core_main_window_fallback"
        n_ref[j] <- length(idx_ref_core)
        fallback_used[j] <- TRUE
      } else {
        # Conservative fallback 2: pooled rows from main window only.
        mu[j] <- global_mu[j]
        sdv[j] <- global_sd[j]
        source[j] <- "global_main_window_fallback"
        n_ref[j] <- 0L
        fallback_used[j] <- TRUE
      }
    }

    sdv[!is.finite(sdv) | sdv == 0] <- 1
    x_core <- x[idx_core, , drop = FALSE]
    out[idx_core, ] <- sweep(sweep(x_core, 2, mu, "-"), 2, sdv, "/")

    param_rows[[co]] <- data.table(
      core = co,
      feature = feature_cols,
      mu = as.numeric(mu),
      sd = as.numeric(sdv),
      source = source,
      n_reference_samples = as.integer(n_ref),
      fallback_used = as.logical(fallback_used),
      variant = "corewise_z_main_params_imported"
    )
  }

  list(
    scaled = out,
    params = rbindlist(param_rows, use.names = TRUE)
  )
}

scale_with_mode <- function(d, mode) {
  if (mode == "corewise_z") {
    scaled_obj <- build_extended_scaling_from_imported_main(
      d = d,
      feature_cols = feature_cols,
      main_params = main_scaling_params,
      main_max_age_kyr = PARAMS$main_max_age_kyr
    )
    return(scaled_obj)
  }
  stop("Unknown mode: ", mode)
}

decode_for_mode <- function(mode_name) {
  scaling_obj <- scale_with_mode(dt, mode_name)
  x_scaled <- scaling_obj$scaled
  pca <- art$pca
  k <- art$selected_k
  hmm_model <- art$hmm[[paste0("K", k)]]

  sc <- as.data.table(predict(pca, newdata = x_scaled)[, 1:3, drop = FALSE])
  setnames(sc, c("PC1", "PC2", "PC3"))
  dd <- cbind(dt[, .(sample, core, age_kyr, d18O, sst)], sc)

  dd[, core_ord := match(core, PARAMS$all_cores)]
  setorder(dd, core_ord, age_kyr)
  seq_sizes <- dd[, .N, by = core][order(match(core, PARAMS$all_cores)), N]

  mod <- setpars(
    depmix(
      list(PC1 ~ 1, PC2 ~ 1, PC3 ~ 1),
      data = dd[, .(PC1, PC2, PC3)],
      nstates = k,
      ntimes = seq_sizes,
      family = list(gaussian(), gaussian(), gaussian())
    ),
    getpars(hmm_model)
  )

  states_new <- tryCatch({
    pv <- posterior(mod, type = "viterbi")
    as.integer(pv$state)
  }, error = function(e) {
    tr <- merge(
      dd[age_kyr <= PARAMS$main_max_age_kyr, .(sample, PC1, PC2, PC3)],
      main_states[, .(sample, state)],
      by = "sample",
      all = FALSE
    )
    if (nrow(tr) == 0) stop("Fallback decode failed: no overlap with main states.")

    pars <- rbindlist(lapply(sort(unique(tr$state)), function(s) {
      x <- tr[state == s]
      data.table(
        state = s,
        mu1 = mean(x$PC1, na.rm = TRUE), sd1 = pmax(sd(x$PC1, na.rm = TRUE), 1e-6),
        mu2 = mean(x$PC2, na.rm = TRUE), sd2 = pmax(sd(x$PC2, na.rm = TRUE), 1e-6),
        mu3 = mean(x$PC3, na.rm = TRUE), sd3 = pmax(sd(x$PC3, na.rm = TRUE), 1e-6)
      )
    }))
    ll <- sapply(seq_len(nrow(pars)), function(i) {
      p <- pars[i]
      dnorm(dd$PC1, p$mu1, p$sd1, log = TRUE) +
        dnorm(dd$PC2, p$mu2, p$sd2, log = TRUE) +
        dnorm(dd$PC3, p$mu3, p$sd3, log = TRUE)
    })
    if (is.null(dim(ll))) ll <- matrix(ll, ncol = 1)
    pars$state[max.col(ll, ties.method = "first")]
  })

  dd[, state_inferred := states_new]
  dd[, state := state_inferred]
  dd[, harmonization_mode := mode_name]

  harmonized_eigengenes <- cbind(
    dt[, .(sample, core, age_kyr, d18O, sst)],
    {
      x_scaled_mat <- x_scaled
      if (is.null(colnames(x_scaled_mat)) || !identical(colnames(x_scaled_mat), feature_cols)) {
        colnames(x_scaled_mat) <- feature_cols
      }
      x_scaled_dt <- as.data.table(x_scaled_mat)
      setcolorder(x_scaled_dt, feature_cols)
      if (!identical(names(x_scaled_dt), feature_cols)) {
        stop("Scaled eigengene columns are misaligned with expected feature order.")
      }
      x_scaled_dt
    }
  )
  harmonized_eigengenes[, harmonization_mode := mode_name]

  list(states = dd, scaling_params = scaling_obj$params, harmonized_eigengenes = harmonized_eigengenes)
}

modes <- c("corewise_z")
mode_states <- lapply(modes, decode_for_mode)
states_all_modes <- rbindlist(lapply(mode_states, function(x) x$states), use.names = TRUE)
scaling_params_all_modes <- rbindlist(
  Map(function(mode_name, x) {
    out <- copy(x$scaling_params)
    out[, harmonization_mode := mode_name]
    out
  }, modes, mode_states),
  use.names = TRUE
)

harmonized_eigengenes_all_modes <- rbindlist(
  lapply(mode_states, function(x) x$harmonized_eigengenes),
  use.names = TRUE
)

main_cmp <- merge(
  states_all_modes[age_kyr <= PARAMS$main_max_age_kyr, .(sample, harmonization_mode, state_extended = state)],
  main_states[, .(sample, state_main = state)],
  by = "sample",
  all.x = FALSE
)
main_cmp[, match := state_extended == state_main]

mode_eval <- main_cmp[, .(
  n_overlap = .N,
  n_match = sum(match, na.rm = TRUE),
  pct_match = mean(match, na.rm = TRUE) * 100
), by = harmonization_mode][order(-pct_match, -n_match)]

chosen_mode <- "corewise_z"

chosen_states <- states_all_modes[harmonization_mode == chosen_mode]
chosen_scaling_params <- scaling_params_all_modes[harmonization_mode == chosen_mode]
chosen_harmonized_eigengenes <- harmonized_eigengenes_all_modes[harmonization_mode == chosen_mode]

apply_older_anchor_for_core <- function(states_dt, hmm_model, selected_k, main_states_dt, anchor_core, anchor_threshold_kyr) {
  dt <- copy(states_dt)

  if (!"state_inferred" %in% names(dt)) {
    stop("Anchored older decoding requires state_inferred column in states table.")
  }

  if (!isTRUE(anchor_core %in% dt$core)) {
    return(dt)
  }

  idx_core <- which(dt$core == anchor_core)
  if (length(idx_core) == 0) return(dt)

  idx_old <- idx_core[dt$age_kyr[idx_core] > anchor_threshold_kyr]
  idx_overlap <- idx_core[dt$age_kyr[idx_core] <= anchor_threshold_kyr]
  if (length(idx_old) == 0 || length(idx_overlap) == 0) return(dt)

  dt[, core_ord_anchor := match(core, PARAMS$all_cores)]
  setorder(dt, core_ord_anchor, age_kyr)
  idx_core <- which(dt$core == anchor_core)
  idx_old <- idx_core[dt$age_kyr[idx_core] > anchor_threshold_kyr]
  idx_overlap <- idx_core[dt$age_kyr[idx_core] <= anchor_threshold_kyr]
  if (length(idx_old) == 0 || length(idx_overlap) == 0) return(dt)

  anchor_sample <- dt$sample[max(idx_overlap)]
  anchor_state <- main_states_dt[sample == anchor_sample, state]
  if (length(anchor_state) == 0 || !is.finite(anchor_state[1])) return(dt)
  anchor_state <- as.integer(anchor_state[1])
  if (!(anchor_state %in% seq_len(selected_k))) return(dt)

  transition <- t(sapply(seq_len(selected_k), function(i) as.numeric(dens(hmm_model@transition[[i]]))))
  transition[!is.finite(transition) | transition <= 0] <- 1e-12
  transition <- transition / rowSums(transition)

  mu <- matrix(NA_real_, nrow = selected_k, ncol = 3)
  sdv <- matrix(NA_real_, nrow = selected_k, ncol = 3)
  for (s in seq_len(selected_k)) {
    for (d in 1:3) {
      pars <- hmm_model@response[[s]][[d]]@parameters
      mu[s, d] <- as.numeric(pars$coefficients[[1]])
      sdv[s, d] <- as.numeric(pars$sd[[1]])
    }
  }
  sdv[!is.finite(sdv) | sdv <= 0] <- 1e-6

  x_old <- as.matrix(dt[idx_old, .(PC1, PC2, PC3)])
  n_old <- nrow(x_old)
  if (n_old == 0) return(dt)

  ll <- matrix(0, nrow = n_old, ncol = selected_k)
  for (s in seq_len(selected_k)) {
    ll[, s] <-
      dnorm(x_old[, 1], mu[s, 1], sdv[s, 1], log = TRUE) +
      dnorm(x_old[, 2], mu[s, 2], sdv[s, 2], log = TRUE) +
      dnorm(x_old[, 3], mu[s, 3], sdv[s, 3], log = TRUE)
  }

  delta <- matrix(-Inf, nrow = n_old, ncol = selected_k)
  psi <- matrix(NA_integer_, nrow = n_old, ncol = selected_k)
  delta[1, anchor_state] <- ll[1, anchor_state]

  if (n_old > 1) {
    for (tt in 2:n_old) {
      for (s in seq_len(selected_k)) {
        vals <- delta[tt - 1, ] + log(transition[, s])
        psi[tt, s] <- which.max(vals)
        delta[tt, s] <- max(vals) + ll[tt, s]
      }
    }
  }

  path <- integer(n_old)
  path[n_old] <- which.max(delta[n_old, ])
  if (n_old > 1) {
    for (tt in (n_old - 1):1) path[tt] <- psi[tt + 1, path[tt + 1]]
  }

  dt[idx_old, state_inferred := as.integer(path)]
  dt[, core_ord_anchor := NULL]
  dt
}

if (isTRUE(PARAMS$extended_anchor_older_validation_core)) {
  chosen_states <- apply_older_anchor_for_core(
    states_dt = chosen_states,
    hmm_model = art$hmm[[paste0("K", art$selected_k)]],
    selected_k = art$selected_k,
    main_states_dt = main_states,
    anchor_core = PARAMS$extended_anchor_core,
    anchor_threshold_kyr = PARAMS$main_max_age_kyr
  )
}

# ensure final pre-lock inferred state reflects any optional anchored re-decoding updates
chosen_states[, state := state_inferred]

# lock main-window states to the main analysis and infer only older-period states
main_state_lock <- unique(main_states[, .(sample, state_main = state)])
chosen_states <- merge(chosen_states, main_state_lock, by = "sample", all.x = TRUE)

chosen_states[, state_origin := fifelse(
  age_kyr <= PARAMS$main_max_age_kyr,
  "inherited_main_window",
  "inferred_extended_only"
)]

missing_main_state <- chosen_states[
  age_kyr <= PARAMS$main_max_age_kyr & is.na(state_main),
  .N
]
if (missing_main_state > 0) {
  warning(
    missing_main_state,
    " main-window samples in extended table lacked a matching main state; retaining inferred state for those rows."
  )
  chosen_states[
    age_kyr <= PARAMS$main_max_age_kyr & is.na(state_main),
    state_origin := "inferred_due_to_missing_main_state"
  ]
}

chosen_states[
  age_kyr <= PARAMS$main_max_age_kyr & !is.na(state_main),
  state := state_main
]

# map labels from main state climatology for interpretability (after lock)
state_map <- unique(main_states[, .(state, label)])
chosen_states <- merge(chosen_states, state_map, by = "state", all.x = TRUE)

lock_diag <- merge(
  chosen_states[age_kyr <= PARAMS$main_max_age_kyr, .(
    sample,
    harmonization_mode,
    state_final = state,
    state_inferred,
    state_origin
  )],
  main_states[, .(sample, state_main = state)],
  by = "sample",
  all.x = TRUE
)
lock_diag[, final_matches_main := state_final == state_main]
lock_diag[, inferred_matches_main := state_inferred == state_main]

state_origin_summary <- chosen_states[, .(
  n_samples = .N,
  n_unique_states = uniqueN(state),
  min_age_kyr = min(age_kyr, na.rm = TRUE),
  max_age_kyr = max(age_kyr, na.rm = TRUE)
), by = state_origin][order(state_origin)]

older_only <- chosen_states[age_kyr > PARAMS$main_max_age_kyr]
older_diag <- older_only[, .(
  n_samples = .N,
  n_states = uniqueN(state),
  min_age_kyr = min(age_kyr, na.rm = TRUE),
  max_age_kyr = max(age_kyr, na.rm = TRUE)
), by = core]

fwrite(me_ext, file.path(DIRS$extended, "module_eigengenes_extended_projection.tsv"), sep = "\t")
fwrite(states_all_modes, file.path(DIRS$extended, "hmm_states_extended_all_modes.tsv"), sep = "\t")
fwrite(main_cmp, file.path(DIRS$extended, "extended_mode_vs_main_overlap.tsv"), sep = "\t")
fwrite(mode_eval, file.path(DIRS$extended, "extended_harmonization_mode_evaluation.tsv"), sep = "\t")
fwrite(chosen_states, file.path(DIRS$extended, "hmm_states_extended_selected.tsv"), sep = "\t")
fwrite(older_only, file.path(DIRS$extended, "hmm_states_older_than_150kyr_selected.tsv"), sep = "\t")
fwrite(older_diag, file.path(DIRS$extended, "extended_older_period_diagnostics.tsv"), sep = "\t")
fwrite(lock_diag, file.path(DIRS$extended, "extended_locked_overlap_diagnostics.tsv"), sep = "\t")
fwrite(state_origin_summary, file.path(DIRS$extended, "extended_state_origin_summary.tsv"), sep = "\t")
fwrite(chosen_harmonized_eigengenes, file.path(DIRS$extended, "corewise_z_eigengenes_extended_harmonized.tsv"), sep = "\t")
fwrite(chosen_scaling_params, file.path(DIRS$extended, "extended_corewise_z_scaling_parameters_used.tsv"), sep = "\t")

anchor_diag <- data.table(
  anchor_enabled = isTRUE(PARAMS$extended_anchor_older_validation_core),
  anchor_core = as.character(PARAMS$extended_anchor_core),
  anchor_threshold_kyr = PARAMS$main_max_age_kyr,
  n_anchor_core_overlap = chosen_states[core == PARAMS$extended_anchor_core & age_kyr <= PARAMS$main_max_age_kyr, .N],
  n_anchor_core_older = chosen_states[core == PARAMS$extended_anchor_core & age_kyr > PARAMS$main_max_age_kyr, .N],
  n_anchor_core_states_older = chosen_states[core == PARAMS$extended_anchor_core & age_kyr > PARAMS$main_max_age_kyr, uniqueN(state)]
)
fwrite(anchor_diag, file.path(DIRS$extended, "extended_older_anchor_diagnostics.tsv"), sep = "\t")

fallback_diag <- chosen_scaling_params[, .(
  n_features = .N,
  n_features_imported = sum(source == "imported_main_params", na.rm = TRUE),
  n_features_core_fallback = sum(source == "core_main_window_fallback", na.rm = TRUE),
  n_features_global_fallback = sum(source == "global_main_window_fallback", na.rm = TRUE),
  fallback_used = any(fallback_used),
  min_reference_samples_fallback = suppressWarnings(min(n_reference_samples[fallback_used == TRUE], na.rm = TRUE)),
  max_reference_samples_fallback = suppressWarnings(max(n_reference_samples[fallback_used == TRUE], na.rm = TRUE))
), by = .(core, harmonization_mode)][order(match(core, PARAMS$all_cores))]
fallback_diag[!is.finite(min_reference_samples_fallback), min_reference_samples_fallback := NA_integer_]
fallback_diag[!is.finite(max_reference_samples_fallback), max_reference_samples_fallback := NA_integer_]
fwrite(fallback_diag, file.path(DIRS$extended, "extended_corewise_z_scaling_fallback_diagnostics.tsv"), sep = "\t")

write_run_metadata(
  file.path(DIRS$extended, "05_extended_projection_run_metadata.tsv"),
  "05_extended_projection.R",
  extra = list(
    tested_modes = modes,
    selected_mode = chosen_mode,
    selection_criterion = "operational_corewise_z_single_mode",
    overlap_state_policy = "lock_main_window_states_to_main_analysis",
    overlap_lock_threshold_kyr = PARAMS$main_max_age_kyr,
    overlap_state_origin_labels = "inherited_main_window;inferred_extended_only;inferred_due_to_missing_main_state",
    n_inherited_main_window_states = nrow(chosen_states[state_origin == "inherited_main_window"]),
    n_newly_inferred_extended_states = nrow(chosen_states[state_origin == "inferred_extended_only"]),
    n_missing_main_state_in_overlap = missing_main_state,
    older_anchor_enabled = isTRUE(PARAMS$extended_anchor_older_validation_core),
    older_anchor_core = as.character(PARAMS$extended_anchor_core),
    scaling_reference_window = "imported_main_corewise_z_scaling_parameters_main.tsv",
    scaling_fallback_rule = "if_missing_or_invalid_imported_param_use_core_main_window_then_global_main_window",
    n_cores_with_scaling_fallback = uniqueN(chosen_scaling_params[fallback_used == TRUE, core])
  )
)
write_session_info(file.path(DIRS$extended, "05_extended_projection_sessionInfo.txt"))

log_msg("05_extended_projection complete. Selected mode: ", chosen_mode)
