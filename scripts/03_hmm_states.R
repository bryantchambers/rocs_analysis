#!/usr/bin/env Rscript
# 03_hmm_states.R — HMM ecological state discovery from module eigengenes
#
# Pipeline:
#   1. Residualise eigengenes within core (remove core offsets)
#   2. PCA on training cores; project validation/all samples
#   3. Fit multi-sequence HMM (depmixS4) for K = 2..5 on training cores
#   4. Score held-out validation core with fixed parameters
#   5. Select K using train BIC + held-out likelihood/stability
#   6. Refit selected K on all cores for final state outputs
#   7. Save state assignments + fingerprints + validation metrics
#
# Inputs:  results/stage1/wgcna/module_eigengenes.tsv
#          results/stage1/sample_metadata_stage1.tsv
# Outputs: results/hmm/
#   hmm_states.tsv              sample → HMM state label + d18O
#   state_fingerprints.tsv      per-state mean eigengene (ecological fingerprint)
#   state_mis_crosstab.tsv      state × MIS stage cross-tabulation
#   bic_comparison.tsv          BIC per K (training cores)
#   hmm_validation_metrics.tsv  held-out performance/stability per K
#   hmm_model_selection.tsv     hybrid model-selection table
#   pca_loadings.tsv            PCA loadings used for HMM input
#   hmm_model.rds               best-fit depmixS4 model object
#   all_hmm_models.rds          all K=2..5 fitted models (for inspection)

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(depmixS4)
})

source(here("config.R"))
set.seed(PARAMS$seed)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

WGCNA_OUT <- file.path(RESULTS$stage1, "wgcna")

# ── 1. Load eigengenes + metadata ─────────────────────────────────────────────

log_msg("Loading eigengenes and metadata...")

MEs  <- fread(file.path(WGCNA_OUT, "module_eigengenes.tsv"))
meta <- fread(file.path(RESULTS$stage1, "sample_metadata_stage1.tsv"))

# Merge
dt <- merge(MEs, meta[, .(label, core, y_bp, mis)],
            by.x = "sample", by.y = "label")
dt[, age_kyr := y_bp / 1000]
setnames(dt, "mis", "d18O")
setorder(dt, core, age_kyr)

# Module eigengene columns (exclude grey, sample, metadata)
me_cols <- grep("^ME(?!grey)", names(dt), value = TRUE, perl = TRUE)
log_msg(sprintf("  Using %d module eigengenes: %s", length(me_cols),
                paste(me_cols, collapse = ", ")))

# ── 2. Residualise within core ────────────────────────────────────────────────

log_msg("Residualising eigengenes within core...")

for (col in me_cols) {
  core_means <- dt[, .(core_mean = mean(get(col), na.rm = TRUE)), by = core]
  dt[core_means, paste0(col, "_resid") := get(col) - i.core_mean, on = "core"]
}
resid_cols <- paste0(me_cols, "_resid")

# ── 3. PCA on residualised eigengenes ─────────────────────────────────────────

log_msg("PCA on residualised eigengenes...")

resid_mat <- as.matrix(dt[, ..resid_cols])
train_idx <- dt$core %in% PARAMS$stage1_cores
valid_idx <- dt$core == PARAMS$validation_core

pca <- prcomp(resid_mat[train_idx, , drop = FALSE], scale. = TRUE)

# Retain first 3 PCs. In the original analysis these captured ≥ 70% of variance.
# Check the log output above; if PC1–PC3 explain substantially less (< 60%),
# consider increasing n_pcs or inspecting whether modules have changed structure.
n_pcs    <- 3
pc_scores <- as.data.table(predict(pca, newdata = resid_mat)[, 1:n_pcs, drop = FALSE])
names(pc_scores) <- paste0("PC", 1:n_pcs)

pc_var <- summary(pca)$importance[2, 1:n_pcs]
log_msg(sprintf("  PC1–PC3 variance explained: %.1f%% / %.1f%% / %.1f%%",
                pc_var[1]*100, pc_var[2]*100, pc_var[3]*100))

pca_dt <- cbind(dt[, .(sample, core, age_kyr, d18O)], pc_scores)

# Save loadings
loadings_dt <- as.data.table(pca$rotation[, 1:n_pcs], keep.rownames = "eigengene")
fwrite(loadings_dt, file.path(RESULTS$hmm, "pca_loadings.tsv"), sep = "\t")

# ── 4. Fit HMM for K = 2..5 ──────────────────────────────────────────────────

log_msg("Fitting HMMs for K = 2..5...")

# depmixS4: treat each core as a separate sequence
pca_train <- pca_dt[core %in% PARAMS$stage1_cores]
pca_valid <- pca_dt[core == PARAMS$validation_core]

train_lengths <- pca_train[, .N, by = core][order(match(core, PARAMS$stage1_cores)), N]
valid_lengths <- pca_valid[, .N, by = core][order(match(core, PARAMS$validation_core)), N]
all_lengths <- pca_dt[, .N, by = core][order(match(core, PARAMS$all_cores)), N]

decode_metrics <- function(state_vec) {
  n <- length(state_vec)
  if (n <= 1) {
    return(list(self_transition_rate = NA_real_, switches_per_100 = NA_real_))
  }
  transitions <- sum(state_vec[-1] != state_vec[-n], na.rm = TRUE)
  list(
    self_transition_rate = 1 - (transitions / (n - 1)),
    switches_per_100 = transitions / (n - 1) * 100
  )
}

fit_hmm <- function(k, data, ntimes) {
  mod <- depmix(
    list(PC1 ~ 1, PC2 ~ 1, PC3 ~ 1),
    data      = data,
    nstates   = k,
    ntimes    = ntimes,
    family    = list(gaussian(), gaussian(), gaussian())
  )
  tryCatch(
    fit(mod, emcontrol = em.control(maxit = PARAMS$hmm_n_iter, random.start = TRUE),
        verbose = FALSE),
    error = function(e) NULL
  )
}

hmm_fits <- lapply(2:5, function(k) {
  log_msg(sprintf("  K = %d ...", k))
  best <- NULL
  for (i in seq_len(PARAMS$hmm_n_start)) {
    f <- fit_hmm(k, data = pca_train, ntimes = train_lengths)
    if (!is.null(f)) {
      if (is.null(best) || BIC(f) < BIC(best)) best <- f
    }
  }
  best
})
names(hmm_fits) <- paste0("K", 2:5)

# BIC table
bic_dt <- data.table(
  K   = 2:5,
  BIC = sapply(hmm_fits, function(f) if (!is.null(f)) BIC(f) else NA_real_),
  logLik = sapply(hmm_fits, function(f) if (!is.null(f)) logLik(f) else NA_real_)
)
fwrite(bic_dt, file.path(RESULTS$hmm, "bic_comparison.tsv"), sep = "\t")

validation_rows <- lapply(seq_along(hmm_fits), function(i) {
  k <- 1 + i
  fit_k <- hmm_fits[[i]]
  if (is.null(fit_k) || nrow(pca_valid) == 0) {
    return(data.table(
      K = k,
      validation_logLik = NA_real_,
      validation_logLik_per_sample = NA_real_,
      validation_mean_max_posterior = NA_real_,
      validation_self_transition_rate = NA_real_,
      validation_switches_per_100 = NA_real_,
      validation_n_observed_states = NA_integer_
    ))
  }

  valid_mod <- depmix(
    list(PC1 ~ 1, PC2 ~ 1, PC3 ~ 1),
    data    = pca_valid,
    nstates = k,
    ntimes  = valid_lengths,
    family  = list(gaussian(), gaussian(), gaussian())
  )
  valid_mod <- setpars(valid_mod, getpars(fit_k))
  valid_ll <- as.numeric(logLik(valid_mod))
  valid_post <- posterior(valid_mod, type = "viterbi")
  max_post <- valid_post[, grep("^S", names(valid_post)), drop = FALSE]
  decoded <- valid_post$state
  dm <- decode_metrics(decoded)

  data.table(
    K = k,
    validation_logLik = valid_ll,
    validation_logLik_per_sample = valid_ll / nrow(pca_valid),
    validation_mean_max_posterior = mean(apply(max_post, 1, max), na.rm = TRUE),
    validation_self_transition_rate = dm$self_transition_rate,
    validation_switches_per_100 = dm$switches_per_100,
    validation_n_observed_states = uniqueN(decoded)
  )
})
valid_dt <- rbindlist(validation_rows, fill = TRUE)
fwrite(valid_dt, file.path(RESULTS$hmm, "hmm_validation_metrics.tsv"), sep = "\t")

sel_dt <- merge(bic_dt, valid_dt, by = "K", all = TRUE)
best_bic_k <- sel_dt[!is.na(BIC)][which.min(BIC), K]
sel_dt[, delta_bic := BIC - min(BIC, na.rm = TRUE)]
sel_dt[, bic_ambiguous := delta_bic <= 10]

# Hybrid selection: among BIC-ambiguous K, choose best held-out logLik/sample.
# If tied, prefer cleaner transitions (higher self-transition rate).
cand <- sel_dt[bic_ambiguous == TRUE & !is.na(validation_logLik_per_sample)]
if (nrow(cand) == 0) cand <- sel_dt[!is.na(BIC)]
setorder(cand, -validation_logLik_per_sample, -validation_self_transition_rate, BIC)
best_k <- cand[1, K]

fwrite(sel_dt, file.path(RESULTS$hmm, "hmm_model_selection.tsv"), sep = "\t")

log_msg(sprintf("  Selected K = %d", best_k))
log_msg(sprintf("  Best train-BIC K = %d; hybrid selected K = %d", best_bic_k, best_k))

log_msg("Refitting selected K on all cores for final outputs...")
best_fit <- NULL
for (i in seq_len(PARAMS$hmm_n_start)) {
  f <- fit_hmm(best_k, data = pca_dt, ntimes = all_lengths)
  if (!is.null(f)) {
    if (is.null(best_fit) || BIC(f) < BIC(best_fit)) best_fit <- f
  }
}
if (is.null(best_fit)) stop("Failed to fit final HMM on all cores.")

# ── 5. Extract state assignments ──────────────────────────────────────────────

log_msg("Extracting state assignments...")

pca_dt[, state := posterior(best_fit, type = "viterbi")$state]

# Characterise states by d18O (link to glacial/interglacial)
state_d18o <- pca_dt[, .(mean_d18O = mean(d18O, na.rm = TRUE),
                          n = .N), by = state]
glacial_state <- state_d18o[which.max(mean_d18O), state]
log_msg(sprintf("  Glacial state (highest mean d18O): state %d", glacial_state))

# Label states deterministically with unique labels for any selected K.
# Glacial state stays G-A. Remaining states are ordered by mean d18O (high to low)
# and labeled IG-B, IG-C, ...
state_d18o[, label := NA_character_]
state_d18o[state == glacial_state, label := "G-A"]
ig_states <- state_d18o[state != glacial_state][order(-mean_d18O)]
if (nrow(ig_states) > 0) {
  ig_labels <- paste0("IG-", LETTERS[seq.int(2, length.out = nrow(ig_states))])
  state_d18o[ig_states, label := ig_labels, on = "state"]
}

# ── 6. State fingerprints ─────────────────────────────────────────────────────

me_means <- merge(pca_dt[, .(sample, state)], dt[, c("sample", me_cols), with = FALSE],
                  by = "sample")
fingerprints <- me_means[, lapply(.SD, mean, na.rm = TRUE), by = state, .SDcols = me_cols]
fingerprints <- merge(fingerprints, state_d18o[, .(state, label, mean_d18O, n)], by = "state")

# ── 7. MIS cross-tabulation ───────────────────────────────────────────────────

# δ¹⁸O ≥ 4.0 ‰ as glacial threshold follows the LR04 benthic stack convention
# (Lisiecki & Raymo 2005); interglacial periods typically < 3.5 ‰.
pca_dt[, mis_stage := ifelse(d18O >= 4.0, "glacial", "interglacial")]
mis_cross <- pca_dt[, .N, by = .(state, mis_stage)]
mis_wide  <- dcast(mis_cross, state ~ mis_stage, value.var = "N", fill = 0)
mis_wide[, pct_glacial := glacial / (glacial + interglacial) * 100]

# ── 8. Save ───────────────────────────────────────────────────────────────────

log_msg("Saving outputs...")

states_out <- merge(pca_dt[, .(sample, core, age_kyr, d18O, state)],
                    state_d18o[, .(state, label)], by = "state")

fwrite(states_out,  file.path(RESULTS$hmm, "hmm_states.tsv"),      sep = "\t")
fwrite(fingerprints, file.path(RESULTS$hmm, "state_fingerprints.tsv"), sep = "\t")
fwrite(mis_wide,    file.path(RESULTS$hmm, "state_mis_crosstab.tsv"), sep = "\t")
saveRDS(best_fit,   file.path(RESULTS$hmm, "hmm_model.rds"))
saveRDS(hmm_fits,   file.path(RESULTS$hmm, "all_hmm_models.rds"))

log_msg("Done. Outputs in ", RESULTS$hmm)
