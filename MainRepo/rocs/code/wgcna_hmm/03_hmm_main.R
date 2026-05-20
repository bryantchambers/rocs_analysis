#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(depmixS4)
})

source(file.path("code", "wgcna_hmm", "00_config.R"))
set.seed(PARAMS$seed)

meta_main <- fread(file.path(DIRS$main, "sample_metadata_main.tsv"))
MEs <- fread(file.path(DIRS$main, "module_eigengenes_main.tsv"))
dt <- merge(MEs, meta_main[, .(label, core, y_bp, mis, sst)], by.x = "sample", by.y = "label")
dt[, age_kyr := y_bp / 1000]
setnames(dt, "mis", "d18O")

me_cols <- setdiff(grep("^ME", names(dt), value = TRUE), "MEgrey")
train_cores <- PARAMS$training_cores
valid_core <- PARAMS$validation_core

train_dt <- dt[core %in% train_cores]
valid_dt <- dt[core == valid_core]
train_dt[, core_order := match(core, train_cores)]
setorder(train_dt, core_order, age_kyr)
train_dt[, core_order := NULL]
setorder(valid_dt, age_kyr)

train_raw <- as.matrix(train_dt[, ..me_cols])
valid_raw <- as.matrix(valid_dt[, ..me_cols])

train_scaled_obj <- corewise_z_transform(
  mat = train_raw,
  core_ids = train_dt$core,
  feature_names = me_cols,
  variant = "corewise_z"
)
valid_scaled_obj <- corewise_z_transform(
  mat = valid_raw,
  core_ids = valid_dt$core,
  feature_names = me_cols,
  variant = "corewise_z"
)

scaled <- list(
  train = train_scaled_obj$scaled,
  valid = valid_scaled_obj$scaled,
  params = rbindlist(list(train_scaled_obj$params, valid_scaled_obj$params), use.names = TRUE)
)

harmonized_train <- cbind(
  train_dt[, .(sample, core, age_kyr, d18O, sst)],
  as.data.table(scaled$train)
)
harmonized_train[, analysis_partition := "training"]

harmonized_valid <- cbind(
  valid_dt[, .(sample, core, age_kyr, d18O, sst)],
  as.data.table(scaled$valid)
)
harmonized_valid[, analysis_partition := "validation"]

harmonized_main <- rbindlist(list(harmonized_train, harmonized_valid), use.names = TRUE)
harmonized_main[, harmonization_mode := "corewise_z"]

pca <- prcomp(scaled$train, center = TRUE, scale. = TRUE)
n_pcs <- 3
train_scores <- as.data.table(predict(pca, newdata = scaled$train)[, 1:n_pcs, drop = FALSE])
valid_scores <- as.data.table(predict(pca, newdata = scaled$valid)[, 1:n_pcs, drop = FALSE])
setnames(train_scores, paste0("PC", 1:n_pcs))
setnames(valid_scores, paste0("PC", 1:n_pcs))

train_pca <- cbind(train_dt[, .(sample, core, age_kyr, d18O, sst)], train_scores)
valid_pca <- cbind(valid_dt[, .(sample, core, age_kyr, d18O, sst)], valid_scores)

train_ntimes <- train_pca[, .N, by = core][order(match(core, train_cores)), N]

fit_hmm <- function(dat, ntimes, k, n_iter) {
  mod <- depmix(
    list(PC1 ~ 1, PC2 ~ 1, PC3 ~ 1),
    data = dat,
    nstates = k,
    ntimes = ntimes,
    family = list(gaussian(), gaussian(), gaussian())
  )
  tryCatch(fit(mod, emcontrol = em.control(maxit = n_iter, random.start = TRUE), verbose = FALSE), error = function(e) NULL)
}

hmm_fits <- lapply(2:5, function(k) {
  best <- NULL
  for (i in seq_len(PARAMS$hmm_n_start)) {
    fi <- fit_hmm(train_pca, train_ntimes, k, PARAMS$hmm_n_iter)
    if (!is.null(fi) && (is.null(best) || BIC(fi) < BIC(best))) best <- fi
  }
  best
})
names(hmm_fits) <- paste0("K", 2:5)

bic_dt <- data.table(
  K = 2:5,
  BIC = sapply(hmm_fits, function(f) if (is.null(f)) NA_real_ else BIC(f)),
  logLik = sapply(hmm_fits, function(f) if (is.null(f)) NA_real_ else as.numeric(logLik(f)))
)

best_k <- bic_dt[!is.na(BIC)][which.min(BIC), K]
if (!is.na(best_k) && best_k != PARAMS$hmm_k_preferred) {
  delta <- bic_dt[K == PARAMS$hmm_k_preferred, BIC] - bic_dt[K == best_k, BIC]
  if (is.finite(delta) && delta < 10) best_k <- PARAMS$hmm_k_preferred
}
best_fit <- hmm_fits[[paste0("K", best_k)]]
if (is.null(best_fit)) stop("No successful HMM fit for selected K.")

decode_states <- function(model_obj, n_expected) {
  pv <- tryCatch(posterior(model_obj, type = "viterbi"), error = function(e) NULL)
  if (!is.null(pv) && "state" %in% names(pv) && length(pv$state) == n_expected) return(as.integer(pv$state))
  pg <- tryCatch(posterior(model_obj, type = "global"), error = function(e) NULL)
  if (!is.null(pg) && "state" %in% names(pg) && length(pg$state) == n_expected) return(as.integer(pg$state))
  stop("State decoding failed.")
}

train_pca[, state := decode_states(best_fit, nrow(train_pca))]

valid_mod <- setpars(
  depmix(
    list(PC1 ~ 1, PC2 ~ 1, PC3 ~ 1),
    data = rbindlist(list(train_pca[, .(PC1, PC2, PC3)], valid_pca[, .(PC1, PC2, PC3)]), use.names = TRUE),
    nstates = best_k,
    ntimes = c(train_ntimes, nrow(valid_pca)),
    family = list(gaussian(), gaussian(), gaussian())
  ),
  getpars(best_fit)
)

all_states <- decode_states(valid_mod, nrow(train_pca) + nrow(valid_pca))
valid_pca[, state := all_states[(nrow(train_pca) + 1):(nrow(train_pca) + nrow(valid_pca))]]

state_train <- train_pca[, .(mean_d18O = mean(d18O, na.rm = TRUE), mean_sst = mean(sst, na.rm = TRUE), n_train = .N), by = state]
glacial <- state_train[which.max(mean_d18O), state]
state_train[, label := ifelse(state == glacial, "G-A", paste0("IG-", LETTERS[rank(-mean_d18O, ties.method = "first")][state]))]

states_all <- rbindlist(list(train_pca, valid_pca), use.names = TRUE)
states_all <- merge(states_all, state_train[, .(state, label)], by = "state", all.x = TRUE)

state_me <- merge(
  states_all[, .(sample, state, state_label = label, core, age_kyr)],
  dt[, c("sample", me_cols), with = FALSE],
  by = "sample",
  all.x = TRUE
)
state_fingerprints_main <- state_me[, c(
  list(
    n_samples = .N,
    mean_age_kyr = mean(age_kyr, na.rm = TRUE)
  ),
  lapply(.SD, mean, na.rm = TRUE)
), by = .(state, state_label), .SDcols = me_cols]
setorder(state_fingerprints_main, state)

r1 <- states_all[core == "GeoB25202_R1", .(age_kyr, state_r1 = state)]
r2 <- states_all[core == valid_core, .(age_kyr, state_r2 = state)]
if (nrow(r1) >= 2 && nrow(r2) >= 2) {
  setorder(r1, age_kyr)
  setorder(r2, age_kyr)
  r2_interp <- approx(r2$age_kyr, r2$state_r2, xout = r1$age_kyr, method = "constant", f = 0, rule = 2)$y
  concord <- data.table(age_kyr = r1$age_kyr, state_r1 = r1$state_r1, state_r2_interp = as.integer(round(r2_interp)))
  concord[, match := state_r1 == state_r2_interp]
  concord_summary <- data.table(
    metric = c("n_compared", "n_match", "pct_match"),
    value = c(nrow(concord), sum(concord$match), mean(concord$match) * 100)
  )
} else {
  concord <- data.table()
  concord_summary <- data.table(metric = c("n_compared", "n_match", "pct_match"), value = c(NA_real_, NA_real_, NA_real_))
}

fwrite(states_all[, .(sample, core, age_kyr, d18O, sst, PC1, PC2, PC3, state, label)], file.path(DIRS$main, "hmm_states_main.tsv"), sep = "\t")
fwrite(state_fingerprints_main, file.path(DIRS$main, "state_fingerprints_main.tsv"), sep = "\t")
fwrite(bic_dt, file.path(DIRS$main, "hmm_bic_main.tsv"), sep = "\t")
fwrite(concord, file.path(DIRS$main, "state_concordance_R1_vs_R2.tsv"), sep = "\t")
fwrite(concord_summary, file.path(DIRS$main, "state_concordance_R1_vs_R2_summary.tsv"), sep = "\t")
fwrite(harmonized_main, file.path(DIRS$main, "corewise_z_eigengenes_main_harmonized.tsv"), sep = "\t")
fwrite(scaled$params, file.path(DIRS$main, "corewise_z_scaling_parameters_main.tsv"), sep = "\t")
fwrite(as.data.table(pca$rotation[, 1:n_pcs, drop = FALSE], keep.rownames = "eigengene"), file.path(DIRS$main, "pca_loadings_main.tsv"), sep = "\t")
saveRDS(best_fit, file.path(DIRS$main, "hmm_model_main.rds"))
saveRDS(list(pca = pca, hmm = hmm_fits, selected_k = best_k, feature_cols = me_cols), file.path(DIRS$main, "main_hmm_artifacts.rds"))

write_run_metadata(
  file.path(DIRS$main, "03_hmm_main_run_metadata.tsv"),
  "03_hmm_main.R",
  extra = list(
    training_cores = train_cores,
    validation_core = valid_core,
    selected_k = best_k,
    harmonization = "corewise_z_only_no_pre_residualization",
    feature_input = "module_eigengenes_raw"
  )
)
write_session_info(file.path(DIRS$main, "03_hmm_main_sessionInfo.txt"))

log_msg("03_hmm_main complete.")
