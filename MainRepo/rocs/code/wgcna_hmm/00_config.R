#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

env_flag <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  val <- tolower(trimws(raw))
  if (val %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  warning(sprintf("Unrecognized boolean env var %s='%s'; using default=%s", name, raw, default))
  default
}

env_csv <- function(name, default = character()) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  vals <- trimws(unlist(strsplit(raw, ",", fixed = TRUE)))
  vals <- vals[nzchar(vals)]
  unique(vals)
}

env_string <- function(name, default = "") {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  trimws(raw)
}

env_numeric <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val)) {
    warning(sprintf("Unrecognized numeric env var %s='%s'; using default=%s", name, raw, as.character(default)))
    return(default)
  }
  val
}

env_integer <- function(name, default) {
  val <- env_numeric(name, default = default)
  if (!is.finite(val)) return(default)
  as.integer(round(val))
}

env_choice <- function(name, choices, default) {
  raw <- tolower(env_string(name, default = default))
  if (!raw %in% choices) {
    warning(sprintf("Unrecognized %s='%s'; using default='%s'", name, raw, default))
    return(default)
  }
  raw
}

REPO_ROOT <- getwd()

if (!file.exists(file.path(REPO_ROOT, "data", "metadata_v5.tsv"))) {
  stop("Run scripts from repo root so relative paths resolve (expected data/metadata_v5.tsv).")
}

results_root_env <- env_string("WGCNA_HMM_RESULTS_ROOT", default = "")
results_suffix_env <- env_string("WGCNA_HMM_RESULTS_SUFFIX", default = "")
if (nzchar(results_root_env) && nzchar(results_suffix_env)) {
  warning("Both WGCNA_HMM_RESULTS_ROOT and WGCNA_HMM_RESULTS_SUFFIX set; using WGCNA_HMM_RESULTS_ROOT.")
}

RESULTS_ROOT <- if (nzchar(results_root_env)) {
  if (grepl("^/", results_root_env)) results_root_env else file.path(REPO_ROOT, results_root_env)
} else if (nzchar(results_suffix_env)) {
  file.path(REPO_ROOT, "results", paste0("wgcna_hmm_", results_suffix_env))
} else {
  file.path(REPO_ROOT, "results", "wgcna_hmm")
}

DIRS <- list(
  code = file.path(REPO_ROOT, "code", "wgcna_hmm"),
  results = RESULTS_ROOT,
  main = file.path(RESULTS_ROOT, "main"),
  main_tea = file.path(RESULTS_ROOT, "main", "tea"),
  extended = file.path(RESULTS_ROOT, "extended"),
  reports = file.path(RESULTS_ROOT, "reports"),
  logs = file.path(RESULTS_ROOT, "logs")
)

module_method_env <- tolower(env_string("WGCNA_HMM_MODULE_METHOD", default = "wgcna"))
if (!module_method_env %in% c("wgcna", "leiden")) {
  warning(sprintf("Unrecognized WGCNA_HMM_MODULE_METHOD='%s'; using 'wgcna'.", module_method_env))
  module_method_env <- "wgcna"
}

INPUTS <- list(
  tax_damage = "/projects/caeg/people/ngm902/apps/repos/rocs/results/microbial/damage/damage-classification-depositional/dmg-summary-ssp-damage-classification-depositional.tsv.gz",
  metadata = file.path(REPO_ROOT, "data", "metadata_v5.tsv"),
  xrf = file.path(REPO_ROOT, "data", "combined_xrf_not_normalized.tsv"),
  kegg_mods = "/projects/caeg/people/ngm902/apps/repos/rocs/data/functional/kegg-modules-summary-rocs.tsv.gz",
  prokaryote_function = file.path(REPO_ROOT, "results", "microbial", "wgcna", "classification", "prokaryote_function_assigned.tsv")
)

PARAMS <- list(
  seed = 42,
  excluded_samples = c("LV3003046968"),
  all_cores = c("ST8", "ST13", "GeoB25202_R1", "GeoB25202_R2"),
  training_cores = c("ST8", "ST13", "GeoB25202_R1"),
  validation_core = "GeoB25202_R2",
  main_max_age_kyr = 150,
  extended_max_age_kyr = c(
    ST8 = 150,
    ST13 = 300,
    GeoB25202_R1 = 600,
    GeoB25202_R2 = 600
  ),
  prevalence_min_samples = 10,
  clr_pseudocount = 0.5,
  wgcna_min_module_size = 20,
  wgcna_deep_split = 2,
  wgcna_merge_cut_height = 0.15,
  soft_power_target_r2 = 0.80,
  module_method = module_method_env,
  leiden_resolution = env_numeric("WGCNA_HMM_LEIDEN_RESOLUTION", default = 1.0),
  leiden_resolution_strict = env_flag("WGCNA_HMM_LEIDEN_RESOLUTION_STRICT", default = FALSE),
  leiden_graph_k_neighbors = max(5L, env_integer("WGCNA_HMM_LEIDEN_K_NEIGHBORS", default = 25)),
  leiden_graph_symmetrization = env_choice("WGCNA_HMM_LEIDEN_GRAPH_SYMMETRIZATION", choices = c("union", "mutual"), default = "union"),
  leiden_consensus_matrix_interpretation = env_choice("WGCNA_HMM_LEIDEN_CONSENSUS_MATRIX_INTERPRETATION", choices = c("auto", "similarity", "dissimilarity"), default = "auto"),
  leiden_merge_eigengene_cor = env_numeric("WGCNA_HMM_LEIDEN_MERGE_EIGENGENE_COR", default = NA_real_),
  leiden_merge_max_iter = max(1L, env_integer("WGCNA_HMM_LEIDEN_MERGE_MAX_ITER", default = 1000L)),
  hmm_k_preferred = 5,
  hmm_n_iter = 500,
  hmm_n_start = 10,
  extended_anchor_older_validation_core = TRUE,
  extended_anchor_core = "GeoB25202_R2",
  xrf_half_width_cm = 0.5,
  q_threshold = 0.10,
  core_proxy_panel = c("ba_ti_z", "ca_ti_z", "zr_rb_z", "ti_al_z", "mn_fe_z"),
  oap_completeness_threshold = 0.30,
  exclude_deep_sediment_taxa = env_flag("WGCNA_HMM_EXCLUDE_DEEP_SEDIMENT_TAXA", default = FALSE),
  extra_excluded_phyla = env_csv("WGCNA_HMM_EXTRA_EXCLUDED_PHYLA", default = character()),
  deep_sediment_display_category = "Diagenetic",
  deep_sediment_signal_sources = c("diagenetic_in_situ", "benthic_resident"),
  deep_sediment_confidence_min = 0.85
)

TEA <- list(
  oap_modules = list(
    M00155 = list(class = "O2", eo_mv = 820, disc = 1.00),
    M00154 = list(class = "O2", eo_mv = 820, disc = 1.00),
    M00416 = list(class = "O2", eo_mv = 820, disc = 1.00),
    M00417 = list(class = "O2", eo_mv = 820, disc = 1.00),
    M00156 = list(class = "O2", eo_mv = 820, disc = 0.50),
    M00153 = list(class = "O2", eo_mv = 820, disc = 0.25),
    M00529 = list(class = "NO3", eo_mv = 740, disc = 1.00),
    M00530 = list(class = "NO3", eo_mv = 360, disc = 1.00),
    M00973 = list(class = "ANX", eo_mv = 350, disc = 1.00),
    M00596 = list(class = "SO4", eo_mv = -217, disc = 1.00),
    M00567 = list(class = "CH4", eo_mv = -244, disc = 1.00),
    M00357 = list(class = "CH4", eo_mv = -244, disc = 1.00),
    M00356 = list(class = "CH4", eo_mv = -244, disc = 1.00),
    M00563 = list(class = "CH4", eo_mv = -244, disc = 1.00)
  ),
  target_kos = c(
    "K00399", "K00401", "K00402",
    "K14080",
    "K00370",
    "K02305",
    "K00376",
    "K00425",
    "K02274", "K02276",
    "K00394", "K00395",
    "K00958",
    "K11180", "K11181"
  )
)

log_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

ensure_dirs <- function(paths) {
  for (d in paths) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_session_info <- function(path) {
  si <- utils::capture.output(sessionInfo())
  writeLines(si, con = path)
}

write_run_metadata <- function(path, script_name, extra = list()) {
  md <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script = script_name,
    repo_root = REPO_ROOT,
    seed = PARAMS$seed
  )
  if (length(extra) > 0) {
    for (nm in names(extra)) {
      val <- extra[[nm]]
      if (length(val) > 1) val <- paste(val, collapse = ";")
      md[[nm]] <- as.character(val)
    }
  }
  fwrite(md, path, sep = "\t")
}

core_group <- function(x) sub("_R[12]$", "", x)

corewise_z_transform <- function(mat, core_ids, feature_names = colnames(mat), variant = "corewise_z") {
  x <- as.matrix(mat)
  storage.mode(x) <- "double"
  if (length(core_ids) != nrow(x)) {
    stop("core_ids length must match number of rows in matrix.")
  }
  if (is.null(feature_names)) {
    feature_names <- paste0("V", seq_len(ncol(x)))
  }
  colnames(x) <- feature_names

  cores <- unique(as.character(core_ids))
  out <- x
  param_rows <- vector("list", length(cores))

  for (i in seq_along(cores)) {
    co <- cores[i]
    idx <- which(core_ids == co)
    x_sub <- x[idx, , drop = FALSE]
    mu <- colMeans(x_sub, na.rm = TRUE)
    sdv <- apply(x_sub, 2, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    out[idx, ] <- sweep(sweep(x_sub, 2, mu, "-"), 2, sdv, "/")

    param_rows[[i]] <- data.table(
      core = co,
      feature = feature_names,
      mu = as.numeric(mu),
      sd = as.numeric(sdv),
      variant = variant
    )
  }

  list(
    scaled = out,
    params = rbindlist(param_rows, use.names = TRUE)
  )
}

compute_clr_train_centered <- function(count_mat_taxa_by_samples, train_sample_ids, pseudocount = 0.5) {
  log_mat <- log(count_mat_taxa_by_samples + pseudocount)
  train_ids <- intersect(train_sample_ids, colnames(log_mat))
  if (length(train_ids) == 0) stop("No training sample IDs found in count matrix.")
  train_row_means <- rowMeans(log_mat[, train_ids, drop = FALSE])
  clr <- sweep(log_mat, 1, train_row_means, "-")
  keep <- apply(clr, 1, function(v) stats::var(v, na.rm = TRUE) > 0)
  clr <- clr[keep, , drop = FALSE]
  list(
    clr = t(clr),
    train_row_means = train_row_means[keep],
    kept_taxa = rownames(clr)
  )
}
