#!/usr/bin/env Rscript

# =============================================================================
# Eukaryote distribution vs XRF proxies
#
# Goal:
# - Community-level associations (Hellinger/PCA + PERMANOVA when available)
# - Taxon-wise associations (presence and count models)
# - Primary rank: genus; family as stability check; species exploratory
#
# Outputs are kept in memory only (no files written).
# =============================================================================

setwd("/maps/projects/caeg/people/ngm902/apps/repos/rocs")
repo_root <- getwd()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

pkg_ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)

analysis_notes <- c(
  "Primary analyses restricted to interval-quality XRF matches.",
  "Community models use one proxy at a time to reduce multicollinearity inflation.",
  "Taxon-wise models are exploratory and compositional caveats apply to read counts."
)

settings <- list(
  repo_root = repo_root,
  primary_rank = "genus",
  secondary_rank = "family",
  exploratory_rank = "species",
  preferred_proxies = c("ba_ti_z", "ca_ti_z", "zr_rb_z", "ti_al_z", "mn_fe_z"),
  min_sample_reads = 100,
  min_sample_taxa = 5,
  min_taxon_prevalence = 0.05,
  min_taxon_total_reads = 100,
  max_taxa_models = 36,
  min_positive_for_count = 10,
  n_perm = 999
)

required_files <- c(
  "results/eukaryotes/taxonomy/genus_exploration.rds",
  "results/eukaryotes/taxonomy/family_exploration.rds",
  "results/eukaryotes/taxonomy/species_exploration.rds",
  "results/common/metadata_xrf_matched.tsv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n - ",
    paste(missing_files, collapse = "\n - "),
    call. = FALSE
  )
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

coalesce_chr <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  ifelse(is.na(a) | a == "", b, a)
}

first_present_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

extract_proxy_row <- function(model, proxy_name, model_kind = c("glm", "glmer", "glmmTMB", "gam")) {
  model_kind <- match.arg(model_kind)
  out <- tibble(term = proxy_name, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p_value = NA_real_)

  if (model_kind %in% c("glm", "glmer", "gam")) {
    sm <- summary(model)
    cf <- sm$coefficients
    if (is.null(cf) && !is.null(sm$p.table)) cf <- sm$p.table
    if (!is.null(cf) && proxy_name %in% rownames(cf)) {
      row <- cf[proxy_name, , drop = FALSE]
      out$estimate <- as.numeric(row[1, 1])
      if (ncol(row) >= 2) out$std_error <- as.numeric(row[1, 2])
      if (ncol(row) >= 3) out$statistic <- as.numeric(row[1, 3])
      if (ncol(row) >= 4) out$p_value <- as.numeric(row[1, 4])
    }
    return(out)
  }

  if (model_kind == "glmmTMB") {
    sm <- summary(model)
    cf <- sm$coefficients$cond
    if (!is.null(cf) && proxy_name %in% rownames(cf)) {
      row <- cf[proxy_name, , drop = FALSE]
      out$estimate <- as.numeric(row[1, "Estimate"])
      out$std_error <- as.numeric(row[1, "Std. Error"])
      out$statistic <- as.numeric(row[1, "z value"])
      out$p_value <- as.numeric(row[1, "Pr(>|z|)"])
    }
    return(out)
  }

  out
}

choose_rank_table <- function(exploration_obj, rank_name) {
  candidates <- exploration_obj$data
  if ("dmg" %in% names(candidates) && nrow(candidates$dmg) > 0) {
    analysis_notes <<- c(
      analysis_notes,
      paste0("Using `", rank_name, "$data$dmg` as primary damage-supported table.")
    )
    return(candidates$dmg)
  }
  analysis_notes <<- c(
    analysis_notes,
    paste0("No non-empty `", rank_name, "$data$dmg`; using full table as fallback.")
  )
  candidates$full
}

read_rank_inputs <- function() {
  list(
    genus = readRDS("results/eukaryotes/taxonomy/genus_exploration.rds"),
    family = readRDS("results/eukaryotes/taxonomy/family_exploration.rds"),
    species = readRDS("results/eukaryotes/taxonomy/species_exploration.rds")
  )
}

metadata_xrf <- read_tsv("results/common/metadata_xrf_matched.tsv", show_col_types = FALSE) %>%
  mutate(
    across(any_of(c("y_bp", "depth_in_core_cm", "match_n_points", "match_min_abs_depth_cm")), safe_numeric),
    across(any_of(settings$preferred_proxies), safe_numeric),
    core = as.character(core),
    core_group = as.character(core_group),
    match_quality = as.character(match_quality)
  )

rank_inputs <- read_rank_inputs()

infer_taxon_column <- function(df, rank_name) {
  if (rank_name %in% names(df)) return(rank_name)
  if ("taxon" %in% names(df)) return("taxon")
  first_present_col(df, c("genus", "family", "species"))
}

preprocess_rank <- function(df_raw, rank_name, metadata_tbl, primary_only_interval = TRUE) {
  taxon_col <- infer_taxon_column(df_raw, rank_name)
  if (is.na(taxon_col)) stop("No taxon column found for rank: ", rank_name, call. = FALSE)

  keep_cols <- unique(c(
    "label", "core", "y_bp", "depth_in_core_cm", "is_dmg", "ecological_group", "n_reads",
    "kingdom", "phylum", "class", "order", "family", "genus", "species", taxon_col
  ))

  dat <- df_raw %>%
    mutate(across(any_of(c("n_reads", "y_bp", "depth_in_core_cm")), safe_numeric)) %>%
    select(any_of(keep_cols)) %>%
    rename(taxon = all_of(taxon_col)) %>%
    mutate(
      taxon = as.character(taxon),
      taxon = ifelse(is.na(taxon) | taxon == "", "Unclassified", taxon),
      label = as.character(label),
      core = as.character(core),
      ecological_group = if ("ecological_group" %in% names(.)) as.character(ecological_group) else NA_character_
    )

  meta_keep <- metadata_tbl %>%
    select(
      label, core_meta = core, y_bp_meta = y_bp, depth_meta = depth_in_core_cm,
      core_group, match_quality, match_n_points, match_min_abs_depth_cm,
      any_of(settings$preferred_proxies)
    )

  joined <- dat %>%
    left_join(meta_keep, by = "label") %>%
    mutate(
      core = coalesce_chr(core, core_meta),
      y_bp = ifelse(is.na(y_bp), y_bp_meta, y_bp),
      depth_in_core_cm = ifelse(is.na(depth_in_core_cm), depth_meta, depth_in_core_cm),
      core_group = coalesce_chr(core_group, core),
      n_reads = ifelse(is.na(n_reads), 0, n_reads)
    ) %>%
    select(-core_meta, -y_bp_meta, -depth_meta)

  if (!"is_dmg" %in% names(joined)) {
    analysis_notes <<- c(analysis_notes, paste0("Rank `", rank_name, "` has no `is_dmg`; limitation for damage filtering."))
  }

  if (primary_only_interval) {
    joined <- joined %>% filter(match_quality == "interval")
  }

  sample_info <- joined %>%
    group_by(label, core, core_group, y_bp, depth_in_core_cm, match_quality) %>%
    summarise(
      total_rank_reads = sum(n_reads, na.rm = TRUE),
      n_taxa_detected = n_distinct(taxon[n_reads > 0]),
      .groups = "drop"
    )

  min_reads <- settings$min_sample_reads
  min_taxa <- settings$min_sample_taxa

  if (nrow(sample_info %>% filter(total_rank_reads >= min_reads, n_taxa_detected >= min_taxa)) < 20) {
    min_taxa <- max(3, min_taxa - 2)
    analysis_notes <<- c(
      analysis_notes,
      paste0("Rank `", rank_name, "`: relaxed min_sample_taxa to ", min_taxa, " due to low retained sample count.")
    )
  }

  keep_samples <- sample_info %>%
    filter(total_rank_reads >= min_reads, n_taxa_detected >= min_taxa) %>%
    pull(label)

  dat_s <- joined %>% filter(label %in% keep_samples)

  taxon_summary <- dat_s %>%
    group_by(taxon) %>%
    summarise(
      total_reads = sum(n_reads, na.rm = TRUE),
      prevalence = mean(n_reads > 0, na.rm = TRUE),
      n_positive = sum(n_reads > 0, na.rm = TRUE),
      ecological_group = first(na.omit(ecological_group)),
      .groups = "drop"
    )

  taxon_keep <- taxon_summary %>%
    filter(
      total_reads >= settings$min_taxon_total_reads,
      prevalence >= settings$min_taxon_prevalence,
      !str_detect(tolower(taxon), "unclassified|uncultured|unknown")
    )

  if ("ecological_group" %in% names(taxon_keep)) {
    taxon_keep <- taxon_keep %>%
      filter(is.na(ecological_group) | ecological_group != "Ambiguous/non-marine")
  }

  if (nrow(taxon_keep) > settings$max_taxa_models) {
    taxon_keep <- taxon_keep %>%
      arrange(desc(prevalence), desc(total_reads)) %>%
      slice_head(n = settings$max_taxa_models)
    analysis_notes <<- c(
      analysis_notes,
      paste0("Rank `", rank_name, "`: truncated taxa to top ", settings$max_taxa_models, " for model stability/runtime.")
    )
  }

  dat_f <- dat_s %>% filter(taxon %in% taxon_keep$taxon)

  list(
    joined = joined,
    sample_info = sample_info,
    filtered = dat_f,
    taxon_summary = taxon_keep,
    thresholds = list(min_reads = min_reads, min_taxa = min_taxa)
  )
}

pick_available_proxies <- function(df, preferred) {
  present <- preferred[preferred %in% names(df)]
  present[vapply(present, function(p) sum(is.finite(df[[p]])) >= 20, logical(1))]
}

community_analysis <- function(dat_f, rank_name, proxies) {
  if (nrow(dat_f) == 0 || length(proxies) == 0) {
    return(list(notes = paste0("Rank `", rank_name, "`: insufficient data for community analysis.")))
  }

  first_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    x[[1]]
  }

  sample_dat <- dat_f %>%
    group_by(label, core, core_group, y_bp, depth_in_core_cm) %>%
    summarise(
      total_rank_reads = sum(n_reads, na.rm = TRUE),
      across(any_of(proxies), first_finite),
      .groups = "drop"
    ) %>%
    distinct()

  mat <- dat_f %>%
    group_by(label, taxon) %>%
    summarise(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = taxon, values_from = n_reads, values_fill = 0)

  rn <- mat$label
  X <- as.matrix(mat[, setdiff(names(mat), "label"), drop = FALSE])
  rownames(X) <- rn
  rs <- rowSums(X)
  X <- X[rs > 0, , drop = FALSE]

  sample_dat <- sample_dat %>% filter(label %in% rownames(X))
  sample_dat <- sample_dat[match(rownames(X), sample_dat$label), , drop = FALSE]

  p <- X / rowSums(X)
  p[!is.finite(p)] <- 0
  hell <- sqrt(p)

  pca <- prcomp(hell, center = TRUE, scale. = FALSE)
  ve <- (pca$sdev^2) / sum(pca$sdev^2)

  axis_df <- as_tibble(pca$x[, seq_len(min(5, ncol(pca$x))), drop = FALSE], rownames = "label") %>%
    left_join(sample_dat, by = "label")

  proxy_axis <- map_dfr(proxies, function(px) {
    map_dfr(seq_len(min(3, ncol(pca$x))), function(ax) {
      x <- axis_df[[paste0("PC", ax)]]
      y <- axis_df[[px]]
      keep <- is.finite(x) & is.finite(y)
      ct <- if (sum(keep) >= 5) suppressWarnings(cor.test(x[keep], y[keep], method = "spearman")) else NULL
      tibble(
        rank = rank_name,
        proxy = px,
        axis = paste0("PC", ax),
        rho = if (!is.null(ct)) unname(ct$estimate) else NA_real_,
        p_value = if (!is.null(ct)) ct$p.value else NA_real_,
        n = sum(keep)
      )
    })
  }) %>% mutate(q_value = p.adjust(p_value, method = "BH"))

  adonis_tbl <- tibble()
  dispersion_tbl <- tibble()

  if (pkg_ok("vegan")) {
    D <- vegan::vegdist(hell, method = "euclidean")
    adonis_tbl <- map_dfr(proxies, function(px) {
      dd <- sample_dat %>%
        mutate(
          proxy = .data[[px]],
          age = y_bp,
          core_group = as.factor(core_group)
        ) %>%
        select(proxy, age, core_group)

      keep <- is.finite(dd$proxy) & is.finite(dd$age) & !is.na(dd$core_group)
      if (sum(keep) < 20) {
        return(tibble(rank = rank_name, proxy = px, n = sum(keep), F = NA_real_, R2 = NA_real_, p_value = NA_real_, model = "adonis2"))
      }

      Dk <- as.dist(as.matrix(D)[keep, keep])
      ddk <- dd[keep, , drop = FALSE]
      strata <- if (length(unique(ddk$core_group)) > 1) ddk$core_group else NULL

      fit <- tryCatch(
        vegan::adonis2(Dk ~ proxy + age + core_group, data = ddk, permutations = settings$n_perm, strata = strata),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(tibble(rank = rank_name, proxy = px, n = nrow(ddk), F = NA_real_, R2 = NA_real_, p_value = NA_real_, model = "adonis2_failed"))
      }

      tt <- as.data.frame(fit)
      rr <- rownames(tt)
      prox_row <- which(rr == "proxy")
      if (length(prox_row) != 1) {
        return(tibble(rank = rank_name, proxy = px, n = nrow(ddk), F = NA_real_, R2 = NA_real_, p_value = NA_real_, model = "adonis2_no_proxy_row"))
      }

      tibble(
        rank = rank_name,
        proxy = px,
        n = nrow(ddk),
        F = tt$F[prox_row],
        R2 = tt$R2[prox_row],
        p_value = tt$`Pr(>F)`[prox_row],
        model = "adonis2"
      )
    }) %>%
      mutate(q_value = p.adjust(p_value, method = "BH"))

    dispersion_tbl <- map_dfr(proxies, function(px) {
      dd <- sample_dat %>% mutate(proxy = .data[[px]])
      keep <- is.finite(dd$proxy)
      if (sum(keep) < 20) {
        return(tibble(rank = rank_name, proxy = px, n = sum(keep), F = NA_real_, p_value = NA_real_, method = "betadisper"))
      }
      dd2 <- dd[keep, , drop = FALSE]
      dd2$proxy_bin <- dplyr::ntile(dd2$proxy, 3)
      if (dplyr::n_distinct(dd2$proxy_bin) < 2) {
        return(tibble(rank = rank_name, proxy = px, n = nrow(dd2), F = NA_real_, p_value = NA_real_, method = "betadisper"))
      }
      D2 <- as.dist(as.matrix(D)[keep, keep])
      bd <- tryCatch(vegan::betadisper(D2, group = factor(dd2$proxy_bin)), error = function(e) NULL)
      if (is.null(bd)) {
        return(tibble(rank = rank_name, proxy = px, n = nrow(dd2), F = NA_real_, p_value = NA_real_, method = "betadisper_failed"))
      }
      av <- anova(bd)
      tibble(rank = rank_name, proxy = px, n = nrow(dd2), F = av$`F value`[1], p_value = av$`Pr(>F)`[1], method = "betadisper")
    })
  } else {
    analysis_notes <<- c(analysis_notes, "Package `vegan` not available; PERMANOVA/dispersion checks skipped.")
  }

  proxy_depth_age <- map_dfr(proxies, function(px) {
    dd <- sample_dat %>% mutate(proxy = .data[[px]])
    keep_age <- is.finite(dd$proxy) & is.finite(dd$y_bp)
    keep_dep <- is.finite(dd$proxy) & is.finite(dd$depth_in_core_cm)
    ca <- if (sum(keep_age) >= 5) suppressWarnings(cor.test(dd$proxy[keep_age], dd$y_bp[keep_age], method = "spearman")) else NULL
    cd <- if (sum(keep_dep) >= 5) suppressWarnings(cor.test(dd$proxy[keep_dep], dd$depth_in_core_cm[keep_dep], method = "spearman")) else NULL
    tibble(
      rank = rank_name,
      proxy = px,
      rho_age = if (!is.null(ca)) unname(ca$estimate) else NA_real_,
      p_age = if (!is.null(ca)) ca$p.value else NA_real_,
      rho_depth = if (!is.null(cd)) unname(cd$estimate) else NA_real_,
      p_depth = if (!is.null(cd)) cd$p.value else NA_real_,
      n_age = sum(keep_age),
      n_depth = sum(keep_dep)
    )
  })

  list(
    rank = rank_name,
    n_samples = nrow(sample_dat),
    n_taxa = ncol(hell),
    pca_variance_explained = tibble(axis = paste0("PC", seq_along(ve)), variance_explained = ve),
    proxy_axis_assoc = proxy_axis,
    permanova = adonis_tbl,
    dispersion = dispersion_tbl,
    proxy_age_depth_correlation = proxy_depth_age,
    ordination = list(pca = pca, hellinger = hell)
  )
}

fit_taxon_models <- function(dat_f, rank_name, proxies, exploratory = FALSE, max_taxa_override = NULL) {
  if (nrow(dat_f) == 0 || length(proxies) == 0) {
    return(list(
      presence = tibble(),
      count = tibble(),
      notes = paste0("Rank `", rank_name, "`: insufficient data for taxon-wise models.")
    ))
  }

  first_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    x[[1]]
  }

  taxa <- dat_f %>%
    group_by(taxon) %>%
    summarise(
      n_positive = n_distinct(label[n_reads > 0]),
      total_reads = sum(n_reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_positive >= ifelse(exploratory, 8, settings$min_positive_for_count)) %>%
    arrange(desc(n_positive), desc(total_reads)) %>%
    pull(taxon)

  if (length(taxa) == 0) {
    return(list(
      presence = tibble(),
      count = tibble(),
      notes = paste0("Rank `", rank_name, "`: no taxa passed positivity threshold.")
    ))
  }

  max_taxa_use <- ifelse(is.null(max_taxa_override), settings$max_taxa_models, max_taxa_override)
  if (length(taxa) > max_taxa_use) {
    taxa <- taxa[seq_len(max_taxa_use)]
  }

  sample_cov <- dat_f %>%
    group_by(label) %>%
    summarise(
      core_group = as.factor(first(core_group)),
      y_bp = first_finite(y_bp),
      total_rank_reads = sum(n_reads, na.rm = TRUE),
      across(any_of(proxies), first_finite),
      .groups = "drop"
    )

  counts <- dat_f %>%
    filter(taxon %in% taxa) %>%
    group_by(label, taxon) %>%
    summarise(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop")

  base_df <- tidyr::expand_grid(label = sample_cov$label, taxon = taxa) %>%
    left_join(counts, by = c("label", "taxon")) %>%
    left_join(sample_cov, by = "label") %>%
    mutate(
      n_reads = replace_na(n_reads, 0),
      core_group = as.factor(core_group),
      y_bp = safe_numeric(y_bp),
      log_total_offset = log(pmax(total_rank_reads, 1)),
      pa = as.integer(n_reads > 0)
    )

  presence_rows <- list()
  count_rows <- list()

  for (tx in taxa) {
    dtx <- base_df %>% filter(taxon == tx)

    for (px in proxies) {
      dmod <- dtx %>%
        transmute(
          pa = pa,
          n_reads = n_reads,
          proxy = .data[[px]],
          y_bp = y_bp,
          core_group = core_group,
          log_total_offset = log_total_offset,
          total_rank_reads = total_rank_reads
        ) %>%
        filter(is.finite(proxy), is.finite(y_bp), is.finite(log_total_offset), !is.na(core_group))

      if (nrow(dmod) < 20 || dplyr::n_distinct(dmod$pa) < 2) {
        next
      }

      # Presence model
      pres_fit <- NULL
      pres_method <- NA_character_
      pres_note <- NA_character_

      try_mixed <- pkg_ok("lme4") && !exploratory && rank_name %in% c("genus", "family") && dplyr::n_distinct(dmod$core_group) > 1
      if (try_mixed) {
        pres_fit <- tryCatch(
          suppressWarnings(
            lme4::glmer(
              pa ~ proxy + splines::ns(y_bp, 3) + (1 | core_group),
              data = dmod,
              family = binomial(),
              control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e4))
            )
          ),
          error = function(e) NULL
        )
        pres_method <- "glmer_binomial"
      }

      if (is.null(pres_fit)) {
        pres_fit <- tryCatch(
          glm(pa ~ proxy + splines::ns(y_bp, 3) + core_group, data = dmod, family = binomial()),
          error = function(e) NULL
        )
        pres_method <- "glm_binomial"
      }

      if (!is.null(pres_fit)) {
        kind <- if (pres_method == "glmer_binomial") "glmer" else "glm"
        ext <- extract_proxy_row(pres_fit, "proxy", model_kind = kind)
        singular <- if (pres_method == "glmer_binomial" && pkg_ok("lme4")) {
          isTRUE(lme4::isSingular(pres_fit, tol = 1e-04))
        } else {
          NA
        }
        if (isTRUE(singular)) {
          pres_note <- "singular_fit_refit_glm"
          pres_fit2 <- tryCatch(
            glm(pa ~ proxy + splines::ns(y_bp, 3) + core_group, data = dmod, family = binomial()),
            error = function(e) NULL
          )
          if (!is.null(pres_fit2)) {
            pres_fit <- pres_fit2
            pres_method <- "glm_binomial"
            kind <- "glm"
            ext <- extract_proxy_row(pres_fit, "proxy", model_kind = kind)
          }
        }

        presence_rows[[length(presence_rows) + 1]] <- tibble(
          rank = rank_name,
          taxon = tx,
          proxy = px,
          method = pres_method,
          n = nrow(dmod),
          n_positive = sum(dmod$pa),
          estimate = ext$estimate,
          std_error = ext$std_error,
          statistic = ext$statistic,
          p_value = ext$p_value,
          convergence_note = pres_note
        )
      }

      # Count model
      cnt_fit <- NULL
      cnt_method <- NA_character_
      cnt_note <- NA_character_

      if (pkg_ok("MASS")) {
        cnt_fit <- tryCatch(
          MASS::glm.nb(n_reads ~ proxy + splines::ns(y_bp, 3) + core_group + offset(log_total_offset), data = dmod),
          error = function(e) NULL
        )
        cnt_method <- "glm_nb"
      }

      if (is.null(cnt_fit) && pkg_ok("glmmTMB") && !exploratory && rank_name %in% c("genus", "family") && dplyr::n_distinct(dmod$core_group) > 1) {
        cnt_fit <- tryCatch(
          suppressWarnings(
            glmmTMB::glmmTMB(
              n_reads ~ proxy + splines::ns(y_bp, 3) + (1 | core_group) + offset(log_total_offset),
              data = dmod,
              family = glmmTMB::nbinom2()
            )
          ),
          error = function(e) NULL
        )
        cnt_method <- "glmmTMB_nbinom2"
      }

      if (is.null(cnt_fit)) {
        cnt_fit <- tryCatch(
          glm(n_reads ~ proxy + splines::ns(y_bp, 3) + core_group + offset(log_total_offset), data = dmod, family = quasipoisson()),
          error = function(e) NULL
        )
        cnt_method <- "glm_quasipoisson"
      }

      if (!is.null(cnt_fit)) {
        kind <- if (cnt_method == "glmmTMB_nbinom2") "glmmTMB" else "glm"
        ext <- extract_proxy_row(cnt_fit, "proxy", model_kind = kind)
        if (cnt_method == "glmmTMB_nbinom2") {
          hess_ok <- tryCatch(isTRUE(cnt_fit$sdr$pdHess), error = function(e) NA)
          if (isFALSE(hess_ok)) cnt_note <- "non_pd_hessian"
        }

        count_rows[[length(count_rows) + 1]] <- tibble(
          rank = rank_name,
          taxon = tx,
          proxy = px,
          method = cnt_method,
          n = nrow(dmod),
          n_positive = sum(dmod$pa),
          estimate = ext$estimate,
          std_error = ext$std_error,
          statistic = ext$statistic,
          p_value = ext$p_value,
          convergence_note = cnt_note
        )
      }
    }
  }

  presence_tbl <- bind_rows(presence_rows)
  count_tbl <- bind_rows(count_rows)

  if (nrow(presence_tbl) > 0) {
    presence_tbl <- presence_tbl %>%
      group_by(rank, proxy) %>%
      mutate(q_value = p.adjust(p_value, method = "BH")) %>%
      ungroup()
  }

  if (nrow(count_tbl) > 0) {
    count_tbl <- count_tbl %>%
      group_by(rank, proxy) %>%
      mutate(q_value = p.adjust(p_value, method = "BH")) %>%
      ungroup()
  }

  presence_comp <- if (nrow(presence_tbl) > 0) {
    presence_tbl %>% select(rank, taxon, proxy, estimate_pa = estimate, q_pa = q_value)
  } else {
    tibble(rank = character(), taxon = character(), proxy = character(), estimate_pa = numeric(), q_pa = numeric())
  }

  count_comp <- if (nrow(count_tbl) > 0) {
    count_tbl %>% select(rank, taxon, proxy, estimate_count = estimate, q_count = q_value)
  } else {
    tibble(rank = character(), taxon = character(), proxy = character(), estimate_count = numeric(), q_count = numeric())
  }

  robustness <- full_join(presence_comp, count_comp, by = c("rank", "taxon", "proxy")) %>%
    mutate(
      sign_consistent = sign(estimate_pa) == sign(estimate_count),
      strong_both = !is.na(q_pa) & !is.na(q_count) & q_pa < 0.10 & q_count < 0.10 & sign_consistent
    )

  list(
    presence = presence_tbl,
    count = count_tbl,
    robustness = robustness,
    tested_taxa = taxa
  )
}

run_rank_distribution <- function(rank_name, exploration_obj, metadata_tbl, exploratory = FALSE) {
  dat_primary <- choose_rank_table(exploration_obj, rank_name)
  pre <- preprocess_rank(dat_primary, rank_name, metadata_tbl, primary_only_interval = TRUE)
  proxies <- pick_available_proxies(pre$filtered, settings$preferred_proxies)

  if (length(proxies) == 0) {
    analysis_notes <<- c(analysis_notes, paste0("Rank `", rank_name, "`: no preferred z-scored proxies available with sufficient data."))
  }

  proxy_corr <- tibble()
  if (length(proxies) >= 2) {
    corr_mat <- suppressWarnings(cor(pre$filtered[, proxies, drop = FALSE], use = "pairwise.complete.obs", method = "spearman"))
    proxy_corr <- as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE) %>%
      as_tibble() %>%
      rename(proxy1 = Var1, proxy2 = Var2, rho = Freq) %>%
      filter(proxy1 != proxy2)
  }

  community <- community_analysis(pre$filtered, rank_name, proxies)
  max_taxa_rank <- dplyr::case_when(
    rank_name == "genus" ~ 36,
    rank_name == "family" ~ 24,
    rank_name == "species" ~ 12,
    TRUE ~ settings$max_taxa_models
  )

  taxon_models <- fit_taxon_models(pre$filtered, rank_name, proxies, exploratory = exploratory, max_taxa_override = max_taxa_rank)

  # sensitivity: include fallback_nearest matches
  pre_sens <- preprocess_rank(dat_primary, rank_name, metadata_tbl, primary_only_interval = FALSE)
  sens_counts <- pre_sens$joined %>%
    count(match_quality, name = "n_rows")

  list(
    rank = rank_name,
    settings = pre$thresholds,
    proxies_used = proxies,
    input_summary = list(
      n_rows_primary_joined = nrow(pre$joined),
      n_rows_filtered = nrow(pre$filtered),
      n_samples_filtered = n_distinct(pre$filtered$label),
      n_taxa_filtered = n_distinct(pre$filtered$taxon),
      match_quality_sensitivity_counts = sens_counts
    ),
    preprocessed = pre,
    proxy_collinearity = proxy_corr,
    community = community,
    taxon_models = taxon_models,
    sensitivity = list(include_fallback = pre_sens$sample_info)
  )
}

optional_hmsc <- function(rank_result) {
  if (!pkg_ok("Hmsc")) {
    return(list(ran = FALSE, note = "Package `Hmsc` not available; skipped."))
  }

  dat <- rank_result$preprocessed$filtered
  if (nrow(dat) == 0) {
    return(list(ran = FALSE, note = "No data after filtering; skipped."))
  }

  n_samples <- n_distinct(dat$label)
  n_taxa <- n_distinct(dat$taxon)
  if (n_samples > 250 || n_taxa > 80) {
    return(list(ran = FALSE, note = paste0("Data too large for lightweight Hmsc block (", n_samples, " samples, ", n_taxa, " taxa).")))
  }

  if (!("ba_ti_z" %in% rank_result$proxies_used)) {
    return(list(ran = FALSE, note = "`ba_ti_z` not available for Hmsc test; skipped."))
  }

  Y <- dat %>%
    group_by(label, taxon) %>%
    summarise(pa = as.integer(sum(n_reads, na.rm = TRUE) > 0), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = pa, values_fill = 0)

  XData <- dat %>%
    group_by(label, core_group, y_bp, ba_ti_z) %>%
    summarise(.groups = "drop") %>%
    distinct() %>%
    filter(label %in% Y$label)

  Ym <- as.matrix(Y[, setdiff(names(Y), "label"), drop = FALSE])
  rownames(Ym) <- Y$label
  XData <- XData[match(rownames(Ym), XData$label), , drop = FALSE]
  keep <- is.finite(XData$ba_ti_z) & is.finite(XData$y_bp)

  Ymk <- Ym[keep, , drop = FALSE]
  Xk <- XData[keep, , drop = FALSE]

  if (nrow(Ymk) < 30) {
    return(list(ran = FALSE, note = "Too few complete samples for Hmsc block after NA filtering."))
  }

  out <- tryCatch({
    XFormula <- ~ ba_ti_z + splines::ns(y_bp, 3)
    studyDesign <- data.frame(sample = Xk$label)
    rL <- Hmsc::HmscRandomLevel(units = studyDesign$sample)
    m <- Hmsc::Hmsc(
      Y = Ymk,
      XData = Xk,
      XFormula = XFormula,
      distr = "probit",
      studyDesign = studyDesign,
      ranLevels = list(sample = rL)
    )
    mfit <- Hmsc::sampleMcmc(m, samples = 40, thin = 2, transient = 20, nChains = 1, verbose = 0)
    list(ran = TRUE, note = "Hmsc lightweight test completed.", model = mfit)
  }, error = function(e) {
    list(ran = FALSE, note = paste("Hmsc block failed:", conditionMessage(e)))
  })

  out
}

genus_distribution <- run_rank_distribution("genus", rank_inputs$genus, metadata_xrf, exploratory = FALSE)
family_distribution <- run_rank_distribution("family", rank_inputs$family, metadata_xrf, exploratory = FALSE)
species_distribution <- run_rank_distribution("species", rank_inputs$species, metadata_xrf, exploratory = TRUE)

hmsc_results <- optional_hmsc(genus_distribution)

diagnostics <- list(
  package_availability = c(
    vegan = pkg_ok("vegan"),
    lme4 = pkg_ok("lme4"),
    glmmTMB = pkg_ok("glmmTMB"),
    MASS = pkg_ok("MASS"),
    mgcv = pkg_ok("mgcv"),
    Hmsc = pkg_ok("Hmsc")
  ),
  metadata_match_quality = metadata_xrf %>% count(match_quality, name = "n")
)

distribution_results <- list(
  settings = settings,
  notes = unique(analysis_notes),
  input_summaries = list(
    metadata_rows = nrow(metadata_xrf),
    metadata_samples = n_distinct(metadata_xrf$label),
    genus_rows_dmg = nrow(rank_inputs$genus$data$dmg),
    family_rows_dmg = nrow(rank_inputs$family$data$dmg),
    species_rows_dmg = nrow(rank_inputs$species$data$dmg)
  ),
  preprocessed_data = list(
    genus = genus_distribution$preprocessed,
    family = family_distribution$preprocessed,
    species = species_distribution$preprocessed
  ),
  community_results = list(
    genus = genus_distribution$community,
    family = family_distribution$community,
    species = species_distribution$community
  ),
  taxon_model_results = list(
    genus = genus_distribution$taxon_models,
    family = family_distribution$taxon_models,
    species = species_distribution$taxon_models
  ),
  optional_hmsc_results = hmsc_results,
  diagnostics = diagnostics,
  sensitivity = list(
    genus = genus_distribution$sensitivity,
    family = family_distribution$sensitivity,
    species = species_distribution$sensitivity
  )
)

message("Created in-memory objects:")
message(" - distribution_results")
message(" - genus_distribution")
message(" - family_distribution")
message(" - species_distribution")

message("Methods summary:")
message(" - Genus community samples/taxa: ", genus_distribution$community$n_samples, "/", genus_distribution$community$n_taxa)
message(" - Family community samples/taxa: ", family_distribution$community$n_samples, "/", family_distribution$community$n_taxa)
message(" - Species community samples/taxa: ", species_distribution$community$n_samples, "/", species_distribution$community$n_taxa)
message(" - HMSC: ", hmsc_results$note)
