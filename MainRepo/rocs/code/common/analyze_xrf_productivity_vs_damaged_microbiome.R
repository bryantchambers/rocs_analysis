#!/usr/bin/env Rscript

# Reworked XRF productivity-state vs damaged-only microbiome analysis.
#
# Guardrails for manuscript-safe interpretation:
# - Damaged-only filtering changes the ecological target (preserved/authentic-leaning
#   fraction), and may not represent whole original communities.
# - Taxonomic community profiles are compositional, sparse, and partly unresolved.
# - XRF-derived axes are indirect environmental summaries, not direct mechanisms.
# - Associations are observational and cannot prove causal productivity control.
# - Core identity and age/depth structure remain major confounders.
# - Some taxa may reflect persistence/preservation/contamination risk rather than
#   clear ecological function.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(splines)
  library(forcats)
  library(vegan)
})

# ------------------------------------------------------------------
# Paths and run contract
# ------------------------------------------------------------------
damage_path <- file.path(
  "results", "microbial", "damage", "damage-classification-depositional",
  "dmg-summary-ssp-damage-classification-depositional.tsv.gz"
)
xrf_path <- file.path("data", "combined_xrf_geochemistry_curated.csv")
meta_path <- file.path("data", "metadata_v5.tsv")
foram_path <- file.path("data", "combined_foraminifera_geochem.tsv")
sst_path <- file.path("data", "combined_sst_proxies_separate_columns.tsv")
geob_path <- file.path("data", "geob25202_clean_proxies.tsv")

out_dir <- file.path("results", "common", "proxies", "final")
plot_dir <- file.path("plots", "proxies", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(damage_path, xrf_path, meta_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required inputs (run from repo root): ", paste(missing, collapse = ", "))
}

params <- list(
  damaged_filter = "is_dmg == 'Damaged'",
  count_field_forced = "abundance",
  primary_domains = "Bacteria + Archaea",
  virus_sensitivity_run = TRUE,
  rank_priority = "genus_with_family_fallback",
  min_layer_total_abundance = 20,
  min_taxon_prevalence_fraction = 0.10,
  min_taxon_prevalence_layers = 8,
  clr_pseudocount_rule = "half_min_positive_count",
  trend_spline_df = 4,
  permutations = 999,
  scenario_split_rule = "axis_quartiles (Q1 low, Q4 high)",
  scenario_min_n = 12,
  random_seed = 13
)

set.seed(params$random_seed)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
normalize_core <- function(x) {
  str_replace(as.character(x), "_R[0-9]+$", "")
}

safe_ratio <- function(num, den) {
  ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
}

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

safe_spearman <- function(x, y, min_n = 8) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_n) return(list(n = sum(ok), estimate = NA_real_, p = NA_real_))
  x2 <- x[ok]; y2 <- y[ok]
  if (dplyr::n_distinct(x2) < 3 || dplyr::n_distinct(y2) < 3) {
    return(list(n = length(x2), estimate = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x2, y2, method = "spearman", exact = FALSE))
  list(n = length(x2), estimate = unname(ct$estimate), p = ct$p.value)
}

layer_intervals <- function(depths) {
  d <- sort(unique(depths[is.finite(depths)]))
  if (length(d) == 0) {
    return(tibble(depth_in_core_cm = numeric(), lower = numeric(), upper = numeric()))
  }
  if (length(d) == 1) {
    return(tibble(depth_in_core_cm = d, lower = -Inf, upper = Inf))
  }
  mids <- (d[-1] + d[-length(d)]) / 2
  tibble(depth_in_core_cm = d, lower = c(-Inf, mids), upper = c(mids, Inf))
}

aggregate_dense_to_layers <- function(layer_df, dense_df, value_cols) {
  intervals <- layer_intervals(layer_df$depth_in_core_cm)
  if (nrow(intervals) == 0 || nrow(dense_df) == 0) {
    out <- intervals %>% select(depth_in_core_cm)
    out$n_dense_rows <- integer(nrow(out))
    for (nm in value_cols) out[[nm]] <- NA_real_
    return(out)
  }
  out <- intervals %>% select(depth_in_core_cm)
  out$n_dense_rows <- NA_integer_
  for (nm in value_cols) out[[nm]] <- NA_real_
  for (i in seq_len(nrow(intervals))) {
    lo <- intervals$lower[i]
    hi <- intervals$upper[i]
    idx <- which(is.finite(dense_df$depth_in_core_cm) & dense_df$depth_in_core_cm > lo & dense_df$depth_in_core_cm <= hi)
    out$n_dense_rows[i] <- length(idx)
    if (length(idx) == 0) next
    for (nm in value_cols) {
      vv <- dense_df[[nm]][idx]
      vv <- vv[is.finite(vv)]
      if (length(vv) > 0) out[[nm]][i] <- median(vv)
    }
  }
  out
}

choose_trend <- function(df) {
  if ("y_bp" %in% names(df)) {
    ok <- is.finite(df$y_bp)
    if (sum(ok) >= 10 && dplyr::n_distinct(df$y_bp[ok]) >= 5) return("y_bp")
  }
  if ("depth_in_core_cm" %in% names(df)) {
    ok <- is.finite(df$depth_in_core_cm)
    if (sum(ok) >= 10 && dplyr::n_distinct(df$depth_in_core_cm[ok]) >= 5) return("depth_in_core_cm")
  }
  NA_character_
}

shannon_from_counts <- function(v) {
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0 || sum(v) <= 0) return(NA_real_)
  p <- v / sum(v)
  -sum(p * log(p))
}

clr_transform <- function(mat_counts, pseudocount = NULL) {
  mat <- as.matrix(mat_counts)
  storage.mode(mat) <- "double"
  min_pos <- suppressWarnings(min(mat[is.finite(mat) & mat > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1
  if (is.null(pseudocount)) pseudocount <- min_pos / 2
  mat_pc <- mat + pseudocount
  mat_comp <- mat_pc / rowSums(mat_pc)
  lmat <- log(mat_comp)
  clr <- lmat - rowMeans(lmat)
  list(clr = clr, pseudocount = pseudocount, min_positive = min_pos)
}

build_axis <- function(df, vars, axis_name, orientation_priority) {
  vars <- vars[vars %in% names(df)]
  zmat <- as.data.frame(lapply(df[vars], safe_scale))
  names(zmat) <- vars
  usable <- vars[vapply(vars, function(v) {
    x <- zmat[[v]]
    sum(is.finite(x)) >= 12 && dplyr::n_distinct(x[is.finite(x)]) >= 4
  }, logical(1))]

  method <- "z_mean"
  axis_score <- rep(NA_real_, nrow(df))
  loadings <- tibble(
    axis = axis_name,
    variable = usable,
    method = "z_mean",
    raw_weight = 1,
    oriented_weight = 1
  )

  if (length(usable) >= 2) {
    zz <- zmat[, usable, drop = FALSE]
    cc <- complete.cases(zz)
    if (sum(cc) >= max(15, length(usable) + 3)) {
      p <- tryCatch(prcomp(zz[cc, , drop = FALSE], center = FALSE, scale. = FALSE), error = function(e) NULL)
      if (!is.null(p) && ncol(p$x) >= 1) {
        method <- "pca_pc1"
        axis_score[cc] <- p$x[, 1]
        loadings <- tibble(
          axis = axis_name,
          variable = usable,
          method = method,
          raw_weight = as.numeric(p$rotation[usable, 1]),
          oriented_weight = as.numeric(p$rotation[usable, 1])
        )
      }
    }
  }

  if (method == "z_mean") {
    axis_score <- rowMeans(zmat[, usable, drop = FALSE], na.rm = TRUE)
    axis_score[!is.finite(axis_score)] <- NA_real_
  }

  anchor <- orientation_priority[orientation_priority %in% usable][1]
  if (is.na(anchor) && length(usable) > 0) anchor <- usable[1]
  rho_before <- NA_real_
  if (!is.na(anchor)) {
    rho_before <- suppressWarnings(cor(axis_score, df[[anchor]], method = "spearman", use = "pairwise.complete.obs"))
    if (is.finite(rho_before) && rho_before < 0) {
      axis_score <- -axis_score
      loadings$oriented_weight <- -loadings$oriented_weight
    }
  }

  list(score = axis_score, method = method, anchor = anchor, rho_before = rho_before,
       usable = usable, loadings = loadings)
}

pick_taxon_rank <- function(df, rank_candidates = c("genus", "family")) {
  rank_candidates <- rank_candidates[rank_candidates %in% names(df)]
  if (length(rank_candidates) == 0) stop("Neither genus nor family column available.")
  rank_eval <- map_dfr(rank_candidates, function(rk) {
    tx <- as.character(df[[rk]])
    tx_ok <- !is.na(tx) & tx != "" & !str_detect(tolower(tx), "^unclassified|^unknown|^uncultured")
    tmp <- df %>%
      mutate(taxon = ifelse(tx_ok, tx, NA_character_)) %>%
      filter(!is.na(taxon)) %>%
      group_by(layer_id, taxon) %>%
      summarise(ab = sum(abundance, na.rm = TRUE), .groups = "drop")
    layer_taxa <- if (nrow(tmp) > 0) tmp %>% count(layer_id, name = "n_taxa") else tibble(layer_id = character(), n_taxa = integer())
    wide <- if (nrow(tmp) > 0) {
      tmp %>% pivot_wider(names_from = taxon, values_from = ab, values_fill = 0)
    } else {
      tibble(layer_id = unique(df$layer_id))
    }
    sparsity <- NA_real_
    if (ncol(wide) > 1) {
      m <- as.matrix(wide[, -1, drop = FALSE])
      sparsity <- mean(m == 0)
    }
    tibble(
      table_type = "taxonomic_rank_candidate",
      option = rk,
      n_layers = n_distinct(df$layer_id),
      n_classified_rows = sum(tx_ok),
      pct_classified_rows = mean(tx_ok),
      n_taxa_total = ifelse(nrow(tmp) > 0, n_distinct(tmp$taxon), 0),
      median_taxa_per_layer = ifelse(nrow(layer_taxa) > 0, median(layer_taxa$n_taxa), 0),
      zero_fraction = sparsity
    )
  })

  rank_ok <- rank_eval %>%
    mutate(stable = pct_classified_rows >= 0.35 & median_taxa_per_layer >= 5)

  if ("genus" %in% rank_ok$option) {
    g <- rank_ok %>% filter(option == "genus")
    chosen <- ifelse(nrow(g) == 1 && isTRUE(g$stable), "genus", "family")
  } else {
    chosen <- rank_candidates[1]
  }
  if (!(chosen %in% rank_candidates)) chosen <- rank_candidates[1]
  list(chosen_rank = chosen, rank_eval = rank_eval, rank_ok = rank_ok)
}

build_community <- function(df, chosen_rank, domain_set = c("d__Bacteria", "d__Archaea"), label = "bacteria_archaea") {
  use <- df %>% filter(domain %in% domain_set)
  tax <- as.character(use[[chosen_rank]])
  tax_ok <- !is.na(tax) & tax != "" & !str_detect(tolower(tax), "^unclassified|^unknown|^uncultured")

  lay <- use %>%
    mutate(taxon = ifelse(tax_ok, tax, NA_character_)) %>%
    filter(!is.na(taxon)) %>%
    group_by(layer_id, core_clean, depth_in_core_cm, y_bp, taxon) %>%
    summarise(count = sum(abundance, na.rm = TRUE), .groups = "drop")

  if (nrow(lay) == 0) stop("No classified taxa after filters for domain set: ", label)

  total_layer <- lay %>% group_by(layer_id) %>% summarise(total_abundance = sum(count), .groups = "drop")
  keep_layer <- total_layer %>% filter(total_abundance >= params$min_layer_total_abundance)

  lay <- lay %>% semi_join(keep_layer, by = "layer_id")
  total_layer <- total_layer %>% semi_join(keep_layer, by = "layer_id")

  n_layers <- n_distinct(lay$layer_id)
  prev <- lay %>% count(taxon, name = "prev_layers") %>%
    mutate(prev_frac = prev_layers / n_layers)
  prev_min <- max(params$min_taxon_prevalence_layers, ceiling(params$min_taxon_prevalence_fraction * n_layers))
  keep_taxa <- prev %>% filter(prev_layers >= prev_min) %>% pull(taxon)

  lay2 <- lay %>% filter(taxon %in% keep_taxa)
  wide <- lay2 %>% select(layer_id, taxon, count) %>%
    pivot_wider(names_from = taxon, values_from = count, values_fill = 0)

  mat <- as.matrix(wide[, setdiff(names(wide), "layer_id"), drop = FALSE])
  rownames(mat) <- wide$layer_id

  meta <- lay2 %>%
    distinct(layer_id, core_clean, depth_in_core_cm, y_bp) %>%
    left_join(total_layer, by = "layer_id")

  alpha <- tibble(
    layer_id = rownames(mat),
    richness = rowSums(mat > 0),
    shannon = apply(mat, 1, shannon_from_counts)
  )

  list(
    label = label,
    domain_set = paste(domain_set, collapse = ";"),
    lay = lay2,
    matrix = mat,
    meta = meta,
    alpha = alpha,
    prev = prev,
    prev_min = prev_min,
    n_layers = nrow(meta),
    n_taxa = ncol(mat)
  )
}

# ------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------
dmg <- read_tsv(damage_path, show_col_types = FALSE)
xrf <- read_csv(xrf_path, show_col_types = FALSE)
meta <- read_tsv(meta_path, show_col_types = FALSE)

foram <- if (file.exists(foram_path)) read_tsv(foram_path, show_col_types = FALSE) else tibble()
sst <- if (file.exists(sst_path)) read_tsv(sst_path, show_col_types = FALSE) else tibble()
geob <- if (file.exists(geob_path)) read_tsv(geob_path, show_col_types = FALSE) else tibble()

# ------------------------------------------------------------------
# Damaged-only + explicit abundance usage
# ------------------------------------------------------------------
if (!("abundance" %in% names(dmg))) stop("Required field 'abundance' not found.")

dmg <- dmg %>%
  mutate(
    core_clean = normalize_core(core),
    depth_in_core_cm = suppressWarnings(as.numeric(depth_in_core_cm)),
    y_bp = suppressWarnings(as.numeric(y_bp)),
    abundance = suppressWarnings(as.numeric(abundance))
  ) %>%
  filter(is_dmg == "Damaged") %>%
  filter(is.finite(abundance), abundance > 0) %>%
  filter(!is.na(core_clean), is.finite(depth_in_core_cm), is.finite(y_bp)) %>%
  mutate(layer_id = paste(core_clean, depth_in_core_cm, y_bp, sep = "|"))

if (nrow(dmg) == 0) stop("No rows after damaged-only + abundance > 0 filtering.")

# ------------------------------------------------------------------
# Rank choice (genus preferred, family fallback)
# ------------------------------------------------------------------
rank_res <- pick_taxon_rank(dmg, c("genus", "family"))
chosen_rank <- rank_res$chosen_rank

# ------------------------------------------------------------------
# Community tables: primary Bacteria+Archaea, sensitivity includes viruses
# ------------------------------------------------------------------
comm_primary <- build_community(
  df = dmg,
  chosen_rank = chosen_rank,
  domain_set = c("d__Bacteria", "d__Archaea"),
  label = "bacteria_archaea"
)

comm_sensitivity <- build_community(
  df = dmg,
  chosen_rank = chosen_rank,
  domain_set = c("d__Bacteria", "d__Archaea", "d__Viruses"),
  label = "bacteria_archaea_viruses"
)

# ------------------------------------------------------------------
# Layer metadata + XRF aggregation + productivity axes
# ------------------------------------------------------------------
layer_meta <- comm_primary$meta %>%
  left_join(comm_primary$alpha, by = "layer_id") %>%
  left_join(
    meta %>%
      mutate(core_clean = normalize_core(core)) %>%
      group_by(core_clean, depth_in_core_cm, y_bp) %>%
      summarise(n_metadata_rows = n(), .groups = "drop"),
    by = c("core_clean", "depth_in_core_cm", "y_bp")
  )

xrf <- xrf %>%
  mutate(
    core_clean = normalize_core(core),
    depth_in_core_cm = suppressWarnings(as.numeric(depth_in_core_cm)),
    y_bp = suppressWarnings(as.numeric(y_bp)),
    ba_ti = if ("ba_ti" %in% names(.)) ba_ti else safe_ratio(ba, ti),
    p_ti = if ("p_ti" %in% names(.)) p_ti else safe_ratio(p, ti),
    br_ti = if ("br_ti" %in% names(.)) br_ti else safe_ratio(br, ti),
    si_ti = if ("si_ti" %in% names(.)) si_ti else safe_ratio(si, ti),
    ca_ti = if ("ca_ti" %in% names(.)) ca_ti else safe_ratio(ca, ti),
    rb_sr = if ("rb_sr" %in% names(.)) rb_sr else safe_ratio(rb, sr),
    al_si = if ("al_si" %in% names(.)) al_si else safe_ratio(al, si),
    ti_al = if ("ti_al" %in% names(.)) ti_al else safe_ratio(ti, al),
    mn_fe = if ("mn_fe" %in% names(.)) mn_fe else safe_ratio(mn, fe),
    fe_al = if ("fe_al" %in% names(.)) fe_al else safe_ratio(fe, al)
  )

xrf_keep <- intersect(
  c("ba_ti", "p_ti", "br_ti", "si_ti", "ca_ti", "ca_area", "rb_sr", "al_si", "ti_al", "rb", "zr", "mn_fe", "fe_al", "mn", "fe"),
  names(xrf)
)

xrf_layer <- map_dfr(unique(layer_meta$core_clean), function(cc) {
  lc <- layer_meta %>% filter(core_clean == cc)
  xc <- xrf %>% filter(core_clean == cc, is.finite(depth_in_core_cm))
  if (nrow(lc) == 0) return(tibble())
  agg <- aggregate_dense_to_layers(lc, xc, xrf_keep)
  left_join(lc %>% select(core_clean, depth_in_core_cm, y_bp, layer_id), agg, by = "depth_in_core_cm")
})

layer_meta <- layer_meta %>%
  left_join(xrf_layer %>% select(-core_clean, -y_bp), by = c("layer_id", "depth_in_core_cm"))

biogenic_res <- build_axis(
  df = layer_meta,
  vars = c("ba_ti", "p_ti", "br_ti", "si_ti"),
  axis_name = "biogenic_export_like",
  orientation_priority = c("ba_ti", "p_ti", "br_ti", "si_ti")
)
carbonate_res <- build_axis(
  df = layer_meta,
  vars = c("ca_ti", "ca_area"),
  axis_name = "carbonate",
  orientation_priority = c("ca_ti", "ca_area")
)
terrig_res <- build_axis(
  df = layer_meta,
  vars = c("rb_sr", "al_si", "ti_al", "rb", "zr", "si_ti"),
  axis_name = "terrigenous_mineralogical",
  orientation_priority = c("rb_sr", "al_si", "ti_al", "rb", "zr")
)
redox_res <- build_axis(
  df = layer_meta,
  vars = c("mn_fe", "fe_al", "mn", "fe"),
  axis_name = "redox",
  orientation_priority = c("mn_fe", "fe_al", "mn", "fe")
)

layer_meta <- layer_meta %>%
  mutate(
    axis_biogenic_export = biogenic_res$score,
    axis_carbonate = carbonate_res$score,
    axis_terrigenous = terrig_res$score,
    axis_redox = redox_res$score
  )

# ------------------------------------------------------------------
# Ordination + multivariate testing (vegan)
# ------------------------------------------------------------------
comm_mat <- comm_primary$matrix
common_layers <- intersect(rownames(comm_mat), layer_meta$layer_id)
comm_mat <- comm_mat[common_layers, , drop = FALSE]
layer_meta <- layer_meta %>% filter(layer_id %in% common_layers)

clr_obj <- clr_transform(comm_mat)
clr_mat <- clr_obj$clr

ok_ord <- rowSums(is.finite(clr_mat)) == ncol(clr_mat)
if (sum(ok_ord) < 20) stop("Too few complete CLR rows for robust ordination/multivariate analysis.")

clr_ok <- clr_mat[ok_ord, , drop = FALSE]
meta_ok <- layer_meta %>% filter(layer_id %in% rownames(clr_ok))

pca <- prcomp(clr_ok, center = TRUE, scale. = FALSE)
ord_scores <- as_tibble(pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE], rownames = "layer_id")
names(ord_scores) <- sub("PC", "ord_pc", names(ord_scores))
layer_ord <- meta_ok %>% left_join(ord_scores, by = "layer_id")

# Hellinger sensitivity
hell <- vegan::decostand(comm_mat, method = "hellinger")
hell_ok <- hell[rownames(clr_ok), , drop = FALSE]
pc_hell <- cmdscale(vegdist(hell_ok, method = "bray"), k = 2, eig = TRUE)
hell_scores <- tibble(layer_id = rownames(hell_ok), hell_axis1 = pc_hell$points[, 1], hell_axis2 = pc_hell$points[, 2])
layer_ord <- layer_ord %>% left_join(hell_scores, by = "layer_id")

trend_col <- choose_trend(layer_ord)
if (is.na(trend_col)) {
  layer_ord <- layer_ord %>% mutate(trend_fallback = row_number())
  trend_col <- "trend_fallback"
}

adonis_dat <- layer_ord %>%
  select(layer_id, core_clean, all_of(trend_col), axis_biogenic_export, axis_carbonate) %>%
  filter(is.finite(axis_biogenic_export), is.finite(axis_carbonate), is.finite(.data[[trend_col]]))

dist_clr <- dist(clr_ok[adonis_dat$layer_id, , drop = FALSE], method = "euclidean")
form_adonis <- as.formula(paste0("dist_clr ~ core_clean + splines::ns(", trend_col, ", df=4) + axis_biogenic_export + axis_carbonate"))
ad <- vegan::adonis2(form_adonis, data = adonis_dat, permutations = params$permutations, by = "margin")
ad_df <- as.data.frame(ad)
ad_df$term <- rownames(ad_df)

dispersion <- list()
dist_clr_mat <- as.matrix(dist_clr)
for (v in c("axis_biogenic_export", "axis_carbonate")) {
  q <- quantile(adonis_dat[[v]], probs = c(0.25, 0.75), na.rm = TRUE)
  g <- ifelse(adonis_dat[[v]] <= q[1], "low", ifelse(adonis_dat[[v]] >= q[2], "high", "mid"))
  keep <- g %in% c("low", "high")
  if (sum(keep) < 8 || length(unique(g[keep])) < 2) {
    bd <- NULL
  } else {
    sub_dist <- as.dist(dist_clr_mat[keep, keep, drop = FALSE])
    bd <- tryCatch(vegan::betadisper(sub_dist, group = factor(g[keep], levels = c("low", "high"))), error = function(e) NULL)
  }
  if (is.null(bd)) {
    dispersion[[v]] <- tibble(grouping = v, n = sum(keep), f_value = NA_real_, p_value = NA_real_, note = "betadisper not estimable")
  } else {
    bt <- anova(bd)
    dispersion[[v]] <- tibble(grouping = v, n = sum(keep), f_value = bt$`F value`[1], p_value = bt$`Pr(>F)`[1], note = "low/high quartile groups")
  }
}
dispersion_df <- bind_rows(dispersion)

# envfit on ordination
ord_for_envfit <- pca$x[, 1:2, drop = FALSE]
envfit_vars <- layer_ord %>%
  filter(layer_id %in% rownames(ord_for_envfit)) %>%
  select(layer_id, axis_biogenic_export, axis_carbonate, axis_terrigenous, axis_redox, y_bp, depth_in_core_cm)
envfit_vars <- envfit_vars %>%
  filter(is.finite(axis_biogenic_export), is.finite(axis_carbonate))
row_keep <- envfit_vars$layer_id
ef <- vegan::envfit(ord_for_envfit[row_keep, , drop = FALSE], envfit_vars %>% select(-layer_id), permutations = params$permutations)
ef_tbl <- tibble(
  variable = names(ef$vectors$r),
  r2 = as.numeric(ef$vectors$r),
  p_value = as.numeric(ef$vectors$pvals)
) %>% mutate(q_value_bh = p.adjust(p_value, method = "BH"))

# ------------------------------------------------------------------
# Scenario definition from axes
# ------------------------------------------------------------------
scenario_tbl <- layer_ord %>%
  mutate(
    export_state = case_when(
      axis_biogenic_export <= quantile(axis_biogenic_export, 0.25, na.rm = TRUE) ~ "export_low",
      axis_biogenic_export >= quantile(axis_biogenic_export, 0.75, na.rm = TRUE) ~ "export_high",
      TRUE ~ "export_mid"
    ),
    carbonate_state = case_when(
      axis_carbonate <= quantile(axis_carbonate, 0.25, na.rm = TRUE) ~ "carbonate_low",
      axis_carbonate >= quantile(axis_carbonate, 0.75, na.rm = TRUE) ~ "carbonate_high",
      TRUE ~ "carbonate_mid"
    ),
    scenario = case_when(
      export_state == "export_high" & carbonate_state == "carbonate_low" ~ "high_export_low_carb",
      export_state == "export_low" & carbonate_state == "carbonate_high" ~ "low_export_high_carb",
      export_state == "export_high" & carbonate_state == "carbonate_high" ~ "high_export_high_carb",
      export_state == "export_low" & carbonate_state == "carbonate_low" ~ "low_export_low_carb",
      TRUE ~ "intermediate"
    )
  )

scenario_counts <- scenario_tbl %>%
  count(scenario, core_clean, name = "n_layers") %>%
  arrange(desc(n_layers))

scenario_global_counts <- scenario_tbl %>% count(scenario, name = "n_layers_total")
scenario_ok <- scenario_global_counts %>%
  filter(n_layers_total >= params$scenario_min_n, scenario != "intermediate") %>%
  pull(scenario)

# ------------------------------------------------------------------
# Assemblage / indicator analysis
# ------------------------------------------------------------------
clr_tax <- as_tibble(clr_ok, rownames = "layer_id") %>%
  pivot_longer(cols = -layer_id, names_to = "taxon", values_to = "clr_value")

tax_prev <- comm_primary$prev %>%
  rename(taxon = taxon, prevalence_layers = prev_layers, prevalence_frac = prev_frac)

tax_mod_df <- scenario_tbl %>%
  select(layer_id, core_clean, y_bp, depth_in_core_cm, scenario, export_state, carbonate_state,
         axis_biogenic_export, axis_carbonate) %>%
  left_join(clr_tax, by = "layer_id") %>%
  left_join(tax_prev, by = "taxon") %>%
  filter(is.finite(clr_value), prevalence_layers >= comm_primary$prev_min)

trend_use <- choose_trend(tax_mod_df)
if (is.na(trend_use)) {
  tax_mod_df <- tax_mod_df %>% mutate(trend_fallback = row_number())
  trend_use <- "trend_fallback"
}

fit_taxon_contrast <- function(df, contrast_name, contrast_col, positive_label) {
  taxa <- sort(unique(df$taxon))
  res <- map_dfr(taxa, function(tx) {
    d <- df %>% filter(taxon == tx, .data[[contrast_col]] %in% c("other", positive_label))
    if (nrow(d) < 25 || dplyr::n_distinct(d$core_clean) < 2) {
      return(tibble(
        contrast = contrast_name,
        taxon = tx,
        n = nrow(d),
        beta = NA_real_, p_value = NA_real_,
        note = "insufficient_data"
      ))
    }
    d <- d %>% mutate(contrast_flag = ifelse(.data[[contrast_col]] == positive_label, 1, 0))
    fml <- as.formula(paste0("clr_value ~ contrast_flag + core_clean + splines::ns(", trend_use, ", df=4)"))
    fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
    if (is.null(fit) || !("contrast_flag" %in% rownames(summary(fit)$coefficients))) {
      return(tibble(
        contrast = contrast_name,
        taxon = tx,
        n = nrow(d),
        beta = NA_real_, p_value = NA_real_,
        note = "model_not_estimable"
      ))
    }
    co <- summary(fit)$coefficients["contrast_flag", ]
    tibble(
      contrast = contrast_name,
      taxon = tx,
      n = nrow(d),
      beta = co["Estimate"],
      p_value = co["Pr(>|t|)"],
      note = "core+trend_adjusted"
    )
  })
  res %>% mutate(q_value_bh = p.adjust(p_value, method = "BH"))
}

tax_mod_df <- tax_mod_df %>%
  mutate(
    export_contrast = ifelse(export_state == "export_high", "export_high", "other"),
    carbonate_contrast = ifelse(carbonate_state == "carbonate_high", "carbonate_high", "other"),
    he_lc_contrast = ifelse(scenario == "high_export_low_carb", "he_lc", "other")
  )

tax_export <- fit_taxon_contrast(tax_mod_df, "export_high_vs_other", "export_contrast", "export_high")
tax_carb <- fit_taxon_contrast(tax_mod_df, "carbonate_high_vs_other", "carbonate_contrast", "carbonate_high")
tax_he_lc <- fit_taxon_contrast(tax_mod_df, "high_export_low_carb_vs_other", "he_lc_contrast", "he_lc")

tax_assoc_all <- bind_rows(tax_export, tax_carb, tax_he_lc) %>%
  left_join(tax_prev, by = "taxon") %>%
  arrange(contrast, q_value_bh, desc(abs(beta)))

top_tax <- tax_assoc_all %>%
  group_by(contrast) %>%
  filter(!is.na(q_value_bh)) %>%
  slice_head(n = 25) %>%
  ungroup()

# Core-specific direction check
core_specific <- tax_mod_df %>%
  filter(core_clean %in% names(which(table(core_clean) >= 20))) %>%
  group_by(core_clean, taxon) %>%
  summarise(
    n = n(),
    rho_export = safe_spearman(clr_value, axis_biogenic_export, min_n = 10)$estimate,
    p_export = safe_spearman(clr_value, axis_biogenic_export, min_n = 10)$p,
    rho_carb = safe_spearman(clr_value, axis_carbonate, min_n = 10)$estimate,
    p_carb = safe_spearman(clr_value, axis_carbonate, min_n = 10)$p,
    .groups = "drop"
  ) %>%
  group_by(core_clean) %>%
  mutate(
    q_export = p.adjust(p_export, method = "BH"),
    q_carb = p.adjust(p_carb, method = "BH")
  ) %>%
  ungroup()

core_dir_summary <- core_specific %>%
  filter(taxon %in% unique(top_tax$taxon)) %>%
  group_by(taxon) %>%
  summarise(
    cores_tested = n_distinct(core_clean),
    n_pos_export = sum(rho_export > 0, na.rm = TRUE),
    n_neg_export = sum(rho_export < 0, na.rm = TRUE),
    n_sig_export = sum(q_export < 0.10, na.rm = TRUE),
    n_pos_carb = sum(rho_carb > 0, na.rm = TRUE),
    n_neg_carb = sum(rho_carb < 0, na.rm = TRUE),
    n_sig_carb = sum(q_carb < 0.10, na.rm = TRUE),
    .groups = "drop"
  )

scenario_assemblages <- top_tax %>%
  mutate(direction = ifelse(beta > 0, "enriched", "depleted")) %>%
  left_join(core_dir_summary, by = "taxon") %>%
  select(contrast, taxon, prevalence_layers, prevalence_frac, n, beta, p_value, q_value_bh,
         direction, cores_tested, n_pos_export, n_neg_export, n_sig_export, n_pos_carb, n_neg_carb, n_sig_carb, note)

# ------------------------------------------------------------------
# Alpha associations
# ------------------------------------------------------------------
alpha_assoc <- crossing(metric = c("richness", "shannon"), axis = c("axis_biogenic_export", "axis_carbonate")) %>%
  mutate(
    result = map2(metric, axis, ~safe_spearman(scenario_tbl[[.x]], scenario_tbl[[.y]], min_n = 12)),
    n = map_int(result, "n"),
    estimate = map_dbl(result, "estimate"),
    p_value = map_dbl(result, "p")
  ) %>%
  mutate(
    model = "pooled_spearman",
    q_value_bh = p.adjust(p_value, method = "BH")
  ) %>%
  select(metric, axis, n, estimate, p_value, q_value_bh, model)

alpha_adj <- map_dfr(c("richness", "shannon"), function(resp) {
  map_dfr(c("axis_biogenic_export", "axis_carbonate"), function(ax) {
    d <- scenario_tbl %>%
      select(core_clean, all_of(resp), all_of(ax), y_bp, depth_in_core_cm) %>%
      filter(is.finite(.data[[resp]]), is.finite(.data[[ax]]), is.finite(.data[[trend_col]]))
    if (nrow(d) < 20 || dplyr::n_distinct(d$core_clean) < 2) {
      return(tibble(metric = resp, axis = ax, n = nrow(d), beta = NA_real_, p_value = NA_real_, note = "insufficient_data"))
    }
    f <- as.formula(paste0(resp, " ~ ", ax, " + core_clean + splines::ns(", trend_col, ", df=4)"))
    fit <- tryCatch(lm(f, data = d), error = function(e) NULL)
    if (is.null(fit) || !(ax %in% rownames(summary(fit)$coefficients))) {
      return(tibble(metric = resp, axis = ax, n = nrow(d), beta = NA_real_, p_value = NA_real_, note = "model_not_estimable"))
    }
    co <- summary(fit)$coefficients[ax, ]
    tibble(metric = resp, axis = ax, n = nrow(d), beta = co["Estimate"], p_value = co["Pr(>|t|)"], note = "core+trend_adjusted")
  })
}) %>% mutate(q_value_bh = p.adjust(p_value, method = "BH"))

alpha_out <- alpha_assoc %>%
  left_join(alpha_adj %>% select(metric, axis, beta, p_value_adj = p_value, q_value_adj_bh = q_value_bh, note_adj = note),
            by = c("metric", "axis"))

# ------------------------------------------------------------------
# Hypothesis generation + conservative tests with contextual proxies
# ------------------------------------------------------------------
build_context_tbl <- function(base_tbl) {
  out <- base_tbl %>% select(layer_id, core_clean, depth_in_core_cm, y_bp, scenario,
                             axis_biogenic_export, axis_carbonate, axis_terrigenous, axis_redox)

  if (nrow(foram) > 0) {
    f2 <- foram %>% mutate(core_clean = normalize_core(core), depth_in_core_cm = suppressWarnings(as.numeric(depth_in_core_cm)), y_bp = suppressWarnings(as.numeric(y_bp))) %>%
      group_by(core_clean, depth_in_core_cm, y_bp) %>%
      summarise(
        p_to_b_ratio = if ("p_to_b_ratio" %in% names(.)) median(p_to_b_ratio, na.rm = TRUE) else NA_real_,
        benthic_foram = if ("benthic_foraminifera_concentration_n_per_gram_gt_150um" %in% names(.)) median(benthic_foraminifera_concentration_n_per_gram_gt_150um, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
    out <- out %>% left_join(f2, by = c("core_clean", "depth_in_core_cm", "y_bp"))
  }

  if (nrow(sst) > 0) {
    s2 <- sst %>% mutate(core_clean = normalize_core(core), depth_in_core_cm = suppressWarnings(as.numeric(depth_in_core_cm))) %>%
      group_by(core_clean, depth_in_core_cm) %>%
      summarise(across(starts_with("sst_"), ~median(as.numeric(.x), na.rm = TRUE), .names = "{.col}"), .groups = "drop")
    out <- out %>% left_join(s2, by = c("core_clean", "depth_in_core_cm"))
  }

  if (nrow(geob) > 0) {
    g2 <- geob %>% mutate(core_clean = normalize_core(core), depth_in_core_cm = suppressWarnings(as.numeric(depth_in_core_cm)), y_bp = suppressWarnings(as.numeric(y_bp))) %>%
      group_by(core_clean, depth_in_core_cm, y_bp) %>%
      summarise(
        toc_percent_weight = if ("toc_percent_weight" %in% names(.)) median(toc_percent_weight, na.rm = TRUE) else NA_real_,
        ratio_greenland_to_iceland = if ("ratio_greenland_to_iceland" %in% names(.)) median(ratio_greenland_to_iceland, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
    out <- out %>% left_join(g2, by = c("core_clean", "depth_in_core_cm", "y_bp"))
  }

  out
}

context_tbl <- build_context_tbl(scenario_tbl)

# Functional bins: cautious broad labels based on names only
functional_bins <- tibble(
  pattern = c("Nitrosopum", "Nitrosopelag", "Nitrospin", "Thioglob", "Pelagibacter", "Polaribacter", "Woeseia", "Desulf", "Methan", "Akkerman"),
  func_bin = c("nitrifier_like_oxic", "nitrifier_like_oxic", "nitrite_oxidizer_like", "sulfur_oxidizer_like", "oligotroph_like", "copiotroph_like", "sulfur_nitrogen_cycler_like", "sulfate_reducer_like", "methanogen_like", "anaerobe_like")
)

tax_bin <- tax_prev %>%
  mutate(func_bin = map_chr(taxon, function(tx) {
    hit <- functional_bins$func_bin[str_detect(tx, functional_bins$pattern)]
    if (length(hit) == 0) return("unassigned_broad")
    hit[1]
  }))

bin_scores <- clr_tax %>%
  left_join(tax_bin %>% select(taxon, func_bin), by = "taxon") %>%
  group_by(layer_id, func_bin) %>%
  summarise(bin_clr_mean = mean(clr_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = func_bin, values_from = bin_clr_mean)

hyp_dat <- context_tbl %>% left_join(bin_scores, by = "layer_id")

test_hyp <- function(hypothesis_id, response, predictor, adjust_core_trend = TRUE, note = "") {
  if (!(response %in% names(hyp_dat)) || !(predictor %in% names(hyp_dat))) {
    return(tibble(hypothesis_id = hypothesis_id, response = response, predictor = predictor,
                  n = 0, estimate = NA_real_, p_value = NA_real_, q_value_bh = NA_real_,
                  model = "not_run", support = "not_testable", note = paste(note, "missing_variable")))
  }
  d <- hyp_dat %>%
    select(core_clean, y_bp, depth_in_core_cm, all_of(response), all_of(predictor)) %>%
    filter(is.finite(.data[[response]]), is.finite(.data[[predictor]]), is.finite(.data[[trend_col]]))
  if (nrow(d) < 20) {
    return(tibble(hypothesis_id = hypothesis_id, response = response, predictor = predictor,
                  n = nrow(d), estimate = NA_real_, p_value = NA_real_, q_value_bh = NA_real_,
                  model = "not_run_low_n", support = "not_testable", note = paste(note, "low_n")))
  }
  if (adjust_core_trend && dplyr::n_distinct(d$core_clean) >= 2) {
    f <- as.formula(paste0(response, " ~ ", predictor, " + core_clean + splines::ns(", trend_col, ", df=4)"))
    fit <- tryCatch(lm(f, data = d), error = function(e) NULL)
    if (is.null(fit) || !(predictor %in% rownames(summary(fit)$coefficients))) {
      return(tibble(hypothesis_id = hypothesis_id, response = response, predictor = predictor,
                    n = nrow(d), estimate = NA_real_, p_value = NA_real_, q_value_bh = NA_real_,
                    model = "lm_core_trend_failed", support = "inconclusive", note = paste(note, "model_failure")))
    }
    co <- summary(fit)$coefficients[predictor, ]
    return(tibble(hypothesis_id = hypothesis_id, response = response, predictor = predictor,
                  n = nrow(d), estimate = co["Estimate"], p_value = co["Pr(>|t|)"], q_value_bh = NA_real_,
                  model = "lm_core_plus_trend", support = "to_classify", note = note))
  }
  sp <- safe_spearman(d[[response]], d[[predictor]], min_n = 15)
  tibble(hypothesis_id = hypothesis_id, response = response, predictor = predictor,
         n = sp$n, estimate = sp$estimate, p_value = sp$p, q_value_bh = NA_real_,
         model = "spearman", support = "to_classify", note = paste(note, "no_core_adjustment"))
}

hyp_tests <- bind_rows(
  test_hyp("H1_export_redox", "axis_redox", "axis_biogenic_export", note = "High export expected to co-vary with redox shifts"),
  test_hyp("H2_export_toc", "toc_percent_weight", "axis_biogenic_export", note = "High export expected to align with higher TOC where available"),
  test_hyp("H3_carb_terrigenous_inverse", "axis_terrigenous", "axis_carbonate", note = "Carbonate high expected to oppose terrigenous signal"),
  test_hyp("H4_nitrifier_carb", "nitrifier_like_oxic", "axis_carbonate", note = "Carbonate-rich states may support oxic/nitrifier-like signals"),
  test_hyp("H5_nitrifier_export", "nitrifier_like_oxic", "axis_biogenic_export", note = "Export-rich states may reduce oxic/nitrifier-like signature"),
  test_hyp("H6_copiotroph_export", "copiotroph_like", "axis_biogenic_export", note = "Export-rich states may enrich copiotroph-like taxa"),
  test_hyp("H7_sulfur_export", "sulfur_oxidizer_like", "axis_biogenic_export", note = "Export-rich states may co-vary with sulfur-cycler-like taxa")
)

hyp_tests <- hyp_tests %>%
  mutate(q_value_bh = p.adjust(p_value, method = "BH")) %>%
  mutate(
    support = case_when(
      !is.finite(p_value) ~ "not_testable_or_inconclusive",
      q_value_bh <= 0.10 & estimate > 0 ~ "supported_positive_direction",
      q_value_bh <= 0.10 & estimate < 0 ~ "supported_negative_direction",
      q_value_bh > 0.10 ~ "not_supported_in_this_dataset",
      TRUE ~ "inconclusive"
    )
  )

# ------------------------------------------------------------------
# Optional sensitivity: viruses included
# ------------------------------------------------------------------
sens_note <- {
  m <- comm_sensitivity$matrix
  common <- intersect(rownames(m), rownames(comm_primary$matrix))
  m <- m[common, , drop = FALSE]
  p <- comm_primary$matrix[common, , drop = FALSE]
  clr_m <- clr_transform(m)$clr
  clr_p <- clr_transform(p)$clr
  pca_m <- prcomp(clr_m, center = TRUE, scale. = FALSE)
  pca_p <- prcomp(clr_p, center = TRUE, scale. = FALSE)
  sc <- suppressWarnings(cor(pca_m$x[, 1], pca_p$x[, 1], method = "spearman", use = "pairwise.complete.obs"))
  paste0("Virus-included sensitivity: Spearman correlation between primary and sensitivity CLR-PC1 = ", signif(sc, 4), ".")
}

# ------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------
rank_choice_out <- bind_rows(
  tibble(
    table_type = "count_field_candidate",
    option = "abundance",
    n_non_missing = sum(is.finite(dmg$abundance)),
    pct_non_missing = mean(is.finite(dmg$abundance)),
    pct_non_negative = mean(dmg$abundance >= 0),
    pct_integerish = mean(abs(dmg$abundance - round(dmg$abundance)) < 1e-8),
    n_positive = sum(dmg$abundance > 0),
    median_positive = median(dmg$abundance[dmg$abundance > 0]),
    chosen = TRUE,
    reason = "Forced by corrected requirement: use abundance (not n_reads)."
  ),
  rank_res$rank_ok %>%
    mutate(
      chosen = option == chosen_rank,
      reason = ifelse(chosen, "Chosen for interpretability + stability criteria", "Not selected")
    )
)

layer_metadata_out <- scenario_tbl %>%
  select(
    core_clean, depth_in_core_cm, y_bp, layer_id,
    total_abundance,
    n_metadata_rows,
    n_dense_rows,
    ba_ti, p_ti, br_ti, si_ti, ca_ti, ca_area,
    rb_sr, al_si, ti_al, mn_fe, fe_al,
    axis_biogenic_export, axis_carbonate, axis_terrigenous, axis_redox,
    richness, shannon,
    export_state, carbonate_state, scenario
  ) %>%
  arrange(core_clean, depth_in_core_cm, y_bp)

axis_load <- bind_rows(
  biogenic_res$loadings,
  carbonate_res$loadings,
  terrig_res$loadings,
  redox_res$loadings
)

axis_methods <- tibble(
  axis = c("biogenic_export_like", "carbonate", "terrigenous_mineralogical", "redox"),
  method = c(biogenic_res$method, carbonate_res$method, terrig_res$method, redox_res$method),
  anchor_variable = c(biogenic_res$anchor, carbonate_res$anchor, terrig_res$anchor, redox_res$anchor),
  orientation_rho_before = c(biogenic_res$rho_before, carbonate_res$rho_before, terrig_res$rho_before, redox_res$rho_before)
)

prod_axes_out <- scenario_tbl %>%
  select(core_clean, depth_in_core_cm, y_bp, layer_id,
         axis_biogenic_export, axis_carbonate, axis_terrigenous, axis_redox,
         ba_ti, p_ti, br_ti, si_ti, ca_ti, ca_area, rb_sr, al_si, ti_al, mn_fe, fe_al)

axis_methods_wide <- axis_methods %>%
  mutate(tmp = 1) %>%
  pivot_wider(names_from = axis, values_from = c(method, anchor_variable, orientation_rho_before)) %>%
  select(-tmp)

prod_axes_out <- bind_cols(
  prod_axes_out,
  axis_methods_wide[rep(1, nrow(prod_axes_out)), , drop = FALSE]
)

multivar_summary <- bind_rows(
  ad_df %>%
    transmute(
      section = "adonis2",
      term = term,
      df = Df,
      statistic = F,
      r2 = R2,
      p_value = `Pr(>F)`,
      note = paste0("Model: core + spline(", trend_col, ") + export axis + carbonate axis")
    ),
  dispersion_df %>%
    transmute(
      section = "betadisper",
      term = grouping,
      df = NA_real_,
      statistic = f_value,
      r2 = NA_real_,
      p_value = p_value,
      note = note
    ),
  ef_tbl %>%
    transmute(
      section = "envfit",
      term = variable,
      df = NA_real_,
      statistic = r2,
      r2 = r2,
      p_value = p_value,
      note = paste0("envfit q=" , signif(q_value_bh, 3))
    )
)

scenario_definition <- scenario_tbl %>%
  select(layer_id, core_clean, depth_in_core_cm, y_bp,
         axis_biogenic_export, axis_carbonate,
         export_state, carbonate_state, scenario)

scenario_q <- scenario_tbl %>%
  summarise(
    export_q25 = quantile(axis_biogenic_export, 0.25, na.rm = TRUE),
    export_q75 = quantile(axis_biogenic_export, 0.75, na.rm = TRUE),
    carbonate_q25 = quantile(axis_carbonate, 0.25, na.rm = TRUE),
    carbonate_q75 = quantile(axis_carbonate, 0.75, na.rm = TRUE)
  )

scenario_definition <- bind_cols(
  scenario_definition,
  scenario_q[rep(1, nrow(scenario_definition)), , drop = FALSE]
)

top_tax_out <- tax_assoc_all %>%
  mutate(direction = ifelse(beta > 0, "positive", "negative")) %>%
  arrange(contrast, q_value_bh, desc(abs(beta)))

# Methods + interpretation notes
methods_notes <- c(
  "xrf_microbiome methods notes",
  "",
  paste0("Primary count field: ", params$count_field_forced, " (explicitly forced by user correction)."),
  paste0("Damaged-only filter: ", params$damaged_filter, "."),
  paste0("Primary domains: ", params$primary_domains, ". Viruses run as sensitivity only."),
  paste0("Taxonomic rank chosen by stability rule: ", chosen_rank, "."),
  paste0("Layer minimum total abundance: ", params$min_layer_total_abundance, "."),
  paste0("Rare-taxon filter: prevalence >= ", comm_primary$prev_min, " layers."),
  paste0("CLR pseudocount: ", signif(clr_obj$pseudocount, 6), " (", params$clr_pseudocount_rule, ")."),
  "Export/biogenic axis built from Ba/Ti, P/Ti, Br/Ti, Si/Ti (PCA PC1 if possible, else z-mean).",
  "Carbonate axis built from Ca/Ti (+ Ca_area if available).",
  "Optional contextual axes built: terrigenous/mineralogical and redox.",
  paste0("Multivariate workflow (vegan): adonis2 + betadisper + envfit; permutations=", params$permutations, "."),
  "Scenario definition: quartile states (Q1=low, Q4=high) for export and carbonate axes.",
  "Taxon scenario contrasts: CLR taxon ~ scenario contrast + core + spline(trend), BH-FDR controlled.",
  "Hypothesis tests: conservative core+trend-adjusted models linking axes/functional bins/contextual proxies."
)

interpretation_notes <- c(
  "xrf_microbiome interpretation notes",
  "",
  "Guardrails:",
  "- Damaged-only filtering shifts interpretation toward preserved/authentic-leaning fractions.",
  "- Associations are compositional and observational, not causal demonstrations.",
  "- XRF axis interpretation is indirect and partly confounded by lithology and core structure.",
  "- Time/core structure remains a major confounder despite adjustment.",
  "- Some taxa may reflect persistence, preservation, or contamination risk.",
  "",
  "Literature-aligned ecological framing used here (without overclaiming):",
  "- Nitrifier-like/oxic tendencies: taxa names matching Nitrosopumilus/Nitrosopelagicus/Nitrospin-like groups.",
  "- Heterotroph/copiotroph-like tendencies: Polaribacter and broad Flavobacteriaceae-like signals.",
  "- Sulfur-cycle-like tendencies: Thioglobaceae/Pseudothioglobus/Woeseia-like patterns.",
  "- Deep-biosphere/reducing labels remain speculative at broad genus/family taxonomic resolution.",
  "",
  "Caution on overinterpretation:",
  "- High-level taxonomic matches do not guarantee specific metabolic function in ancient sediments.",
  "- Functional bins are coarse and should be treated as hypothesis scaffolds, not direct proof.",
  sens_note
)

# Write tables
write_tsv(layer_metadata_out, file.path(out_dir, "xrf_microbiome_layer_metadata.tsv"))
write_tsv(rank_choice_out, file.path(out_dir, "xrf_microbiome_rank_choice_and_counts.tsv"))
write_tsv(prod_axes_out, file.path(out_dir, "xrf_microbiome_productivity_axes.tsv"))
write_tsv(alpha_out, file.path(out_dir, "xrf_microbiome_alpha_associations.tsv"))
write_tsv(multivar_summary, file.path(out_dir, "xrf_microbiome_multivariate_test_summary.tsv"))
write_tsv(scenario_definition, file.path(out_dir, "xrf_microbiome_scenario_definition.tsv"))
write_tsv(scenario_assemblages, file.path(out_dir, "xrf_microbiome_scenario_assemblages.tsv"))
write_tsv(top_tax_out, file.path(out_dir, "xrf_microbiome_top_taxon_associations.tsv"))
write_tsv(hyp_tests, file.path(out_dir, "xrf_microbiome_hypothesis_tests.tsv"))

writeLines(interpretation_notes, con = file.path(out_dir, "xrf_microbiome_interpretation_notes.txt"))
writeLines(methods_notes, con = file.path(out_dir, "xrf_microbiome_methods_notes.txt"))

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
plot_ord <- layer_ord %>% filter(is.finite(ord_pc1), is.finite(ord_pc2))

p1 <- ggplot(plot_ord, aes(ord_pc1, ord_pc2, color = axis_biogenic_export, shape = core_clean)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_viridis_c(option = "C", na.value = "grey70") +
  theme_bw() +
  labs(
    title = "Damaged-only microbiome ordination (CLR-PCA)",
    subtitle = "Colored by export/biogenic axis",
    x = "Ordination PC1", y = "Ordination PC2", color = "Export axis", shape = "Core"
  )
ggsave(file.path(plot_dir, "xrf_microbiome_ordination_biogenic_axis.png"), p1, width = 9, height = 6, dpi = 300)

p2 <- ggplot(plot_ord, aes(ord_pc1, ord_pc2, color = axis_carbonate, shape = core_clean)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_viridis_c(option = "B", na.value = "grey70") +
  theme_bw() +
  labs(
    title = "Damaged-only microbiome ordination (CLR-PCA)",
    subtitle = "Colored by carbonate axis",
    x = "Ordination PC1", y = "Ordination PC2", color = "Carbonate axis", shape = "Core"
  )
ggsave(file.path(plot_dir, "xrf_microbiome_ordination_carbonate_axis.png"), p2, width = 9, height = 6, dpi = 300)

scenario_plot_df <- scenario_tbl %>%
  filter(scenario != "intermediate") %>%
  count(core_clean, scenario, name = "n")

p3 <- ggplot(scenario_plot_df, aes(x = scenario, y = n, fill = core_clean)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "Scenario layer counts by core", x = "Scenario", y = "Number of layers", fill = "Core")
ggsave(file.path(plot_dir, "xrf_microbiome_scenario_assemblage_summary.png"), p3, width = 9, height = 5.8, dpi = 300)

heat_df <- top_tax %>%
  mutate(
    contrast = recode(contrast,
                      export_high_vs_other = "Export high vs other",
                      carbonate_high_vs_other = "Carbonate high vs other",
                      high_export_low_carb_vs_other = "High export + low carb vs other"),
    taxon = fct_reorder(taxon, beta)
  )

p4 <- ggplot(heat_df, aes(x = contrast, y = taxon, fill = beta)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_bw() +
  labs(title = "Top taxa associated with XRF productivity scenarios",
       subtitle = "CLR association coefficients (core + trend adjusted)",
       x = "Scenario contrast", y = paste0("Taxon (", chosen_rank, ")"), fill = "Beta")
ggsave(file.path(plot_dir, "xrf_microbiome_top_taxa_indicator_heatmap.png"), p4, width = 10, height = 8.5, dpi = 300)

ctx_vars <- intersect(c("axis_redox", "axis_terrigenous", "toc_percent_weight", "p_to_b_ratio"), names(context_tbl))
if (length(ctx_vars) > 0) {
  ctx_long <- context_tbl %>%
    select(scenario, all_of(ctx_vars)) %>%
    pivot_longer(cols = -scenario, names_to = "proxy", values_to = "value") %>%
    filter(is.finite(value), scenario != "intermediate")
  p5 <- ggplot(ctx_long, aes(x = scenario, y = value, fill = scenario)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.9) +
    facet_wrap(~proxy, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none") +
    labs(title = "Contextual proxy distributions across scenarios", x = "Scenario", y = "Proxy value")
  ggsave(file.path(plot_dir, "xrf_microbiome_scenario_contextual_proxies.png"), p5, width = 12, height = 7, dpi = 300)
}

message("Done. Reworked XRF microbiome analysis complete.")
message("Outputs: ", out_dir)
message("Plots: ", plot_dir)
