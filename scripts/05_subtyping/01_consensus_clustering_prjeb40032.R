#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(sprintf("Missing value for %s", flag))
  args[[idx + 1]]
}

flag_present <- function(flag) {
  any(args == flag)
}

expr_path <- get_arg("--expr-matrix", "data/processed/prjeb40032_expression_matrix_all42.tsv")
meta_path <- get_arg("--metadata", "data/processed/prjeb40032_sample_metadata_all42.tsv")
gene_list_path <- get_arg("--genes", "results/tables/T3.3_human_intersection_genes.csv")

case_label <- get_arg("--case-label", "case")
max_k <- as.integer(get_arg("--max-k", "6"))
reps <- as.integer(get_arg("--reps", "500"))
p_item <- as.numeric(get_arg("--p-item", "0.8"))
p_feature <- as.numeric(get_arg("--p-feature", "1.0"))
seed <- as.integer(get_arg("--seed", "42"))

out_dir_tables <- get_arg("--out-tables", "results/tables")
out_dir_figs <- get_arg("--out-figs", "figures/raw_plots")
prefix <- get_arg("--prefix", "prjeb40032")

out_assignments <- get_arg("--assignments-name", sprintf("T5.1_cluster_assignments_%s.csv", prefix))
out_summary <- get_arg("--summary-name", sprintf("T5.1_consensus_clustering_summary_%s.csv", prefix))
out_fig_heatmap <- get_arg("--heatmap-name", sprintf("Fig6A.png", prefix))
out_fig_cdf <- get_arg("--cdf-name", sprintf("Fig6B.png", prefix))
out_fig_subtype_heatmap <- get_arg("--subtype-heatmap-name", "Fig6C.png")
out_fig_subtype_markers <- get_arg("--subtype-markers-name", "Fig6D.png")

dry_run <- flag_present("--dry-run")

if (!file.exists(expr_path)) stop(sprintf("Expression matrix not found: %s", expr_path))
if (!file.exists(meta_path)) stop(sprintf("Metadata not found: %s", meta_path))
if (!file.exists(gene_list_path)) stop(sprintf("Gene list not found: %s", gene_list_path))
if (!is.finite(max_k) || max_k < 2) stop("--max-k must be >= 2")
if (!is.finite(reps) || reps < 50) stop("--reps must be >= 50")
if (!is.finite(p_item) || p_item <= 0 || p_item > 1) stop("--p-item must be in (0, 1]")
if (!is.finite(p_feature) || p_feature <= 0 || p_feature > 1) stop("--p-feature must be in (0, 1]")

dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figs, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ConsensusClusterPlus)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(seed)

read_gene_symbols <- function(path) {
  df <- readr::read_csv(path, col_types = cols(), progress = FALSE)
  if (!("gene_symbol" %in% colnames(df))) stop("Gene list must contain column: gene_symbol")
  genes <- df$gene_symbol |> as.character() |> trimws()
  genes <- genes[genes != ""]
  unique(genes)
}

read_expression_matrix <- function(path) {
  df <- readr::read_tsv(path, col_types = cols(), progress = FALSE)
  if (!("gene_symbol" %in% colnames(df))) stop("Expression matrix must contain column: gene_symbol")
  gene <- df$gene_symbol |> as.character()
  mat <- as.matrix(df[, setdiff(colnames(df), "gene_symbol"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene
  mat
}

read_metadata <- function(path) {
  readr::read_tsv(path, col_types = cols(), progress = FALSE) |>
    mutate(
      sample_id = as.character(sample_id),
      group = as.character(group),
      phase = as.character(phase)
    )
}

zscore_by_gene <- function(mat) {
  # mat: genes x samples
  out <- mat
  for (i in seq_len(nrow(out))) {
    v <- out[i, ]
    s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      out[i, ] <- 0
    } else {
      out[i, ] <- (v - mean(v, na.rm = TRUE)) / s
    }
  }
  out
}

compute_cdf_auc <- function(consensus_matrix, step = 0.01) {
  vals <- consensus_matrix[upper.tri(consensus_matrix)]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(list(x = numeric(0), y = numeric(0), auc = NA_real_))
  f <- stats::ecdf(vals)
  x <- seq(0, 1, by = step)
  y <- f(x)
  auc <- sum(y) * step
  list(x = x, y = y, auc = auc)
}

plot_consensus_heatmap <- function(cm, classes, meta, out_pdf) {
  stopifnot(length(classes) == nrow(cm), nrow(cm) == ncol(cm))
  if (is.null(rownames(cm)) || is.null(colnames(cm))) {
    if (!is.null(names(classes)) && length(names(classes)) == nrow(cm)) {
      rownames(cm) <- names(classes)
      colnames(cm) <- names(classes)
    }
  }
  ord <- order(classes, names(classes))
  cm2 <- cm[ord, ord, drop = FALSE]
  classes2 <- classes[ord]

  clus_levels <- sort(unique(classes2))
  n_clus <- length(clus_levels)
  # SOTA: Dark2 is good for categorical clusters
  pal_cols <- grDevices::hcl.colors(n_clus, "Dark 2")
  clus_cols <- setNames(pal_cols[seq_len(n_clus)], as.character(clus_levels))

  phase_map <- meta |>
    filter(sample_id %in% names(classes2)) |>
    distinct(sample_id, phase)
  phase_vec <- setNames(phase_map$phase, phase_map$sample_id)
  phase_vec <- phase_vec[names(classes2)]
  phase_vec[is.na(phase_vec)] <- "unknown"
  phase_levels <- sort(unique(phase_vec))

  anno <- data.frame(
    cluster = factor(as.character(classes2), levels = as.character(clus_levels)),
    phase = factor(phase_vec, levels = phase_levels)
  )
  rownames(anno) <- names(classes2)

  phase_cols <- NULL
  if (length(phase_levels) > 0) {
    pal_phase <- grDevices::hcl.colors(max(3, length(phase_levels)), "Dark 3")
    phase_cols <- setNames(pal_phase[seq_len(length(phase_levels))], phase_levels)
  }

  anno_colors <- list(cluster = clus_cols)
  if (!is.null(phase_cols)) {
    anno_colors$phase <- phase_cols
  }

  # SOTA: Consensus Matrix should be mono-hue (White -> Dark Blue) to represent probability 0-1
  heat_cols <- grDevices::colorRampPalette(c("white", "#08306B"))(100)

  pheatmap::pheatmap(
    cm2,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_row = anno,
    annotation_col = anno,
    annotation_colors = anno_colors,
    color = heat_cols,
    border_color = NA,
    legend = TRUE,
    main = "",
    filename = out_pdf,
    width = 7,
    height = 6
  )
  invisible(NULL)
}

genes <- read_gene_symbols(gene_list_path)
meta <- read_metadata(meta_path)
mat <- read_expression_matrix(expr_path)

case_samples <- meta |>
  filter(group == case_label) |>
  pull(sample_id) |>
  unique()

common_samples <- intersect(colnames(mat), case_samples)
if (length(common_samples) < 10) stop("Too few case samples after matching expression matrix columns (need >=10).")

genes_present <- intersect(genes, rownames(mat))
if (length(genes_present) < 5) stop("Too few genes present in expression matrix after filtering (need >=5).")
genes_missing <- setdiff(genes, genes_present)

if (dry_run) {
  cat(paste0(
    "Dry run OK.\n",
    "- expr_matrix: ", expr_path, "\n",
    "- metadata: ", meta_path, "\n",
    "- gene_list: ", gene_list_path, "\n",
    "- total_genes_requested: ", length(genes), "\n",
    "- genes_present: ", length(genes_present), "\n",
    "- genes_missing: ", length(genes_missing), "\n",
    "- case_samples: ", length(case_samples), "\n",
    "- matched_case_samples: ", length(common_samples), "\n",
    "- max_k: ", max_k, "\n",
    "- reps: ", reps, "\n"
  ))
  quit(status = 0)
}

expr_case <- mat[genes_present, common_samples, drop = FALSE]
expr_case <- zscore_by_gene(expr_case)

cc <- ConsensusClusterPlus(
  d = expr_case,
  maxK = max_k,
  reps = reps,
  pItem = p_item,
  pFeature = p_feature,
  clusterAlg = "hc",
  distance = "pearson",
  seed = seed,
  plot = NULL,
  title = tempdir()
)

ks <- 2:max_k
cdf_rows <- lapply(ks, function(k) {
  cm <- cc[[k]]$consensusMatrix
  cdf <- compute_cdf_auc(cm)
  tibble::tibble(k = k, auc = cdf$auc)
})
cdf_df <- bind_rows(cdf_rows) |>
  arrange(k) |>
  mutate(delta_auc = auc - dplyr::lag(auc))

best_k <- cdf_df |>
  filter(!is.na(delta_auc)) |>
  arrange(desc(delta_auc), k) |>
  slice(1) |>
  pull(k)
best_k <- as.integer(best_k)

classes <- cc[[best_k]]$consensusClass
classes <- setNames(as.integer(classes), names(classes))

assignments <- tibble::tibble(
  sample_id = names(classes),
  cluster = unname(classes)
) |>
  left_join(meta |> select(sample_id, group, phase), by = "sample_id") |>
  arrange(cluster, sample_id)

summary <- cdf_df |>
  mutate(selected_k = (k == best_k)) |>
  mutate(
    n_samples = length(common_samples),
    n_genes = length(genes_present),
    gene_list_path = gene_list_path,
    expr_matrix_path = expr_path,
    metadata_path = meta_path,
    seed = seed,
    reps = reps,
    p_item = p_item,
    p_feature = p_feature
  )

readr::write_csv(assignments, file.path(out_dir_tables, out_assignments))
readr::write_csv(summary, file.path(out_dir_tables, out_summary))

# Fig6B: CDF curves (SOTA: theme_bw)
cdf_plot_rows <- lapply(ks, function(k) {
  cm <- cc[[k]]$consensusMatrix
  cdf <- compute_cdf_auc(cm)
  tibble::tibble(k = k, x = cdf$x, y = cdf$y)
})
cdf_plot_df <- bind_rows(cdf_plot_rows) |>
  mutate(
    k = factor(k, levels = ks),
    selected = (as.integer(as.character(k)) == best_k)
  )

p_cdf <- ggplot(cdf_plot_df, aes(x = x, y = y, group = k)) +
  geom_line(aes(color = k, linewidth = selected, alpha = selected)) +
  scale_linewidth_manual(values = c(`FALSE` = 0.7, `TRUE` = 1.3), guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.75, `TRUE` = 1.0), guide = "none") +
  scale_color_manual(values = setNames(grDevices::hcl.colors(length(ks), "Dark 2"), levels(cdf_plot_df$k))) +
  labs(
    x = "Consensus Index",
    y = "CDF",
    subtitle = sprintf("Selected k = %d (Max Delta AUC)", best_k),
    color = "k"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_dir_figs, out_fig_cdf),
  plot = p_cdf,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 300
)

# Fig6A: consensus matrix heatmap (selected k)
plot_consensus_heatmap(
  cm = cc[[best_k]]$consensusMatrix,
  classes = classes,
  meta = meta,
  out_pdf = file.path(out_dir_figs, out_fig_heatmap)
)

plot_subtype_expression_heatmap <- function(expr_case, classes, meta, out_png) {
  stopifnot(ncol(expr_case) == length(classes))
  clus <- factor(classes)
  ord <- order(clus)
  expr_ord <- expr_case[, ord, drop = FALSE]
  clus_ord <- clus[ord]

  phase_map <- meta |>
    filter(sample_id %in% colnames(expr_ord)) |>
    distinct(sample_id, phase)
  phase_vec <- setNames(phase_map$phase, phase_map$sample_id)
  phase_vec <- phase_vec[colnames(expr_ord)]
  phase_vec[is.na(phase_vec)] <- "unknown"

  clus_levels <- sort(unique(clus_ord))
  pal_cols <- grDevices::hcl.colors(length(clus_levels), "Dark 2")
  clus_cols <- setNames(pal_cols, as.character(clus_levels))

  phase_levels <- sort(unique(phase_vec))
  pal_phase <- grDevices::hcl.colors(max(3, length(phase_levels)), "Dark 3")
  phase_cols <- setNames(pal_phase[seq_len(length(phase_levels))], phase_levels)

  anno <- data.frame(
    cluster = factor(clus_ord, levels = clus_levels),
    phase = factor(phase_vec, levels = phase_levels)
  )
  rownames(anno) <- colnames(expr_ord)

  # SOTA: Navy-White-Firebrick for Z-scores
  heat_cols <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)
  
  pheatmap::pheatmap(
    expr_ord,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = FALSE, # Clean look
    show_colnames = FALSE,
    annotation_col = anno,
    annotation_colors = list(cluster = clus_cols, phase = phase_cols),
    color = heat_cols,
    border_color = NA,
    main = "",
    filename = out_png,
    width = 8,
    height = 6
  )
}

plot_subtype_marker_panels <- function(expr_case, classes, out_png, n_markers = 6) {
  clus <- factor(classes)
  pvals <- sapply(seq_len(nrow(expr_case)), function(i) {
    df <- data.frame(expr = as.numeric(expr_case[i, ]), cluster = clus)
    stats::anova(stats::aov(expr ~ cluster, data = df))[["Pr(>F)"]][1]
  })
  names(pvals) <- rownames(expr_case)
  pvals <- sort(pvals)
  top_genes <- names(pvals)[seq_len(min(n_markers, length(pvals)))]

  df_long <- data.frame(
    gene_symbol = rep(top_genes, each = ncol(expr_case)),
    expr = as.numeric(expr_case[top_genes, , drop = FALSE]),
    cluster = rep(as.character(clus), times = length(top_genes))
  )

  pal_cols <- grDevices::hcl.colors(length(unique(df_long$cluster)), "Dark 2")
  names(pal_cols) <- sort(unique(df_long$cluster))

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = cluster, y = expr, fill = cluster)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.2, width = 0.8, color = NA) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.4, linewidth = 0.4, color = "black") +
    ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, shape = 21, color = "white", stroke = 0.2) +
    ggplot2::facet_wrap(~gene_symbol, scales = "free_y", ncol = 3) +
    ggplot2::scale_fill_manual(values = pal_cols) +
    ggplot2::labs(x = "Subtype (Cluster)", y = "Z-scored Expression") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(filename = out_png, plot = p, width = 8, height = 6, units = "in", dpi = 300)
}

# Fig6C: subtype expression heatmap (cases only)
plot_subtype_expression_heatmap(
  expr_case = expr_case,
  classes = classes,
  meta = meta,
  out_png = file.path(out_dir_figs, out_fig_subtype_heatmap)
)

# Fig6D: subtype marker panels (cases only)
plot_subtype_marker_panels(
  expr_case = expr_case,
  classes = classes,
  out_png = file.path(out_dir_figs, out_fig_subtype_markers),
  n_markers = 6
)

cat(paste0(
  "Consensus clustering completed.\n",
  "- best_k: ", best_k, "\n",
  "- cases: ", length(common_samples), "\n",
  "- genes: ", length(genes_present), "\n",
  "- assignments: ", file.path(out_dir_tables, out_assignments), "\n",
  "- summary: ", file.path(out_dir_tables, out_summary), "\n",
  "- Fig6A: ", file.path(out_dir_figs, out_fig_heatmap), "\n",
  "- Fig6B: ", file.path(out_dir_figs, out_fig_cdf), "\n",
  "- Fig6C: ", file.path(out_dir_figs, out_fig_subtype_heatmap), "\n",
  "- Fig6D: ", file.path(out_dir_figs, out_fig_subtype_markers), "\n"
))
