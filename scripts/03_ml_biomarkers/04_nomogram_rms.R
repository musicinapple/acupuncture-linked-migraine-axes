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

input_train_matrix <- get_arg("--train-matrix", "data/processed/human_expression_matrix.tsv")
input_train_meta <- get_arg("--train-meta", "data/processed/human_sample_metadata.tsv")
input_external_matrix <- get_arg("--external-matrix", "data/processed/validation_PRJEB40032_expression_matrix.tsv")
input_external_meta <- get_arg("--external-meta", "data/processed/validation_PRJEB40032_sample_metadata.tsv")

signature_genes_path <- get_arg("--signature", "results/tables/T3.3_human_intersection_genes.csv")
deg_table_path <- get_arg("--deg-table", "results/tables/T2.6_human_deg_emtab13397.csv")

out_dir_tables <- get_arg("--out-tables", "results/tables")
out_dir_figs <- get_arg("--out-figs", "figures/raw_plots")

out_metrics_name <- get_arg("--metrics-name", "T3.2_nomogram_metrics.csv")
out_config_name <- get_arg("--config-name", "T3.2_nomogram_config.json")

roc_name <- get_arg("--roc-name", "Fig4E_nomogram_ROC.pdf")
nomogram_name <- get_arg("--nomogram-name", "Fig4F_nomogram.pdf")
calibration_name <- get_arg("--calibration-name", "Fig4G_calibration.pdf")
dca_name <- get_arg("--dca-name", "Fig4H_dca.pdf")

bootstrap_B <- as.integer(get_arg("--bootstrap", "200"))
dca_step <- as.numeric(get_arg("--dca-step", "0.01"))
seed <- as.integer(get_arg("--seed", "42"))

if (!file.exists(input_train_matrix)) stop(sprintf("Train matrix not found: %s", input_train_matrix))
if (!file.exists(input_train_meta)) stop(sprintf("Train meta not found: %s", input_train_meta))
if (!file.exists(signature_genes_path)) stop(sprintf("Signature genes file not found: %s", signature_genes_path))
if (!file.exists(deg_table_path)) stop(sprintf("DEG table not found: %s", deg_table_path))
if (!file.exists(input_external_matrix)) stop(sprintf("External matrix not found: %s", input_external_matrix))
if (!file.exists(input_external_meta)) stop(sprintf("External meta not found: %s", input_external_meta))

dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figs, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(rms)
  library(ggplot2)
})

utc_now_iso <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}

read_expression <- function(path) {
  expr <- readr::read_tsv(path, col_types = cols(), progress = FALSE) %>% as.data.frame()
  gene_col <- colnames(expr)[[1]]
  rownames(expr) <- expr[[gene_col]]
  expr[[gene_col]] <- NULL
  expr <- as.matrix(expr)
  list(expr = expr, gene_col = gene_col)
}

read_meta <- function(path) {
  meta <- readr::read_tsv(path, col_types = cols(), progress = FALSE) %>% as.data.frame()
  if (!all(c("sample_id", "group") %in% colnames(meta))) {
    stop("Metadata must contain columns: sample_id, group")
  }
  meta$sample_id <- as.character(meta$sample_id)
  meta$group <- as.character(meta$group)
  meta
}

read_signature_genes <- function(path) {
  if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) {
    df <- readr::read_tsv(path, col_types = cols(), progress = FALSE) %>% as.data.frame()
  } else {
    df <- readr::read_csv(path, col_types = cols(), progress = FALSE) %>% as.data.frame()
  }
  if (!("gene_symbol" %in% colnames(df))) stop("Signature file must contain column: gene_symbol")
  genes <- df$gene_symbol %>% as.character() %>% trimws()
  genes <- genes[genes != ""]
  unique(genes)
}

read_deg_signs <- function(path) {
  deg <- readr::read_csv(path, col_types = cols(), progress = FALSE) %>% as.data.frame()
  if (!all(c("gene_symbol", "logFC") %in% colnames(deg))) stop("DEG table must contain columns: gene_symbol, logFC")
  deg <- deg %>% filter(!is.na(gene_symbol), !is.na(logFC))
  signs <- ifelse(deg$logFC >= 0, 1, -1)
  stats::setNames(signs, as.character(deg$gene_symbol))
}

compute_signature_score <- function(expr, genes, signs) {
  genes_present <- intersect(genes, rownames(expr))
  if (length(genes_present) == 0) stop("No signature genes present in expression matrix.")

  m <- expr[genes_present, , drop = FALSE]
  # Per-gene z-score across samples, preserving sample names.
  m <- t(scale(t(m)))
  m[is.na(m)] <- 0
  colnames(m) <- colnames(expr)
  rownames(m) <- genes_present

  gsign <- sapply(genes_present, function(g) {
    if (!is.null(signs[[g]])) signs[[g]] else 1
  })
  score <- colMeans(m * gsign, na.rm = TRUE)
  names(score) <- colnames(expr)

  missing <- setdiff(genes, genes_present)
  list(score = score, genes_present = genes_present, genes_missing = missing)
}

auc_rank <- function(y, score) {
  # y: 0/1, score numeric; compute AUC via Mann-Whitney U (ties handled by average ranks)
  y <- as.integer(y)
  score <- as.numeric(score)
  r <- rank(score, ties.method = "average")
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  sum_r_pos <- sum(r[y == 1])
  u <- sum_r_pos - n_pos * (n_pos + 1) / 2
  u / (n_pos * n_neg)
}

roc_curve_points <- function(y, score) {
  y <- as.integer(y)
  score <- as.numeric(score)
  th <- sort(unique(score), decreasing = TRUE)
  th <- c(Inf, th, -Inf)
  tpr <- numeric(length(th))
  fpr <- numeric(length(th))

  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  for (i in seq_along(th)) {
    pred_pos <- score >= th[[i]]
    tp <- sum(pred_pos & (y == 1))
    fp <- sum(pred_pos & (y == 0))
    tpr[[i]] <- if (n_pos == 0) NA_real_ else tp / n_pos
    fpr[[i]] <- if (n_neg == 0) NA_real_ else fp / n_neg
  }
  data.frame(threshold = th, fpr = fpr, tpr = tpr)
}

decision_curve <- function(y, prob, thresholds) {
  y <- as.integer(y)
  prob <- as.numeric(prob)
  n <- length(y)
  prev <- mean(y == 1)
  out <- lapply(thresholds, function(pt) {
    pred_pos <- prob >= pt
    tp <- sum(pred_pos & (y == 1))
    fp <- sum(pred_pos & (y == 0))
    w <- pt / (1 - pt)
    nb_model <- (tp / n) - (fp / n) * w
    nb_all <- prev - (1 - prev) * w
    nb_none <- 0
    data.frame(threshold = pt, nb_model = nb_model, nb_all = nb_all, nb_none = nb_none)
  })
  dplyr::bind_rows(out)
}

set.seed(seed)

sig_genes <- read_signature_genes(signature_genes_path)
deg_signs <- read_deg_signs(deg_table_path)

train <- read_expression(input_train_matrix)
train_expr <- train$expr
train_meta <- read_meta(input_train_meta)

missing_meta <- setdiff(colnames(train_expr), train_meta$sample_id)
if (length(missing_meta) > 0) {
  stop(sprintf("Train metadata missing %s sample(s) from matrix (first 5: %s)", length(missing_meta), paste(head(missing_meta, 5), collapse = ", ")))
}
train_meta <- train_meta[match(colnames(train_expr), train_meta$sample_id), ]

train_score_res <- compute_signature_score(train_expr, sig_genes, deg_signs)
train_df <- train_meta %>%
  mutate(
    signature_score = as.numeric(train_score_res$score[sample_id]),
    y = ifelse(group == "case", 1, 0)
  )

if (any(is.na(train_df$signature_score))) stop("NA signature_score in training set after alignment.")

dd <- rms::datadist(train_df)
options(datadist = "dd")
fit <- rms::lrm(y ~ signature_score, data = train_df, x = TRUE, y = TRUE)

train_prob <- stats::predict(fit, type = "fitted")
train_auc <- auc_rank(train_df$y, train_prob)

external <- read_expression(input_external_matrix)
ext_expr <- external$expr
ext_meta <- read_meta(input_external_meta)
missing_meta_ext <- setdiff(colnames(ext_expr), ext_meta$sample_id)
if (length(missing_meta_ext) > 0) {
  stop(sprintf("External metadata missing %s sample(s) from matrix (first 5: %s)", length(missing_meta_ext), paste(head(missing_meta_ext, 5), collapse = ", ")))
}
ext_meta <- ext_meta[match(colnames(ext_expr), ext_meta$sample_id), ]

ext_score_res <- compute_signature_score(ext_expr, sig_genes, deg_signs)
ext_df <- ext_meta %>%
  mutate(
    signature_score = as.numeric(ext_score_res$score[sample_id]),
    y = ifelse(group == "case", 1, 0)
  )
if (any(is.na(ext_df$signature_score))) stop("NA signature_score in external set after alignment.")

ext_prob <- stats::predict(fit, newdata = ext_df, type = "fitted")
ext_auc <- auc_rank(ext_df$y, ext_prob)

roc_train <- roc_curve_points(train_df$y, train_prob) %>% mutate(dataset = "Train (E-MTAB-13397)")
roc_ext <- roc_curve_points(ext_df$y, ext_prob) %>% mutate(dataset = "External (PRJEB40032)")
roc_plot_df <- bind_rows(roc_train, roc_ext) %>% filter(!is.na(fpr), !is.na(tpr))

roc_title <- "Signature Performance (ROC)"
# roc_subtitle <- sprintf("Train AUC=%.3f (n=%d) | External AUC=%.3f (n=%d)", 
#                        train_auc, nrow(train_df), ext_auc, nrow(ext_df))

roc_plot <- ggplot(roc_plot_df, aes(x = fpr, y = tpr, color = dataset)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("Train (E-MTAB-13397)" = "#E41A1C", "External (PRJEB40032)" = "#377EB8")) +
  coord_equal() +
  labs(title = roc_title, x = "False Positive Rate", y = "True Positive Rate") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.7, 0.15),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Add distinct AUC label using geom_label to simulate bbox
auc_label_df <- data.frame(
  fpr = 0.5, tpr = 0.1,
  label = sprintf("Train AUC: %.3f\nExternal AUC: %.3f", train_auc, ext_auc)
)
roc_plot <- roc_plot + geom_label(
  data = auc_label_df, 
  aes(x = fpr, y = tpr, label = label), 
  inherit.aes = FALSE,
  color = "black", 
  fill = "white", 
  alpha = 0.9,
  fontface = "bold",
  size = 3.5,
  hjust = 0
)

ggsave(file.path(out_dir_figs, roc_name), roc_plot, width = 6, height = 6)

pdf(file.path(out_dir_figs, nomogram_name), width = 9, height = 6.5)
par(mar = c(5, 5, 4, 2), cex = 0.9)
nom <- rms::nomogram(fit, fun = plogis, funlabel = "Migraine probability")
plot(nom, xfrac = 0.4)
dev.off()

pdf(file.path(out_dir_figs, calibration_name), width = 6, height = 6)
par(mar = c(5, 5, 4, 2), cex.lab = 1.1, cex.axis = 1.0, las = 1)
cal <- rms::calibrate(fit, method = "boot", B = bootstrap_B)
plot(cal, xlab = "Predicted Probability", ylab = "Observed Proportion", 
     subtitles = FALSE, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lty = 2, col = "gray50", lwd = 1.5)
title(sprintf("Calibration (B=%d)", bootstrap_B), adj = 0)
grid(col = "gray90")
dev.off()

thresholds <- seq(dca_step, 1 - dca_step, by = dca_step)
dca_train <- decision_curve(train_df$y, train_prob, thresholds) %>% mutate(dataset = "Train")
dca_ext <- decision_curve(ext_df$y, ext_prob, thresholds) %>% mutate(dataset = "External")
dca_plot_df <- bind_rows(dca_train, dca_ext) %>%
  pivot_longer(cols = c(nb_model, nb_all, nb_none), names_to = "strategy", values_to = "net_benefit")

strategy_label <- c(nb_model = "Model", nb_all = "Treat all", nb_none = "Treat none")
dca_plot_df$strategy <- factor(dca_plot_df$strategy, levels = names(strategy_label), labels = strategy_label)

dca_plot <- ggplot(dca_plot_df, aes(x = threshold, y = net_benefit, color = strategy, linetype = dataset)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Decision Curve Analysis (Net Benefit)",
    x = "Threshold Probability",
    y = "Net Benefit"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir_figs, dca_name), dca_plot, width = 7, height = 5)

metrics <- tibble::tibble(
  dataset = c("Train (E-MTAB-13397)", "External (PRJEB40032)"),
  n_total = c(nrow(train_df), nrow(ext_df)),
  n_cases = c(sum(train_df$y == 1), sum(ext_df$y == 1)),
  n_controls = c(sum(train_df$y == 0), sum(ext_df$y == 0)),
  auc = c(train_auc, ext_auc),
  signature_genes_total = length(sig_genes),
  signature_genes_used_train = length(train_score_res$genes_present),
  signature_genes_missing_train = length(train_score_res$genes_missing),
  signature_genes_used_external = length(ext_score_res$genes_present),
  signature_genes_missing_external = length(ext_score_res$genes_missing)
)
readr::write_csv(metrics, file.path(out_dir_tables, out_metrics_name))

config <- list(
  created_at = utc_now_iso(),
  inputs = list(
    train_matrix = input_train_matrix,
    train_meta = input_train_meta,
    external_matrix = input_external_matrix,
    external_meta = input_external_meta,
    signature = signature_genes_path,
    deg_table = deg_table_path
  ),
  params = list(
    bootstrap_B = bootstrap_B,
    dca_step = dca_step,
    seed = seed
  ),
  outputs = list(
    roc = file.path(out_dir_figs, roc_name),
    nomogram = file.path(out_dir_figs, nomogram_name),
    calibration = file.path(out_dir_figs, calibration_name),
    dca = file.path(out_dir_figs, dca_name),
    metrics = file.path(out_dir_tables, out_metrics_name)
  ),
  notes = list(
    "Model: logistic regression (rms::lrm) with a single predictor (direction-aware signature_score).",
    "Signature score: mean(sign * per-gene z-score), where sign is from the training DEG logFC table.",
    "External evaluation uses the same fixed signature and score definition; no tuning on external data."
  )
)
writeLines(jsonlite::toJSON(config, pretty = TRUE, auto_unbox = TRUE), file.path(out_dir_tables, out_config_name))

message("Done.")
