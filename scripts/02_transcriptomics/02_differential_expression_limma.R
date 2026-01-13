#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(sprintf("Missing value for %s", flag))
  args[[idx + 1]]
}

input_matrix <- get_arg("--matrix", "data/processed/expression_matrix_combat.tsv")
input_meta <- get_arg("--meta", "data/processed/sample_metadata.tsv")
out_dir_tables <- get_arg("--out-tables", "results/tables")
out_dir_figs <- get_arg("--out-figs", "figures/raw_plots")
logfc_cutoff <- as.numeric(get_arg("--logfc", "0.585"))
adjp_cutoff <- as.numeric(get_arg("--adjp", "0.05"))
block_col <- get_arg("--block-col", "")
covariates_raw <- get_arg("--covariates", "")
include_phases_raw <- get_arg("--include-phases", "")
phase_col <- get_arg("--phase-col", "phase")
expr_unit <- get_arg("--expr-unit", "Normalized Expression (log2)")
id_map_path <- get_arg("--id-map", "data/processed/human_id_to_symbol.tsv")

deg_name <- get_arg("--deg-name", "T2.2_deg_results.csv")
volcano_name <- get_arg("--volcano-name", "Fig3A.png")
heatmap_name <- get_arg("--heatmap-name", "Fig3B.png")
boxplot_name <- get_arg("--boxplot-name", "Fig3C.png")

if (!file.exists(input_matrix)) stop(sprintf("Matrix not found: %s", input_matrix))
if (!file.exists(input_meta)) stop(sprintf("Metadata not found: %s", input_meta))

dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figs, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

id_map <- NULL
if (!is.null(id_map_path) && nzchar(id_map_path) && file.exists(id_map_path)) {
  id_map <- readr::read_tsv(
    id_map_path,
    col_names = c("ensembl_id", "gene_symbol_mapped"),
    col_types = cols(),
    progress = FALSE
  )
}

strip_ensembl_version <- function(x) {
  sub("\\\\.[0-9]+$", "", x)
}

map_gene_labels <- function(ids) {
  if (is.null(id_map)) return(ids)
  ids_in <- strip_ensembl_version(ids)
  out <- ids
  ensg_idx <- grepl("^ENSG", ids_in)
  if (!any(ensg_idx)) return(ids)
  m <- match(ids_in[ensg_idx], id_map$ensembl_id)
  syms <- id_map$gene_symbol_mapped[m]
  missing <- is.na(syms) | syms == ""
  syms[missing] <- ids[ensg_idx][missing]
  out[ensg_idx] <- syms
  out
}

make_unique_labels <- function(labels, ids) {
  dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[dup] <- paste0(labels[dup], " (", ids[dup], ")")
  labels
}

expr <- readr::read_tsv(input_matrix, col_types = cols(), progress = FALSE)
expr <- as.data.frame(expr)
gene_col <- colnames(expr)[[1]]
rownames(expr) <- expr[[gene_col]]
expr[[gene_col]] <- NULL
expr <- as.matrix(expr)

meta <- readr::read_tsv(input_meta, col_types = cols(), progress = FALSE)
meta <- as.data.frame(meta)
if (!all(c("sample_id", "group") %in% colnames(meta))) {
  stop("Metadata must contain columns: sample_id, group")
}

missing_meta <- setdiff(colnames(expr), meta$sample_id)
if (length(missing_meta) > 0) {
  stop(sprintf("Metadata missing %s sample(s) from matrix (first 5: %s)", length(missing_meta), paste(head(missing_meta, 5), collapse = ", ")))
}

meta <- meta[match(colnames(expr), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(expr)))

if (!is.null(include_phases_raw) && nzchar(include_phases_raw)) {
  if (!phase_col %in% colnames(meta)) {
    stop(sprintf("Phase filter requested but metadata is missing column: %s", phase_col))
  }
  phases <- unlist(strsplit(include_phases_raw, ",", fixed = TRUE))
  phases <- trimws(phases)
  phases <- phases[nzchar(phases)]
  if (length(phases) == 0) stop("include-phases provided but no phases parsed.")
  keep <- meta[[phase_col]] %in% phases
  if (!any(keep)) stop(sprintf("No samples matched include-phases=%s using column %s.", include_phases_raw, phase_col))
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  meta <- meta[match(colnames(expr), meta$sample_id), , drop = FALSE]
  stopifnot(all(meta$sample_id == colnames(expr)))
  message(sprintf("Phase filter: %s (column=%s) -> N=%d samples", paste(phases, collapse = ","), phase_col, nrow(meta)))
}

group <- factor(meta$group, levels = c("control", "case"))
if (any(is.na(group))) stop("Group must be 'control' or 'case' for all samples.")
meta$group <- group

covariates <- character(0)
if (!is.null(covariates_raw) && nzchar(covariates_raw)) {
  covariates <- unlist(strsplit(covariates_raw, ",", fixed = TRUE))
  covariates <- trimws(covariates)
  covariates <- covariates[nzchar(covariates)]
}

for (cov in covariates) {
  if (!cov %in% colnames(meta)) {
    stop(sprintf("Covariate not found in metadata: %s", cov))
  }
  if (!is.numeric(meta[[cov]])) {
    meta[[cov]] <- factor(meta[[cov]])
  }
}

terms <- c()
if ("batch" %in% colnames(meta) && length(unique(meta$batch)) > 1) {
  batch <- factor(meta$batch)
  meta$batch <- batch
  terms <- c(terms, "batch")
}
terms <- c(terms, covariates, "group")
terms <- unique(terms)

design_formula <- as.formula(paste("~", paste(terms, collapse = " + ")))
message("Design formula: ", deparse(design_formula))
design <- model.matrix(design_formula, data = meta)
coef_name <- paste0("group", levels(meta$group)[2])

block <- NULL
if (!is.null(block_col) && nzchar(block_col)) {
  if (!block_col %in% colnames(meta)) {
    stop(sprintf("block column not found in metadata: %s", block_col))
  }
  block <- factor(meta[[block_col]])
  if (length(unique(block)) >= nrow(meta)) {
    message(sprintf("Block column %s has no repeats; ignoring block.", block_col))
    block <- NULL
  }
}

if (!is.null(block)) {
  corfit <- duplicateCorrelation(expr, design, block = block)
  message(sprintf("duplicateCorrelation consensus correlation: %.4f", corfit$consensus))
  fit <- lmFit(expr, design, block = block, correlation = corfit$consensus)
} else {
  fit <- lmFit(expr, design)
}
fit <- eBayes(fit)

deg <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
deg <- deg %>%
  rownames_to_column("gene_symbol") %>%
  as_tibble()

deg_path <- file.path(out_dir_tables, deg_name)
readr::write_csv(deg, deg_path)

deg_filtered <- deg %>%
  filter(!is.na(adj.P.Val)) %>%
  mutate(sig = (abs(logFC) > logfc_cutoff) & (adj.P.Val < adjp_cutoff))
deg_filtered <- deg_filtered %>%
  mutate(gene_label = map_gene_labels(gene_symbol))

case_label <- "Case" # Define case label for plotting
volcano_path <- file.path(out_dir_figs, volcano_name)
deg_filtered <- deg_filtered %>%
  mutate(neglog10_adjP = -log10(pmax(adj.P.Val, 1e-300)))
deg_filtered$status <- ifelse(deg_filtered$sig & deg_filtered$logFC > 0, "Up",
                              ifelse(deg_filtered$sig & deg_filtered$logFC < 0, "Down", "Not Sig"))
sig_n <- deg_filtered %>% filter(sig) %>% nrow()

# SOTA Volcano Plot (EnhancedVolcano Style)
p_volcano <- ggplot(deg_filtered, aes(x = logFC, y = neglog10_adjP)) +
  geom_point(aes(color = status), alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c("Up" = "firebrick3", "Down" = "navy", "Not Sig" = "grey85"),
    labels = c("Up" = "Up-regulated", "Down" = "Down-regulated", "Not Sig" = "NS")
  ) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(adjp_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = sprintf("%s vs Control (n=%d)", case_label, nrow(deg_filtered)),
    subtitle = sprintf("Sig: %d (|logFC|>%.2f, FDR<%.2f)", sig_n, logfc_cutoff, adjp_cutoff),
    x = bquote(log[2] ~ "Fold Change"),
    y = bquote(-log[10] ~ "FDR")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Label top genes (SOTA: use ggrepel)
if ("gene_label" %in% colnames(deg_filtered)) {
  # Label top 10 significant genes, EXCLUDING those that failed mapping (still ENSG)
  top_labels <- deg_filtered %>%
    filter(status != "Not Sig") %>%
    filter(!grepl("^ENSG", gene_label)) %>% # HARD FILTER: No ENSG IDs in labels
    arrange(adj.P.Val) %>%
    head(10)
  
  if (nrow(top_labels) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p_volcano <- p_volcano +
        ggrepel::geom_text_repel(
          data = top_labels,
          aes(label = gene_label),
          size = 3,
          box.padding = 0.5,
          max.overlaps = Inf
        )
    } else {
      p_volcano <- p_volcano + geom_text(data = top_labels, aes(label = gene_label), vjust = 1.5)
    }
  }
}

ggsave(volcano_path, p_volcano, width = 6, height = 5, units = "in", dpi = 300)

heatmap_path <- file.path(out_dir_figs, heatmap_name)
top_n <- 50
deg_ranked <- deg_filtered %>%
  arrange(adj.P.Val) %>%
  filter(is.finite(adj.P.Val))
deg_ranked_pref <- deg_ranked %>% filter(!grepl("^ENSG", gene_label))
if (nrow(deg_ranked_pref) >= top_n) {
  top_genes <- deg_ranked_pref %>% slice_head(n = top_n) %>% pull(gene_symbol)
} else {
  top_genes <- deg_ranked %>% slice_head(n = top_n) %>% pull(gene_symbol)
}
top_labels <- make_unique_labels(map_gene_labels(top_genes), top_genes)

heat_mat <- expr[top_genes, , drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))
heat_mat[is.na(heat_mat)] <- 0

anno <- data.frame(Group = group)
rownames(anno) <- meta$sample_id
# SOTA Colors: Navy/White/Firebrick for Z-scores
anno_colors <- list(Group = c(control = "#4C72B0", case = "#C44E52"))
heat_colors <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(
  heat_mat,
  annotation_col = anno,
  annotation_colors = anno_colors,
  color = heat_colors,
  show_colnames = FALSE,
  labels_row = top_labels,
  fontsize_row = 8,
  clustering_method = "ward.D2", # Better clusters than 'complete' usually
  border_color = NA,
  main = sprintf("Top %d DEGs (Z-score)", length(top_genes)),
  filename = heatmap_path,
  width = 8,
  height = 6
)

boxplot_path <- file.path(out_dir_figs, boxplot_name)
box_pool <- deg_filtered %>% filter(sig)
if (nrow(box_pool) < 1) {
  box_pool <- deg_filtered
}
box_pool <- box_pool %>% arrange(adj.P.Val)
box_pool_pref <- box_pool %>% filter(!grepl("^ENSG", gene_label))
if (nrow(box_pool_pref) >= 6) {
  box_genes <- box_pool_pref %>% slice_head(n = 6) %>% pull(gene_symbol)
} else {
  box_genes <- box_pool %>% slice_head(n = 6) %>% pull(gene_symbol)
}
box_labels <- make_unique_labels(map_gene_labels(box_genes), box_genes)
box_label_map <- tibble(gene_symbol = box_genes, gene_label = box_labels)

long_df <- as.data.frame(t(expr[box_genes, , drop = FALSE])) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta %>% select(any_of(c("sample_id", "group", "batch"))), by = "sample_id") %>%
  pivot_longer(cols = all_of(box_genes), names_to = "gene_symbol", values_to = "expression") %>%
  left_join(box_label_map, by = "gene_symbol")

p_box <- ggplot(long_df, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.2, width = 0.8, color = NA) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.4, linewidth = 0.4, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, shape = 21, color = "white", stroke = 0.2) +
  facet_wrap(~gene_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(control = "#4C72B0", case = "#C44E52")) +
  labs(x = NULL, y = expr_unit) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

ggsave(boxplot_path, p_box, width = 8, height = 5, units = "in", dpi = 300)

message("Saved DEG table: ", deg_path)
message("Saved volcano:   ", volcano_path)
message("Saved heatmap:   ", heatmap_path)
message("Saved boxplots:  ", boxplot_path)
message("\nSession info:")
print(sessionInfo())
