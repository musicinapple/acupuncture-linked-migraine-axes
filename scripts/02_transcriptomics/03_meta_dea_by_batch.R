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
out_table <- get_arg("--out-table", "results/tables/T2.3_meta_deg_results.csv")
out_fig <- get_arg("--out-fig", "figures/raw_plots/Fig3E_meta_volcano.pdf")
logfc_cutoff <- as.numeric(get_arg("--logfc", "0.585"))
adjp_cutoff <- as.numeric(get_arg("--adjp", "0.05"))

if (!file.exists(input_matrix)) stop(sprintf("Matrix not found: %s", input_matrix))
if (!file.exists(input_meta)) stop(sprintf("Metadata not found: %s", input_meta))

dir.create(dirname(out_table), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_fig), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

expr <- readr::read_tsv(input_matrix, col_types = cols(), progress = FALSE)
expr <- as.data.frame(expr)
gene_col <- colnames(expr)[[1]]
rownames(expr) <- expr[[gene_col]]
expr[[gene_col]] <- NULL
expr <- as.matrix(expr)

meta <- readr::read_tsv(input_meta, col_types = cols(), progress = FALSE) %>% as.data.frame()
req <- c("sample_id", "group", "batch")
if (!all(req %in% colnames(meta))) stop("Metadata must contain: sample_id, group, batch")

meta <- meta[match(colnames(expr), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(expr)))

meta$group <- factor(meta$group, levels = c("control", "case"))
if (any(is.na(meta$group))) stop("All samples must have group 'control' or 'case'.")

batches <- sort(unique(meta$batch))
per_batch <- list()

for (b in batches) {
  idx <- which(meta$batch == b)
  sub_meta <- meta[idx, , drop = FALSE]
  if (length(unique(sub_meta$group)) < 2) {
    next
  }
  sub_expr <- expr[, idx, drop = FALSE]
  design <- model.matrix(~ sub_meta$group)
  fit <- lmFit(sub_expr, design)
  fit <- eBayes(fit)
  coef_name <- "sub_meta$groupcase"
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
    rownames_to_column("gene_symbol") %>%
    as_tibble() %>%
    transmute(
      gene_symbol,
      logFC = logFC,
      P.Value = P.Value,
      n_samples = ncol(sub_expr),
      batch = b
    )
  per_batch[[b]] <- tt
}

if (length(per_batch) < 2) stop("Need at least 2 batches with both case/control for meta-DE.")

all_tt <- bind_rows(per_batch)

meta_deg <- all_tt %>%
  group_by(gene_symbol) %>%
  summarise(
    n_batches = n(),
    n_total = sum(n_samples),
    logFC_meta = weighted.mean(logFC, w = n_samples, na.rm = TRUE),
    fisher_stat = -2 * sum(log(pmax(P.Value, 1e-300))),
    fisher_p = pchisq(fisher_stat, df = 2 * n(), lower.tail = FALSE),
    .groups = "drop"
  ) %>%
  mutate(adj.P.Val = p.adjust(fisher_p, method = "BH")) %>%
  arrange(adj.P.Val)

readr::write_csv(meta_deg, out_table)

plot_df <- meta_deg %>%
  mutate(
    sig = (abs(logFC_meta) > logfc_cutoff) & (adj.P.Val < adjp_cutoff),
    neglog10_adjP = -log10(pmax(adj.P.Val, 1e-300))
  )

p <- ggplot(plot_df, aes(x = logFC_meta, y = neglog10_adjP)) +
  geom_point(aes(color = sig), alpha = 0.75, size = 1.2) +
  scale_color_manual(values = c(`FALSE` = "#7F7F7F", `TRUE` = "#C44E52")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(adjp_cutoff), linetype = "dashed", linewidth = 0.3) +
  labs(x = "Meta logFC (weighted mean)", y = "-log10(adj.P)", color = "Significant") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

ggsave(out_fig, p, width = 6.5, height = 5.0, units = "in")

message("Saved meta-DE table: ", out_table)
message("Saved meta volcano:  ", out_fig)

