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

input_genes <- get_arg("--genes", "results/tables/T3.3_human_intersection_genes.csv")
outcome_id <- get_arg("--outcome-id", "ukb-d-G43")
outcome_trait_arg <- get_arg("--outcome-trait", "")

opengwas_jwt <- trimws(get_arg("--opengwas-jwt", Sys.getenv("OPENGWAS_JWT", "")))
eqtl_source_pattern <- get_arg("--eqtl-source-pattern", "eQTLGen")
eqtl_mapping_path <- get_arg("--eqtl-mapping", "")
symbol_to_ensembl_path <- get_arg("--symbol-to-ensembl", "data/processed/human_id_to_symbol.tsv")
eqtl_id_prefix <- get_arg("--eqtl-id-prefix", "eqtl-a-")
eqtl_mode <- get_arg("--eqtl-mode", "prefix")
eqtl_tissue_regex <- get_arg("--eqtl-tissue-regex", "")
eqtl_match <- get_arg("--eqtl-match", "ensembl")
eqtl_candidates_out <- get_arg("--eqtl-candidates-out", "")
eqtl_max_hits_per_gene <- as.integer(get_arg("--eqtl-max-hits-per-gene", "10"))
eqtl_discover_only <- flag_present("--eqtl-discover-only")
gwasinfo_cache <- get_arg("--gwasinfo-cache", "")
eqtl_catalogue_dataset_id <- get_arg("--eqtl-catalogue-dataset-id", "")
eqtl_catalogue_metadata_path <- get_arg("--eqtl-catalogue-metadata", "data/raw/eqtl_catalogue/dataset_metadata.tsv")
eqtl_catalogue_api_base <- get_arg("--eqtl-catalogue-api-base", "https://www.ebi.ac.uk/eqtl/api/v2")
eqtl_catalogue_page_size <- as.integer(get_arg("--eqtl-catalogue-page-size", "1000"))
eqtl_catalogue_max_pages <- as.integer(get_arg("--eqtl-catalogue-max-pages", "20"))
eqtl_catalogue_sleep_s <- as.numeric(get_arg("--eqtl-catalogue-sleep-s", "0.25"))

pval_threshold <- as.numeric(get_arg("--pval-threshold", "1e-5"))
clump_r2 <- as.numeric(get_arg("--clump-r2", "0.001"))
clump_kb <- as.integer(get_arg("--clump-kb", "10000"))
clump_pop <- get_arg("--clump-pop", "EUR")
max_instruments_per_gene <- as.integer(get_arg("--max-instruments-per-gene", "200"))

use_proxies <- as.integer(get_arg("--proxies", "0"))
proxy_r2 <- as.numeric(get_arg("--proxy-r2", "0.8"))

maf_ambiguous <- as.numeric(get_arg("--maf-ambiguous", "0.42"))
f_stat_min <- as.numeric(get_arg("--f-stat-min", "10"))

out_dir_tables <- get_arg("--out-tables", "results/tables")
out_dir_figs <- get_arg("--out-figs", "figures/raw_plots")
out_results_name <- get_arg("--results-name", "T4.1_mr_gene_level_results.csv")
out_snp_name <- get_arg("--snp-name", "T4.1_mr_harmonised_snps.tsv")
out_mapping_name <- get_arg("--mapping-name", "T4.1_eqtl_mapping_used.tsv")
out_iv_coverage_name <- get_arg("--iv-coverage-name", "T4.1_iv_coverage_summary.tsv")
out_config_name <- get_arg("--config-name", "T4.1_mr_config.json")
forest_name <- get_arg("--forest-name", "Fig5A_MR_forest.pdf")
scatter_name <- get_arg("--scatter-name", "Fig5B_MR_scatter.pdf")

scatter_top_n <- as.integer(get_arg("--scatter-top-n", "9"))
seed <- as.integer(get_arg("--seed", "42"))
dry_run <- flag_present("--dry-run")
plot_only <- flag_present("--plot-only")
skip_outcome_trait_check <- flag_present("--skip-outcome-trait-check")
expected_outcome_trait_regex <- get_arg("--expected-outcome-trait-regex", "(?i)migraine")
debug_outcome_trait <- flag_present("--debug-outcome-trait")

if (!file.exists(input_genes)) stop(sprintf("Gene list not found: %s", input_genes))
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figs, recursive = TRUE, showWarnings = FALSE)

looks_like_jwt <- function(token) {
  token <- trimws(as.character(token))
  grepl("^[A-Za-z0-9_-]+\\.[A-Za-z0-9_-]+\\.[A-Za-z0-9_-]+$", token)
}

if (!plot_only && !nzchar(opengwas_jwt)) {
  stop(
    paste(
      "Missing OpenGWAS JWT token.",
      "",
      "How to run:",
      "  1) Obtain an OpenGWAS JWT (see https://api.opengwas.io/).",
      "  2) Export it in your shell:  export OPENGWAS_JWT='<token>'",
      "  3) Re-run this script.",
      "",
      "Note: do NOT commit tokens to git; do NOT paste tokens into Implementation Logs.",
      sep = "\n"
    )
  )
}

if (!plot_only && (opengwas_jwt == "PASTE_TOKEN_HERE" || grepl("PASTE_TOKEN_HERE", opengwas_jwt, fixed = TRUE))) {
  stop("OPENGWAS_JWT is still set to the placeholder value 'PASTE_TOKEN_HERE'. Update your `.env.local` and re-run.")
}

if (!plot_only && !looks_like_jwt(opengwas_jwt)) {
  stop(
    paste(
      "OPENGWAS_JWT does not look like a JWT (expected three base64url segments separated by '.').",
      "Fix: update `.env.local` (or your environment) to the exact JWT from https://api.opengwas.io/ and re-run.",
      sep = "\n"
    )
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(jsonlite)
  library(httr2)
})

set.seed(seed)

pub_theme <- function() {
  theme_classic(base_size = 11) +
    theme(
      text = element_text(family = "sans"),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, color = "grey30"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8, color = "grey15"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 8),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "grey96", linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
}

if (plot_only) {
  res_path <- file.path(out_dir_tables, out_results_name)
  snp_path <- file.path(out_dir_tables, out_snp_name)
  if (!file.exists(res_path)) stop(sprintf("--plot-only requires MR results CSV at: %s", res_path))
  if (!file.exists(snp_path)) stop(sprintf("--plot-only requires harmonised SNP TSV at: %s", snp_path))

  mr_df <- suppressMessages(readr::read_csv(res_path, show_col_types = FALSE))
  snp_df <- suppressMessages(readr::read_tsv(snp_path, show_col_types = FALSE))

  primary <- mr_df |> filter(method %in% c("IVW (ratio FE)", "Wald ratio")) |>
    group_by(gene_symbol) |>
    slice(1) |>
    ungroup() |>
    arrange(pval)

  if (nrow(primary) > 0) {
    forest_df <- primary |>
      mutate(
        gene_symbol = factor(gene_symbol, levels = rev(unique(gene_symbol))),
        log_or = beta,
        log_or_l = beta - 1.96 * se,
        log_or_u = beta + 1.96 * se,
        qc_label = dplyr::case_when(
          is.finite(q_pval) & is.finite(egger_intercept_pval) ~ sprintf("Qp=%.2g; Egger p=%.2g", q_pval, egger_intercept_pval),
          is.finite(q_pval) ~ sprintf("Qp=%.2g", q_pval),
          is.finite(egger_intercept_pval) ~ sprintf("Egger p=%.2g", egger_intercept_pval),
          TRUE ~ NA_character_
        )
      )
    span <- max(forest_df$log_or_u, na.rm = TRUE) - min(forest_df$log_or_l, na.rm = TRUE)
    span <- ifelse(is.finite(span) && span > 0, span, 0.02)
    qc_x <- max(forest_df$log_or_u, na.rm = TRUE) + 0.10 * span

    plot_title <- if (nzchar(outcome_trait_arg)) outcome_trait_arg else outcome_id
    p_forest <- ggplot(forest_df, aes(x = log_or, y = gene_symbol)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
      geom_errorbar(aes(xmin = log_or_l, xmax = log_or_u), orientation = "y", width = 0.18, linewidth = 0.6) +
      geom_point(size = 2.2, color = "#1f77b4") +
      geom_text(
        aes(x = qc_x, label = qc_label),
        hjust = 0,
        size = 2.6,
        color = "grey30",
        na.rm = TRUE
      ) +
      labs(
        x = "MR effect (log OR per 1-SD increase in expression)",
        y = NULL,
        title = sprintf("Two-sample MR (Outcome: %s)", plot_title),
        subtitle = sprintf("Primary estimate per gene (IVW/Wald); outcome ID: %s", outcome_id)
      ) +
      pub_theme() +
      expand_limits(x = qc_x + 0.05 * span)

    ggsave(
      filename = file.path(out_dir_figs, forest_name),
      plot = p_forest,
      width = 7.8,
      height = max(3.8, 0.27 * nrow(forest_df) + 1.1),
      units = "in"
    )
  }

  top_genes <- primary$gene_symbol |> as.character() |> head(scatter_top_n)
  scatter_df <- snp_df |> filter(gene_symbol %in% top_genes)
  if (nrow(scatter_df) > 0) {
    scatter_df <- scatter_df |>
      mutate(gene_symbol = factor(gene_symbol, levels = top_genes))

    slopes <- primary |>
      filter(gene_symbol %in% top_genes) |>
      select(gene_symbol, beta) |>
      mutate(gene_symbol = factor(gene_symbol, levels = top_genes))

    plot_title <- if (nzchar(outcome_trait_arg)) outcome_trait_arg else outcome_id
    p_scatter <- ggplot(scatter_df, aes(x = bx, y = by)) +
      geom_hline(yintercept = 0, color = "grey90") +
      geom_vline(xintercept = 0, color = "grey90") +
      geom_point(alpha = 0.75, size = 1.6, color = "grey20") +
      geom_abline(
        data = slopes,
        aes(intercept = 0, slope = beta),
        color = "#2C7FB8",
        linewidth = 0.7
      ) +
      facet_wrap(~gene_symbol, scales = "free") +
      labs(
        x = "SNP effect on expression (beta)",
        y = "SNP effect on migraine (beta)",
        title = sprintf("MR scatter (Outcome: %s)", plot_title),
        subtitle = sprintf("Line = IVW/Wald primary estimate; outcome ID: %s", outcome_id)
      ) +
      pub_theme()

    ggsave(
      filename = file.path(out_dir_figs, scatter_name),
      plot = p_scatter,
      width = 8.8,
      height = 6.8,
      units = "in"
    )
  }

  quit(status = 0)
}

if (!(eqtl_mode %in% c("prefix", "mapping", "gwasinfo", "eqtl_catalogue"))) {
  stop("--eqtl-mode must be one of: prefix, mapping, gwasinfo, eqtl_catalogue")
}
if (!(eqtl_match %in% c("symbol", "ensembl", "both"))) {
  stop("--eqtl-match must be one of: symbol, ensembl, both")
}
if (!is.finite(eqtl_max_hits_per_gene) || eqtl_max_hits_per_gene < 1) {
  stop("--eqtl-max-hits-per-gene must be >= 1")
}
if (!is.finite(eqtl_catalogue_page_size) || eqtl_catalogue_page_size < 1) {
  stop("--eqtl-catalogue-page-size must be >= 1")
}
if (eqtl_catalogue_page_size > 1000) {
  warning(sprintf("--eqtl-catalogue-page-size=%d exceeds API limit; clamping to 1000.", eqtl_catalogue_page_size))
  eqtl_catalogue_page_size <- 1000L
}
if (!is.finite(eqtl_catalogue_max_pages) || eqtl_catalogue_max_pages < 1) {
  stop("--eqtl-catalogue-max-pages must be >= 1")
}
if (!is.finite(eqtl_catalogue_sleep_s) || eqtl_catalogue_sleep_s < 0) {
  stop("--eqtl-catalogue-sleep-s must be >= 0")
}

utc_now_iso <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}

token_fingerprint <- function(token) {
  token <- as.character(token)
  if (!nzchar(token)) return(NA_character_)
  n <- nchar(token)
  if (n <= 8) return(paste0("len=", n))
  paste0(substr(token, 1, 4), "...", substr(token, n - 3, n), " (len=", n, ")")
}

opengwas_base <- "https://api.opengwas.io/api"
opengwas_trial_base <- "https://api.opengwas.io/api"

opengwas_request <- function(endpoint, method = c("POST", "GET"), query = list()) {
  method <- match.arg(method)
  timeout_s <- if (endpoint == "/gwasinfo" && method == "GET") 1200 else 240

  req <- httr2::request(paste0(opengwas_base, endpoint)) |>
    httr2::req_headers(Authorization = paste("Bearer", opengwas_jwt)) |>
    httr2::req_timeout(timeout_s) |>
    httr2::req_retry(max_tries = 3)

  if (length(query) > 0) {
    req <- do.call(httr2::req_url_query, c(list(req), query, list(.multi = "explode")))
  }

  resp <- if (method == "POST") httr2::req_perform(httr2::req_method(req, "POST")) else httr2::req_perform(req)

  if (httr2::resp_status(resp) >= 400) {
    body <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
    stop(sprintf("OpenGWAS API error %s on %s: %s", httr2::resp_status(resp), endpoint, body))
  }

  txt <- httr2::resp_body_string(resp)
  if (!nzchar(txt)) return(NULL)
  jsonlite::fromJSON(txt, simplifyVector = TRUE)
}

normalise_gwasinfo <- function(x) {
  if (is.null(x)) return(tibble::tibble())

  # Case 1: already a data.frame / tibble
  if (is.data.frame(x)) {
    if (!("id" %in% colnames(x))) {
      return(tibble::as_tibble(x))
    }
    return(tibble::as_tibble(x))
  }

  # Case 2: list wrapper with 'data'
  if (is.list(x) && "data" %in% names(x) && is.data.frame(x$data)) {
    return(tibble::as_tibble(x$data))
  }

  # Case 3: OpenGWAS /gwasinfo returns a named list: { "<id>": { ... }, ... }
  if (is.list(x) && !is.null(names(x)) && length(names(x)) > 0) {
    ids <- names(x)
    rows <- lapply(ids, function(id) {
      item <- x[[id]]
      if (is.null(item) || !is.list(item)) return(NULL)
      df <- as.data.frame(item, stringsAsFactors = FALSE)
      df$id <- id
      df
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(tibble::tibble())
    return(dplyr::bind_rows(rows) |> tibble::as_tibble())
  }

  tibble::tibble()
}

match_token_in_text <- function(text, token) {
  if (is.na(text) || !nzchar(text) || is.na(token) || !nzchar(token)) return(FALSE)
  pat <- paste0("(^|[^A-Za-z0-9])", token, "([^A-Za-z0-9]|$)")
  grepl(pat, text, ignore.case = FALSE, perl = TRUE)
}

score_eqtl_candidate <- function(trait, source_pattern = "", tissue_regex = "") {
  s <- 0
  if (!is.na(trait) && nzchar(trait)) {
    if (grepl("cis", trait, ignore.case = TRUE)) s <- s + 3
    if (grepl("trans", trait, ignore.case = TRUE)) s <- s - 1
    if (nzchar(source_pattern) && grepl(source_pattern, trait, ignore.case = TRUE)) s <- s + 5
    if (nzchar(tissue_regex) && grepl(tissue_regex, trait, perl = TRUE)) s <- s + 10
  }
  s
}

discover_eqtl_candidates <- function(
    genes,
    symbol_to_ensembl_path,
    source_pattern = "",
    tissue_regex = "",
    match_mode = "ensembl",
    max_hits_per_gene = 10,
    gwasinfo_cache_path = "") {
  if (!file.exists(symbol_to_ensembl_path)) {
    stop(sprintf("symbol-to-ensembl table not found: %s", symbol_to_ensembl_path))
  }
  gene_tbl <- readr::read_tsv(
    symbol_to_ensembl_path,
    col_names = c("ensembl_id", "gene_symbol"),
    col_types = cols(),
    progress = FALSE
  ) |>
    transmute(gene_symbol = as.character(gene_symbol), ensembl_id = as.character(ensembl_id)) |>
    filter(!is.na(gene_symbol), gene_symbol != "", !is.na(ensembl_id), ensembl_id != "")

  genes <- unique(as.character(genes))
  gene_tbl <- gene_tbl |> filter(gene_symbol %in% genes) |> distinct(gene_symbol, .keep_all = TRUE)

  gwasinfo <- NULL
  if (nzchar(gwasinfo_cache_path) && file.exists(gwasinfo_cache_path)) {
    txt <- readr::read_file(gwasinfo_cache_path)
    gwasinfo <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  } else {
    gwasinfo <- opengwas_request("/gwasinfo", method = "GET")
  }
  gwasinfo <- normalise_gwasinfo(gwasinfo)
  if (nrow(gwasinfo) == 0) {
    stop("OpenGWAS /gwasinfo did not return a table. Token may not have permissions.")
  }

  trait_col <- select_text_col(gwasinfo, c("trait", "phenotype", "name"))
  id_col <- select_text_col(gwasinfo, c("id", "gwas_id"))
  if (is.null(trait_col) || is.null(id_col)) {
    stop(sprintf("OpenGWAS /gwasinfo is missing expected columns. Columns: %s", paste(colnames(gwasinfo), collapse = ", ")))
  }

  gwas_tbl <- gwasinfo |>
    transmute(
      exposure_id = as.character(.data[[id_col]]),
      exposure_trait = as.character(.data[[trait_col]]),
      author = if ("author" %in% colnames(gwasinfo)) as.character(author) else NA_character_,
      consortium = if ("consortium" %in% colnames(gwasinfo)) as.character(consortium) else NA_character_,
      note = if ("note" %in% colnames(gwasinfo)) as.character(note) else NA_character_
    ) |>
    filter(!is.na(exposure_id), exposure_id != "", !is.na(exposure_trait), exposure_trait != "") |>
    mutate(
      meta_text = paste(
        exposure_trait,
        author,
        consortium,
        note,
        sep = " | "
      )
    ) |>
    filter(grepl("^eqtl-", exposure_id, ignore.case = TRUE))

  if (nzchar(source_pattern)) {
    gwas_tbl <- gwas_tbl |> filter(grepl(source_pattern, meta_text, ignore.case = TRUE))
  }
  if (nzchar(tissue_regex)) {
    gwas_tbl <- gwas_tbl |> filter(grepl(tissue_regex, meta_text, perl = TRUE))
  }

  candidates <- lapply(seq_len(nrow(gene_tbl)), function(i) {
    gsym <- gene_tbl$gene_symbol[[i]]
    ens <- gene_tbl$ensembl_id[[i]]
    if (is.na(gsym) || !nzchar(gsym) || is.na(ens) || !nzchar(ens)) return(tibble::tibble())

    m_symbol <- match_mode %in% c("symbol", "both")
    m_ens <- match_mode %in% c("ensembl", "both")

    hits <- gwas_tbl |>
      mutate(
        match_symbol = if (m_symbol) vapply(meta_text, function(t) match_token_in_text(t, gsym), logical(1)) else FALSE,
        match_ensembl = if (m_ens) grepl(ens, exposure_id, fixed = TRUE) | grepl(ens, meta_text, fixed = TRUE) else FALSE
      ) |>
      filter(match_symbol | match_ensembl) |>
      mutate(
        gene_symbol = gsym,
        ensembl_id = ens,
        matched_on = case_when(
          match_symbol & match_ensembl ~ "both",
          match_ensembl ~ "ensembl",
          match_symbol ~ "symbol",
          TRUE ~ "unknown"
        ),
        score = vapply(meta_text, function(t) score_eqtl_candidate(t, source_pattern, tissue_regex), numeric(1))
      ) |>
      arrange(desc(score), exposure_id) |>
      head(max_hits_per_gene) |>
      select(gene_symbol, ensembl_id, exposure_id, exposure_trait, author, consortium, matched_on, score)

    hits
  }) |> bind_rows()

  if (nrow(candidates) == 0) {
    candidates <- tibble::tibble(
      gene_symbol = character(),
      ensembl_id = character(),
      exposure_id = character(),
      exposure_trait = character(),
      matched_on = character(),
      score = numeric()
    )
  }

  candidates
}

read_eqtl_catalogue_metadata <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(tibble::tibble())
  df <- readr::read_tsv(path, col_types = cols(), progress = FALSE)
  tibble::as_tibble(df)
}

build_eqtl_catalogue_mapping <- function(genes, symbol_to_ensembl_path, dataset_id, metadata_path) {
  if (!nzchar(dataset_id)) stop("--eqtl-catalogue-dataset-id must be provided when --eqtl-mode eqtl_catalogue.")
  if (!file.exists(symbol_to_ensembl_path)) {
    stop(sprintf("--eqtl-mode eqtl_catalogue requires --symbol-to-ensembl to exist: %s", symbol_to_ensembl_path))
  }

  md <- read_eqtl_catalogue_metadata(metadata_path)
  md_row <- md |> filter(dataset_id == !!dataset_id) |> slice_head(n = 1)
  tissue_label <- if (nrow(md_row) > 0 && "tissue_label" %in% colnames(md_row)) as.character(md_row$tissue_label[[1]]) else NA_character_
  study_label <- if (nrow(md_row) > 0 && "study_label" %in% colnames(md_row)) as.character(md_row$study_label[[1]]) else NA_character_
  quant_method <- if (nrow(md_row) > 0 && "quant_method" %in% colnames(md_row)) as.character(md_row$quant_method[[1]]) else NA_character_
  sample_size <- if (nrow(md_row) > 0 && "sample_size" %in% colnames(md_row)) as.character(md_row$sample_size[[1]]) else NA_character_

  if (!is.na(quant_method) && nzchar(quant_method) && quant_method != "ge") {
    warning(sprintf("eQTL Catalogue dataset %s quant_method is '%s' (expected 'ge'). Proceeding anyway.", dataset_id, quant_method))
  }

  map_df <- readr::read_tsv(symbol_to_ensembl_path, col_names = c("ensembl_id", "gene_symbol"), col_types = cols(), progress = FALSE)
  map_df <- map_df |>
    mutate(
      ensembl_id = as.character(ensembl_id),
      gene_symbol = as.character(gene_symbol)
    ) |>
    filter(!is.na(ensembl_id), ensembl_id != "", !is.na(gene_symbol), gene_symbol != "")

  picked <- map_df |>
    filter(gene_symbol %in% genes) |>
    group_by(gene_symbol) |>
    slice(1) |>
    ungroup() |>
    mutate(
      exposure_id = paste0("eqtlcat-", dataset_id, "-", ensembl_id),
      exposure_trait = paste0(
        "eQTL Catalogue ",
        ifelse(is.na(study_label) || !nzchar(study_label), "", paste0(study_label, " ")),
        dataset_id,
        ifelse(is.na(tissue_label) || !nzchar(tissue_label), "", paste0(" (", tissue_label, ")")),
        ifelse(is.na(sample_size) || !nzchar(sample_size), "", paste0(" N=", sample_size))
      ),
      eqtl_catalogue_dataset_id = dataset_id
    ) |>
    select(gene_symbol, ensembl_id, exposure_id, exposure_trait, eqtl_catalogue_dataset_id)

  missing <- setdiff(genes, unique(picked$gene_symbol))
  if (length(missing) > 0) {
    picked <- bind_rows(
      picked,
      tibble::tibble(
        gene_symbol = missing,
        ensembl_id = NA_character_,
        exposure_id = NA_character_,
        exposure_trait = NA_character_,
        eqtl_catalogue_dataset_id = dataset_id
      )
    )
  }

  picked |> arrange(gene_symbol)
}

eqtl_catalogue_request <- function(path, query = list()) {
  base <- sub("/+$", "", eqtl_catalogue_api_base)
  url <- paste0(base, "/", sub("^/+", "", path))
  req <- httr2::request(url) |>
    httr2::req_timeout(60) |>
    httr2::req_retry(max_tries = 2)
  if (length(query) > 0) {
    req <- do.call(httr2::req_url_query, c(list(req), query, list(.multi = "explode")))
  }
  resp <- httr2::req_perform(req)
  if (httr2::resp_status(resp) >= 400) {
    body <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
    stop(sprintf("eQTL Catalogue API error %s on %s: %s", httr2::resp_status(resp), url, body))
  }
  txt <- httr2::resp_body_string(resp)
  if (!nzchar(txt)) return(NULL)
  jsonlite::fromJSON(txt, simplifyVector = TRUE)
}

fetch_eqtl_catalogue_hits <- function(dataset_id, gene_id, pval_threshold) {
  if (!nzchar(dataset_id) || !nzchar(gene_id)) return(tibble::tibble())
  nlog10p <- -log10(pval_threshold)
  if (!is.finite(nlog10p)) return(tibble::tibble())

  rows <- list()
  start <- 0
  for (page in seq_len(eqtl_catalogue_max_pages)) {
    if (eqtl_catalogue_sleep_s > 0) Sys.sleep(eqtl_catalogue_sleep_s)
    res <- NULL
    last_err <- NULL
    for (attempt in seq_len(3)) {
      res <- tryCatch(
        eqtl_catalogue_request(
          sprintf("datasets/%s/associations", dataset_id),
          query = list(
            gene_id = gene_id,
            nlog10p = nlog10p,
            start = start,
            size = eqtl_catalogue_page_size
          )
        ),
        error = function(e) {
          last_err <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(res)) break
      Sys.sleep(min(5, 2^(attempt - 1)))
    }

    if (is.null(res)) {
      warning(sprintf("eQTL Catalogue query failed (dataset=%s gene=%s page=%d): %s", dataset_id, gene_id, page, last_err))
      break
    }
    if (is.list(res) && !is.data.frame(res) && "message" %in% names(res)) break

    df <- as.data.frame(res)
    if (nrow(df) == 0) break
    rows[[length(rows) + 1]] <- df

    if (nrow(df) < eqtl_catalogue_page_size) break
    start <- start + eqtl_catalogue_page_size
  }

  if (length(rows) == 0) return(tibble::tibble())

  df <- dplyr::bind_rows(rows) |> tibble::as_tibble()
  df <- df |>
    mutate(
      variant_type = if ("type" %in% colnames(df)) as.character(type) else NA_character_,
      rsid = if ("rsid" %in% colnames(df)) as.character(rsid) else NA_character_,
      beta = if ("beta" %in% colnames(df)) as.numeric(beta) else NA_real_,
      se = if ("se" %in% colnames(df)) as.numeric(se) else NA_real_,
      pval = if ("pvalue" %in% colnames(df)) as.numeric(pvalue) else NA_real_,
      effect_allele = if ("alt" %in% colnames(df)) as.character(alt) else NA_character_,
      other_allele = if ("ref" %in% colnames(df)) as.character(ref) else NA_character_,
      eaf = if ("maf" %in% colnames(df)) as.numeric(maf) else NA_real_
    ) |>
    filter(is.na(variant_type) | variant_type == "SNP") |>
    select(rsid, beta, se, pval, effect_allele, other_allele, eaf) |>
    filter(!is.na(rsid), rsid != "", !is.na(beta), !is.na(se), !is.na(pval), pval <= pval_threshold)

  df
}

opengwas_trial_report_url <- function(gwas_id) {
  url <- paste0(opengwas_trial_base, "/gwasinfo/files/trial")
  payload <- list(id = gwas_id, type = "report")

  req <- httr2::request(url) |>
    httr2::req_headers(`Content-Type` = "application/json") |>
    httr2::req_method("POST") |>
    httr2::req_body_json(payload, auto_unbox = TRUE) |>
    httr2::req_timeout(120) |>
    httr2::req_retry(max_tries = 3)

  resp <- httr2::req_perform(req)
  if (httr2::resp_status(resp) >= 400) {
    body <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
    stop(sprintf("OpenGWAS trial report URL error %s for %s: %s", httr2::resp_status(resp), gwas_id, body))
  }
  txt <- httr2::resp_body_string(resp)
  urls <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  if (length(urls) == 0) {
    stop(sprintf("OpenGWAS trial endpoint returned no report URLs for %s.", gwas_id))
  }
  as.character(urls[[1]])
}

opengwas_trial_trait <- function(gwas_id) {
  url <- opengwas_trial_report_url(gwas_id)
  if (!nzchar(url) || is.na(url)) return(NA_character_)

  req <- httr2::request(url) |>
    httr2::req_headers(Range = "bytes=0-5000000", `Accept-Encoding` = "identity") |>
    httr2::req_timeout(180) |>
    httr2::req_retry(max_tries = 3)
  resp <- httr2::req_perform(req)
  if (httr2::resp_status(resp) >= 400) {
    body <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
    stop(sprintf("OpenGWAS report fetch error %s for %s: %s", httr2::resp_status(resp), gwas_id, body))
  }
  html <- httr2::resp_body_string(resp)
  if (!nzchar(html)) return(NA_character_)

  if (debug_outcome_trait) {
    idx <- regexpr("Trait", html, ignore.case = TRUE)
    idx2 <- regexpr("Trait: <strong>", html, ignore.case = TRUE, fixed = FALSE)
    fixed_has <- grepl("Trait: <strong>", html, fixed = TRUE)
    cat(
      paste0(
        "[debug] outcome_report_url=", url, "\n",
        "[debug] report_status=", httr2::resp_status(resp), "\n",
        "[debug] html_nchar=", nchar(html), "\n",
        "[debug] first_Trait_index=", idx[[1]], "\n"
        ,
        "[debug] has_literal_Trait_strong=", fixed_has, "\n",
        "[debug] first_Trait_strong_index=", idx2[[1]], "\n"
      ),
      file = stderr()
    )
  }

  m <- regexec("Trait:[[:space:]]*<strong>([^<]+)</strong>", html, ignore.case = TRUE)
  mm <- regmatches(html, m)
  if (length(mm) == 0 || length(mm[[1]]) < 2) {
    stop(sprintf("Could not locate 'Trait: <strong>...'</strong> in OpenGWAS report HTML for %s.", gwas_id))
  }
  trimws(mm[[1]][[2]])
}

opengwas_gwasinfo_trait <- function(gwas_id) {
  x <- opengwas_request("/gwasinfo", method = "POST", query = list(id = gwas_id))
  df <- normalise_gwasinfo(x)
  if (nrow(df) == 0) return(NA_character_)
  if (!("trait" %in% colnames(df))) return(NA_character_)
  trait <- df$trait[[1]]
  if (is.null(trait)) return(NA_character_)
  trimws(as.character(trait))
}

read_gene_list <- function(path) {
  df <- readr::read_csv(path, col_types = cols(), progress = FALSE) |> as.data.frame()
  if (!("gene_symbol" %in% colnames(df))) stop("Gene list must contain column: gene_symbol")
  genes <- df$gene_symbol |> as.character() |> trimws()
  genes <- genes[genes != ""]
  unique(genes)
}

build_eqtl_mapping_from_symbol_table <- function(genes, symbol_to_ensembl_path, eqtl_id_prefix) {
  if (!file.exists(symbol_to_ensembl_path)) {
    stop(sprintf("symbol-to-ensembl file not found: %s", symbol_to_ensembl_path))
  }
  map_df <- readr::read_tsv(symbol_to_ensembl_path, col_names = c("ensembl_id", "gene_symbol"), col_types = cols(), progress = FALSE)
  map_df <- map_df |>
    mutate(
      ensembl_id = as.character(ensembl_id),
      gene_symbol = as.character(gene_symbol)
    ) |>
    filter(!is.na(ensembl_id), ensembl_id != "", !is.na(gene_symbol), gene_symbol != "")

  picked <- map_df |>
    filter(gene_symbol %in% genes) |>
    group_by(gene_symbol) |>
    slice(1) |>
    ungroup() |>
    mutate(
      exposure_id = paste0(eqtl_id_prefix, ensembl_id),
      exposure_trait = ensembl_id
    ) |>
    select(gene_symbol, ensembl_id, exposure_id, exposure_trait)

  missing <- setdiff(genes, unique(picked$gene_symbol))
  if (length(missing) > 0) {
    picked <- bind_rows(
      picked,
      tibble::tibble(
        gene_symbol = missing,
        ensembl_id = NA_character_,
        exposure_id = NA_character_,
        exposure_trait = NA_character_
      )
    )
  }

  picked |> arrange(gene_symbol)
}

select_text_col <- function(df, candidates) {
  for (c in candidates) {
    if (c %in% colnames(df)) return(c)
  }
  NULL
}

select_numeric_col <- function(df, candidates) {
  for (c in candidates) {
    if (c %in% colnames(df)) return(c)
  }
  NULL
}

extract_variant_col <- function(df) {
  select_text_col(df, c("rsid", "variant", "SNP"))
}

extract_beta_col <- function(df) {
  select_numeric_col(df, c("beta", "b", "beta_hat"))
}

extract_se_col <- function(df) {
  select_numeric_col(df, c("se", "sebeta", "standard_error"))
}

extract_pval_col <- function(df) {
  select_numeric_col(df, c("p", "pval", "p_value"))
}

extract_ea_col <- function(df) {
  select_text_col(df, c("ea", "effect_allele", "a1"))
}

extract_oa_col <- function(df) {
  select_text_col(df, c("nea", "other_allele", "a2"))
}

extract_eaf_col <- function(df) {
  select_numeric_col(df, c("eaf", "eaf1", "af", "effect_allele_freq"))
}

standardise_assoc <- function(df, prefix) {
  if (is.null(df)) return(tibble::tibble())
  if (is.list(df) && length(df) == 0) return(tibble::tibble())
  df <- as.data.frame(df)
  if (nrow(df) == 0 || ncol(df) == 0) return(tibble::tibble())
  vcol <- extract_variant_col(df)
  bcol <- extract_beta_col(df)
  scol <- extract_se_col(df)
  pcol <- extract_pval_col(df)
  eacol <- extract_ea_col(df)
  oacol <- extract_oa_col(df)
  eafcol <- extract_eaf_col(df)

  required <- c(vcol, bcol, scol, pcol, eacol, oacol)
  if (any(is.null(required))) {
    stop(
      sprintf(
        "Could not standardise %s columns; got columns: %s",
        prefix,
        paste(colnames(df), collapse = ", ")
      )
    )
  }

  out <- tibble::tibble(
    rsid = as.character(df[[vcol]]),
    beta = as.numeric(df[[bcol]]),
    se = as.numeric(df[[scol]]),
    pval = as.numeric(df[[pcol]]),
    effect_allele = as.character(df[[eacol]]),
    other_allele = as.character(df[[oacol]]),
    eaf = if (!is.null(eafcol)) as.numeric(df[[eafcol]]) else NA_real_
  )

  out <- out |>
    mutate(
      rsid = trimws(rsid),
      effect_allele = toupper(trimws(effect_allele)),
      other_allele = toupper(trimws(other_allele))
    ) |>
    filter(!is.na(rsid), rsid != "", !is.na(beta), !is.na(se))

  colnames(out) <- c(
    "rsid",
    paste0("beta_", prefix),
    paste0("se_", prefix),
    paste0("pval_", prefix),
    paste0("effect_allele_", prefix),
    paste0("other_allele_", prefix),
    paste0("eaf_", prefix)
  )

  out
}

allele_complement <- function(a) {
  a <- toupper(a)
  map <- c(A = "T", T = "A", C = "G", G = "C")
  out <- map[a]
  out[is.na(out)] <- NA_character_
  out
}

is_palindromic <- function(a1, a2) {
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  (a1 == "A" & a2 == "T") |
    (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") |
    (a1 == "G" & a2 == "C")
}

harmonise_exposure_outcome <- function(merged, maf_ambiguous = 0.42) {
  merged <- merged |> mutate(
    palindromic = is_palindromic(effect_allele_exposure, other_allele_exposure),
    maf_exposure = pmin(eaf_exposure, 1 - eaf_exposure),
    maf_outcome = pmin(eaf_outcome, 1 - eaf_outcome),
    maf_min = pmin(maf_exposure, maf_outcome, na.rm = TRUE),
    ambiguous_palindrome = palindromic & (is.na(maf_min) | maf_min > maf_ambiguous),
    keep = TRUE
  )

  merged$keep[merged$ambiguous_palindrome] <- FALSE

  flip_needed <- rep(FALSE, nrow(merged))

  a1e <- merged$effect_allele_exposure
  a2e <- merged$other_allele_exposure
  a1o <- merged$effect_allele_outcome
  a2o <- merged$other_allele_outcome

  direct <- a1e == a1o & a2e == a2o
  swapped <- a1e == a2o & a2e == a1o

  a1o_c <- allele_complement(a1o)
  a2o_c <- allele_complement(a2o)
  direct_c <- a1e == a1o_c & a2e == a2o_c
  swapped_c <- a1e == a2o_c & a2e == a1o_c

  merged$keep[!(direct | swapped | direct_c | swapped_c)] <- FALSE
  flip_needed[swapped | swapped_c] <- TRUE

  merged <- merged |> mutate(
    beta_outcome_aligned = ifelse(flip_needed, -beta_outcome, beta_outcome),
    effect_allele_outcome_aligned = ifelse(flip_needed, other_allele_outcome, effect_allele_outcome),
    other_allele_outcome_aligned = ifelse(flip_needed, effect_allele_outcome, other_allele_outcome)
  )

  merged |> filter(keep)
}

mr_wald_ratio <- function(bx, sebx, by, seby) {
  theta <- by / bx
  se_theta <- sqrt((seby^2) / (bx^2) + (by^2) * (sebx^2) / (bx^4))
  z <- theta / se_theta
  p <- 2 * stats::pnorm(-abs(z))
  list(method = "Wald ratio", beta = theta, se = se_theta, pval = p, nsnp = 1, q = NA_real_, q_pval = NA_real_, egger_intercept = NA_real_, egger_intercept_pval = NA_real_)
}

mr_ivw_ratio_fe <- function(bx, sebx, by, seby) {
  theta_i <- by / bx
  var_i <- (seby^2) / (bx^2) + (by^2) * (sebx^2) / (bx^4)
  w <- 1 / var_i
  theta <- sum(w * theta_i) / sum(w)
  se_theta <- sqrt(1 / sum(w))
  z <- theta / se_theta
  p <- 2 * stats::pnorm(-abs(z))
  q <- sum(w * (theta_i - theta)^2)
  q_p <- stats::pchisq(q, df = length(theta_i) - 1, lower.tail = FALSE)
  list(method = "IVW (ratio FE)", beta = theta, se = se_theta, pval = p, nsnp = length(theta_i), q = q, q_pval = q_p, egger_intercept = NA_real_, egger_intercept_pval = NA_real_)
}

weighted_median <- function(theta, w) {
  ord <- order(theta)
  theta <- theta[ord]
  w <- w[ord]
  w <- w / sum(w)
  cw <- cumsum(w)
  idx <- which(cw >= 0.5)[1]
  theta[[idx]]
}

mr_weighted_median <- function(bx, sebx, by, seby) {
  theta_i <- by / bx
  var_i <- (seby^2) / (bx^2) + (by^2) * (sebx^2) / (bx^4)
  w <- 1 / var_i
  theta <- weighted_median(theta_i, w)
  # Approximate SE using weighted MAD / normal approximation (lightweight; report as exploratory).
  mad <- stats::mad(theta_i, constant = 1)
  se_theta <- 1.2533 * mad / sqrt(length(theta_i))
  z <- theta / se_theta
  p <- 2 * stats::pnorm(-abs(z))
  list(method = "Weighted median", beta = theta, se = se_theta, pval = p, nsnp = length(theta_i), q = NA_real_, q_pval = NA_real_, egger_intercept = NA_real_, egger_intercept_pval = NA_real_)
}

mr_egger <- function(bx, by, seby) {
  w <- 1 / (seby^2)
  fit <- stats::lm(by ~ bx, weights = w)
  sm <- summary(fit)
  slope <- sm$coefficients["bx", "Estimate"]
  slope_se <- sm$coefficients["bx", "Std. Error"]
  slope_p <- sm$coefficients["bx", "Pr(>|t|)"]
  int <- sm$coefficients["(Intercept)", "Estimate"]
  int_se <- sm$coefficients["(Intercept)", "Std. Error"]
  int_p <- sm$coefficients["(Intercept)", "Pr(>|t|)"]
  list(method = "MR-Egger", beta = slope, se = slope_se, pval = slope_p, nsnp = length(bx), q = NA_real_, q_pval = NA_real_, egger_intercept = int, egger_intercept_pval = int_p)
}

format_mr_row <- function(gene, res, outcome_id) {
  or <- exp(res$beta)
  ci_l <- exp(res$beta - 1.96 * res$se)
  ci_u <- exp(res$beta + 1.96 * res$se)
  tibble::tibble(
    gene_symbol = gene,
    outcome_id = outcome_id,
    method = res$method,
    nsnp = res$nsnp,
    beta = res$beta,
    se = res$se,
    pval = res$pval,
    or = or,
    ci_lower = ci_l,
    ci_upper = ci_u,
    q = res$q,
    q_pval = res$q_pval,
    egger_intercept = res$egger_intercept,
    egger_intercept_pval = res$egger_intercept_pval
  )
}

genes <- read_gene_list(input_genes)

if (dry_run) {
  trial_error <- NULL
  gwasinfo_error <- NULL
  trait_trial <- tryCatch(opengwas_trial_trait(outcome_id), error = function(e) {
    trial_error <<- conditionMessage(e)
    NA_character_
  })
  trait_gwasinfo <- tryCatch(opengwas_gwasinfo_trait(outcome_id), error = function(e) {
    gwasinfo_error <<- conditionMessage(e)
    NA_character_
  })
  trait_verified <- if (!is.na(trait_trial) && nzchar(trait_trial)) trait_trial else trait_gwasinfo
  trait_source <- if (!is.na(trait_trial) && nzchar(trait_trial)) "trial_report" else "gwasinfo"
  cat(paste0(
    "Dry run OK.\n",
    "- genes: ", length(genes), "\n",
    "- outcome_id: ", outcome_id, "\n",
    "- outcome_trait(trial): ", ifelse(is.na(trait_trial), "(unavailable)", trait_trial), "\n",
    if (!is.null(trial_error)) paste0("- outcome_trait_error: ", trial_error, "\n") else "",
    "- outcome_trait(gwasinfo): ", ifelse(is.na(trait_gwasinfo), "(unavailable)", trait_gwasinfo), "\n",
    if (!is.null(gwasinfo_error)) paste0("- outcome_gwasinfo_error: ", gwasinfo_error, "\n") else "",
    "- outcome_trait(verified): ", ifelse(is.na(trait_verified) || !nzchar(trait_verified), "(unavailable)", trait_verified), "\n",
    "- outcome_trait_source: ", trait_source, "\n",
    "- opengwas_token: ", token_fingerprint(opengwas_jwt), "\n"
  ))
  quit(status = 0)
}

trial_error <- NULL
gwasinfo_error <- NULL
outcome_trait_trial <- tryCatch(opengwas_trial_trait(outcome_id), error = function(e) {
  trial_error <<- conditionMessage(e)
  NA_character_
})
outcome_trait_gwasinfo <- tryCatch(opengwas_gwasinfo_trait(outcome_id), error = function(e) {
  gwasinfo_error <<- conditionMessage(e)
  NA_character_
})
outcome_trait <- if (!is.na(outcome_trait_trial) && nzchar(outcome_trait_trial)) outcome_trait_trial else outcome_trait_gwasinfo
outcome_trait_source <- if (!is.na(outcome_trait_trial) && nzchar(outcome_trait_trial)) "trial_report" else "gwasinfo"
if (!skip_outcome_trait_check) {
  if (is.na(outcome_trait) || !nzchar(outcome_trait)) {
    stop(
      paste(
        "Could not verify outcome trait via OpenGWAS metadata.",
        if (!is.null(trial_error)) paste0("Trial report error: ", trial_error) else NULL,
        if (!is.null(gwasinfo_error)) paste0("GWAS info error: ", gwasinfo_error) else NULL,
        "Re-run with --skip-outcome-trait-check to bypass (not recommended), or provide a different --outcome-id.",
        sep = "\n"
      )
    )
  }
  if (!grepl(expected_outcome_trait_regex, outcome_trait, perl = TRUE)) {
    stop(
      paste(
        sprintf("Outcome ID trait mismatch for %s.", outcome_id),
        sprintf("Public report trait: %s", outcome_trait),
        sprintf("Expected regex: %s", expected_outcome_trait_regex),
        "Fix: choose a migraine outcome GWAS ID and pass --outcome-id <id> (or adjust --expected-outcome-trait-regex).",
        sep = "\n"
      )
    )
  }
}

if (eqtl_mode == "mapping") {
  if (!nzchar(eqtl_mapping_path) || !file.exists(eqtl_mapping_path)) {
    stop("--eqtl-mode mapping requires an existing --eqtl-mapping <file>.")
  }
  map_in <- if (grepl("\\.tsv$|\\.txt$", eqtl_mapping_path, ignore.case = TRUE)) {
    readr::read_tsv(eqtl_mapping_path, col_types = cols(), progress = FALSE)
  } else {
    readr::read_csv(eqtl_mapping_path, col_types = cols(), progress = FALSE)
  }
  if (!("gene_symbol" %in% colnames(map_in)) || !("exposure_id" %in% colnames(map_in))) {
    stop("eqtl-mapping file must contain columns: gene_symbol, exposure_id (optional: exposure_trait).")
  }
  mapping <- map_in |>
    transmute(
      gene_symbol = as.character(gene_symbol),
      ensembl_id = if ("ensembl_id" %in% colnames(map_in)) as.character(ensembl_id) else NA_character_,
      exposure_id = as.character(exposure_id),
      exposure_trait = if ("exposure_trait" %in% colnames(map_in)) as.character(exposure_trait) else NA_character_
    ) |>
    filter(!is.na(gene_symbol), gene_symbol != "") |>
    distinct(gene_symbol, .keep_all = TRUE)
} else if (eqtl_mode == "prefix") {
  if (!nzchar(eqtl_id_prefix)) {
    stop("--eqtl-mode prefix requires a non-empty --eqtl-id-prefix.")
  }
  if (!file.exists(symbol_to_ensembl_path)) {
    stop(sprintf("--eqtl-mode prefix requires --symbol-to-ensembl to exist: %s", symbol_to_ensembl_path))
  }
  # Preferred path: use known OpenGWAS eQTL dataset naming convention: eqtl-a-<ENSG...>
  mapping <- build_eqtl_mapping_from_symbol_table(genes, symbol_to_ensembl_path, eqtl_id_prefix)
} else if (eqtl_mode == "gwasinfo") {
  candidates <- discover_eqtl_candidates(
    genes = genes,
    symbol_to_ensembl_path = symbol_to_ensembl_path,
    source_pattern = eqtl_source_pattern,
    tissue_regex = eqtl_tissue_regex,
    match_mode = eqtl_match,
    max_hits_per_gene = eqtl_max_hits_per_gene,
    gwasinfo_cache_path = gwasinfo_cache
  )

  if (nzchar(eqtl_candidates_out)) {
    if (grepl("\\.tsv$|\\.txt$", eqtl_candidates_out, ignore.case = TRUE)) {
      readr::write_tsv(candidates, eqtl_candidates_out)
    } else {
      readr::write_csv(candidates, eqtl_candidates_out)
    }
  }

  mapping <- candidates |>
    group_by(gene_symbol) |>
    arrange(desc(score), exposure_id) |>
    slice_head(n = 1) |>
    ungroup() |>
    select(gene_symbol, ensembl_id, exposure_id, exposure_trait)

  missing <- setdiff(genes, mapping$gene_symbol)
  if (length(missing) > 0) {
    mapping <- bind_rows(
      mapping,
      tibble::tibble(
        gene_symbol = as.character(missing),
        ensembl_id = NA_character_,
        exposure_id = NA_character_,
        exposure_trait = NA_character_
      )
    )
  }
  mapping <- mapping |>
    arrange(gene_symbol)
} else if (eqtl_mode == "eqtl_catalogue") {
  mapping <- build_eqtl_catalogue_mapping(genes, symbol_to_ensembl_path, eqtl_catalogue_dataset_id, eqtl_catalogue_metadata_path)
} else {
  stop("Unhandled --eqtl-mode.")
}

if ("eqtl_catalogue_dataset_id" %in% colnames(mapping)) {
  readr::write_tsv(mapping |> select(gene_symbol, ensembl_id, exposure_id, exposure_trait, eqtl_catalogue_dataset_id), file.path(out_dir_tables, out_mapping_name))
} else {
  readr::write_tsv(mapping, file.path(out_dir_tables, out_mapping_name))
}

if (eqtl_discover_only) {
  message("eqtl-discover-only requested; stopping after writing mapping.")
  quit(status = 0)
}

genes_with_exposure <- mapping |> filter(!is.na(exposure_id) & exposure_id != "") |> pull(gene_symbol) |> as.character()
if (length(genes_with_exposure) == 0) {
  stop(
    paste(
      "No matching eQTL exposure IDs found for the provided genes.",
      sprintf("Try providing --eqtl-mapping <file> (gene_symbol -> exposure_id)."),
      sprintf("Or provide --symbol-to-ensembl <file> and --eqtl-id-prefix (current prefix: %s).", eqtl_id_prefix),
      sprintf("Fallback legacy mode: adjust --eqtl-source-pattern (current: %s) if using trait-string matching.", eqtl_source_pattern),
      sep = "\n"
    )
  )
}

get_tophits <- function(exposure_id) {
  res <- tryCatch(
    opengwas_request(
      "/tophits",
      method = "POST",
      query = list(
        id = exposure_id,
        pval = pval_threshold,
        preclumped = 0,
        clump = 0,
        pop = clump_pop
      )
    ),
    error = function(e) {
      warning(sprintf("OpenGWAS tophits failed for exposure_id=%s: %s", exposure_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) return(tibble::tibble())
  standardise_assoc(res, "exposure")
}

get_eqtl_catalogue_tophits <- function(ensembl_id) {
  hits <- fetch_eqtl_catalogue_hits(eqtl_catalogue_dataset_id, ensembl_id, pval_threshold)
  if (nrow(hits) == 0) return(tibble::tibble())
  standardise_assoc(hits, "exposure")
}

clump_instruments <- function(rsid, pval) {
  if (length(rsid) == 0) return(character())
  rsid <- as.character(rsid)
  pval <- as.numeric(pval)
  keep <- !is.na(rsid) & nzchar(rsid) & !is.na(pval) & is.finite(pval)
  rsid <- rsid[keep]
  pval <- pval[keep]
  if (length(rsid) == 0) return(character())

  # Some OpenGWAS datasets may emit p=0 in /tophits due to underflow; avoid sending literal 0.
  pval[pval <= 0] <- 1e-300

  # Use form body instead of URL query to avoid HTTP 414 "URI Too Long" for multi-variant clumping.
  res <- tryCatch(
    {
      req <- httr2::request(paste0(opengwas_base, "/ld/clump")) |>
        httr2::req_headers(Authorization = paste("Bearer", opengwas_jwt)) |>
        httr2::req_timeout(240) |>
        httr2::req_retry(max_tries = 3) |>
        httr2::req_method("POST") |>
        httr2::req_body_form(
          rsid = rsid,
          pval = pval,
          pthresh = pval_threshold,
          r2 = clump_r2,
          kb = clump_kb,
          pop = clump_pop
          ,
          .multi = "explode"
        )
      resp <- httr2::req_perform(req)
      if (httr2::resp_status(resp) >= 400) {
        body <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
        stop(sprintf("OpenGWAS API error %s on /ld/clump: %s", httr2::resp_status(resp), body))
      }
      txt <- httr2::resp_body_string(resp)
      if (!nzchar(txt)) return(NULL)
      jsonlite::fromJSON(txt, simplifyVector = TRUE)
    },
    error = function(e) {
      warning(sprintf("OpenGWAS clump failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) return(character())
  if (is.character(res)) return(unique(res))
  if (is.list(res) && length(res) == 0) return(character())
  df <- as.data.frame(res)
  vcol <- extract_variant_col(df)
  if (is.null(vcol)) {
    # Some responses may return a single unnamed vector coerced into a column; handle gracefully.
    return(character())
  }
  unique(as.character(df[[vcol]]))
}

get_outcome_assoc <- function(rsids) {
  if (length(rsids) == 0) return(tibble::tibble())
  res <- tryCatch(
    opengwas_request(
      "/associations",
      method = "POST",
      query = list(
        variant = rsids,
        id = outcome_id,
        proxies = use_proxies,
        r2 = proxy_r2,
        align_alleles = 1,
        palindromes = 1
      )
    ),
    error = function(e) {
      warning(sprintf("OpenGWAS associations failed for outcome_id=%s: %s", outcome_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) return(tibble::tibble())
  standardise_assoc(res, "outcome")
}

all_snp_rows <- list()
mr_rows <- list()
iv_coverage_rows <- list()

for (i in seq_len(nrow(mapping))) {
  gene <- mapping$gene_symbol[[i]]
  exposure_id <- mapping$exposure_id[[i]]
  exposure_trait <- mapping$exposure_trait[[i]]
  ensembl_id <- if ("ensembl_id" %in% colnames(mapping)) mapping$ensembl_id[[i]] else NA_character_

  if (is.na(exposure_id) || !nzchar(exposure_id)) {
    iv_coverage_rows[[length(iv_coverage_rows) + 1]] <- tibble::tibble(
      gene_symbol = gene,
      ensembl_id = ensembl_id,
      exposure_id = NA_character_,
      exposure_trait = NA_character_,
      outcome_id = outcome_id,
      n_tophits = 0L,
      n_after_f = 0L,
      n_after_cap = 0L,
      n_after_clump = 0L,
      n_after_outcome_join = 0L,
      n_after_harmonise = 0L,
      nsnp_used = 0L,
      status = "missing_exposure_id"
    )
    next
  }

  hits0 <- if (eqtl_mode == "eqtl_catalogue") get_eqtl_catalogue_tophits(ensembl_id) else get_tophits(exposure_id)
  n_tophits <- as.integer(nrow(hits0))
  if (n_tophits == 0L) {
    iv_coverage_rows[[length(iv_coverage_rows) + 1]] <- tibble::tibble(
      gene_symbol = gene,
      ensembl_id = ensembl_id,
      exposure_id = exposure_id,
      exposure_trait = exposure_trait,
      outcome_id = outcome_id,
      n_tophits = n_tophits,
      n_after_f = 0L,
      n_after_cap = 0L,
      n_after_clump = 0L,
      n_after_outcome_join = 0L,
      n_after_harmonise = 0L,
      nsnp_used = 0L,
      status = "no_tophits"
    )
    next
  }

  f_stat <- (hits0$beta_exposure^2) / (hits0$se_exposure^2)
  hits1 <- hits0 |>
    mutate(f_stat = f_stat) |>
    filter(is.finite(f_stat), f_stat >= f_stat_min) |>
    arrange(pval_exposure) |>
    slice_head(n = max_instruments_per_gene)
  n_after_f <- as.integer(sum(is.finite(f_stat) & f_stat >= f_stat_min))
  n_after_cap <- as.integer(nrow(hits1))
  if (n_after_cap == 0L) {
    iv_coverage_rows[[length(iv_coverage_rows) + 1]] <- tibble::tibble(
      gene_symbol = gene,
      ensembl_id = ensembl_id,
      exposure_id = exposure_id,
      exposure_trait = exposure_trait,
      outcome_id = outcome_id,
      n_tophits = n_tophits,
      n_after_f = n_after_f,
      n_after_cap = n_after_cap,
      n_after_clump = 0L,
      n_after_outcome_join = 0L,
      n_after_harmonise = 0L,
      nsnp_used = 0L,
      status = "no_strong_instruments_after_f"
    )
    next
  }

  clumped <- clump_instruments(hits1$rsid, hits1$pval_exposure)
  hits2 <- hits1 |> filter(rsid %in% clumped)
  n_after_clump <- as.integer(nrow(hits2))
  if (n_after_clump == 0L) {
    iv_coverage_rows[[length(iv_coverage_rows) + 1]] <- tibble::tibble(
      gene_symbol = gene,
      ensembl_id = ensembl_id,
      exposure_id = exposure_id,
      exposure_trait = exposure_trait,
      outcome_id = outcome_id,
      n_tophits = n_tophits,
      n_after_f = n_after_f,
      n_after_cap = n_after_cap,
      n_after_clump = n_after_clump,
      n_after_outcome_join = 0L,
      n_after_harmonise = 0L,
      nsnp_used = 0L,
      status = "all_dropped_after_clump"
    )
    next
  }

  out <- get_outcome_assoc(hits2$rsid)
  if (nrow(out) == 0) {
    iv_coverage_rows[[length(iv_coverage_rows) + 1]] <- tibble::tibble(
      gene_symbol = gene,
      ensembl_id = ensembl_id,
      exposure_id = exposure_id,
      exposure_trait = exposure_trait,
      outcome_id = outcome_id,
      n_tophits = n_tophits,
      n_after_f = n_after_f,
      n_after_cap = n_after_cap,
      n_after_clump = n_after_clump,
      n_after_outcome_join = 0L,
      n_after_harmonise = 0L,
      nsnp_used = 0L,
      status = "no_outcome_associations"
    )
    next
  }

  merged0 <- hits2 |>
    inner_join(out, by = "rsid") |>
    mutate(
      gene_symbol = gene,
      exposure_id = exposure_id,
      exposure_trait = exposure_trait
    )

  n_after_outcome_join <- as.integer(nrow(merged0))
  merged <- harmonise_exposure_outcome(merged0, maf_ambiguous = maf_ambiguous)
  n_after_harmonise <- as.integer(nrow(merged))
  if (n_after_harmonise == 0L) {
    iv_coverage_rows[[length(iv_coverage_rows) + 1]] <- tibble::tibble(
      gene_symbol = gene,
      ensembl_id = ensembl_id,
      exposure_id = exposure_id,
      exposure_trait = exposure_trait,
      outcome_id = outcome_id,
      n_tophits = n_tophits,
      n_after_f = n_after_f,
      n_after_cap = n_after_cap,
      n_after_clump = n_after_clump,
      n_after_outcome_join = n_after_outcome_join,
      n_after_harmonise = n_after_harmonise,
      nsnp_used = 0L,
      status = "all_dropped_after_harmonise"
    )
    next
  }

  iv_coverage_rows[[length(iv_coverage_rows) + 1]] <- tibble::tibble(
    gene_symbol = gene,
    ensembl_id = ensembl_id,
    exposure_id = exposure_id,
    exposure_trait = exposure_trait,
    outcome_id = outcome_id,
    n_tophits = n_tophits,
    n_after_f = n_after_f,
    n_after_cap = n_after_cap,
    n_after_clump = n_after_clump,
    n_after_outcome_join = n_after_outcome_join,
    n_after_harmonise = n_after_harmonise,
    nsnp_used = as.integer(n_after_harmonise),
    status = ifelse(n_after_harmonise >= 3, "mr_multi_snp_possible", ifelse(n_after_harmonise == 2, "mr_two_snp_limited", "mr_single_snp_only"))
  )

  merged <- merged |>
    mutate(
      bx = beta_exposure,
      sebx = se_exposure,
      by = beta_outcome_aligned,
      seby = se_outcome,
      f_stat = (bx^2) / (sebx^2)
    )

  all_snp_rows[[length(all_snp_rows) + 1]] <- merged

  bx <- merged$bx
  sebx <- merged$sebx
  by <- merged$by
  seby <- merged$seby

  methods <- list()
  if (length(bx) == 1) {
    methods[["primary"]] <- mr_wald_ratio(bx, sebx, by, seby)
  } else {
    methods[["primary"]] <- mr_ivw_ratio_fe(bx, sebx, by, seby)
  }

  if (length(bx) >= 3) {
    methods[["egger"]] <- mr_egger(bx, by, seby)
    methods[["wmed"]] <- mr_weighted_median(bx, sebx, by, seby)
  }

  for (k in names(methods)) {
    res <- methods[[k]]
    row <- format_mr_row(gene, res, outcome_id) |>
      mutate(exposure_id = exposure_id, exposure_trait = exposure_trait)
    mr_rows[[length(mr_rows) + 1]] <- row
  }
}

if (length(all_snp_rows) == 0 || length(mr_rows) == 0) {
  stop(
    paste(
      "No MR results produced.",
      "Common causes:",
      "- No instruments passed the p-value / F-stat threshold",
      "- Outcome associations missing for selected SNPs",
      "- eQTL source pattern did not match your accessible datasets",
      sep = "\n"
    )
  )
}

snp_df <- bind_rows(all_snp_rows) |>
  select(
    gene_symbol, exposure_id, exposure_trait, rsid,
    bx, sebx, by, seby,
    effect_allele_exposure, other_allele_exposure,
    effect_allele_outcome_aligned, other_allele_outcome_aligned,
    eaf_exposure, eaf_outcome, f_stat,
    pval_exposure, pval_outcome
  )

mr_df <- bind_rows(mr_rows) |>
  mutate(
    pval = as.numeric(pval),
    method = as.character(method)
  ) |>
  mutate(pval_fdr_bh = p.adjust(pval, method = "BH")) |>
  relocate(pval_fdr_bh, .after = pval) |>
  arrange(method, pval)

readr::write_tsv(snp_df, file.path(out_dir_tables, out_snp_name))
readr::write_csv(mr_df, file.path(out_dir_tables, out_results_name))

iv_cov_df <- bind_rows(iv_coverage_rows) |>
  arrange(desc(nsnp_used), gene_symbol)
readr::write_tsv(iv_cov_df, file.path(out_dir_tables, out_iv_coverage_name))

primary <- mr_df |> filter(method %in% c("IVW (ratio FE)", "Wald ratio")) |>
  group_by(gene_symbol) |>
  slice(1) |>
  ungroup() |>
  arrange(pval)

if (nrow(primary) > 0) {
  forest_df <- primary |>
    mutate(
      gene_symbol = factor(gene_symbol, levels = rev(unique(gene_symbol))),
      log_or = beta,
      log_or_l = beta - 1.96 * se,
      log_or_u = beta + 1.96 * se,
      qc_label = dplyr::case_when(
        is.finite(q_pval) & is.finite(egger_intercept_pval) ~ sprintf("Qp=%.2g; Egger p=%.2g", q_pval, egger_intercept_pval),
        is.finite(q_pval) ~ sprintf("Qp=%.2g", q_pval),
        is.finite(egger_intercept_pval) ~ sprintf("Egger p=%.2g", egger_intercept_pval),
        TRUE ~ NA_character_
      )
    )
  span <- max(forest_df$log_or_u, na.rm = TRUE) - min(forest_df$log_or_l, na.rm = TRUE)
  span <- ifelse(is.finite(span) && span > 0, span, 0.02)
  qc_x <- max(forest_df$log_or_u, na.rm = TRUE) + 0.10 * span
  plot_title <- if (nzchar(outcome_trait)) outcome_trait else outcome_id

  p_forest <- ggplot(forest_df, aes(x = log_or, y = gene_symbol)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = log_or_l, xmax = log_or_u), orientation = "y", width = 0.18, linewidth = 0.6) +
    geom_point(size = 2.2, color = "#1f77b4") +
    geom_text(
      aes(x = qc_x, label = qc_label),
      hjust = 0,
      size = 2.6,
      color = "grey30",
      na.rm = TRUE
    ) +
    labs(
      x = "MR effect (log OR per 1-SD increase in expression)",
      y = NULL,
      title = sprintf("Two-sample MR (Outcome: %s)", plot_title),
      subtitle = sprintf("Primary estimate per gene; outcome ID: %s (p<=%g, r2<=%g, kb=%d, F>=%g)", outcome_id, pval_threshold, clump_r2, clump_kb, f_stat_min)
    ) +
    pub_theme() +
    expand_limits(x = qc_x + 0.05 * span)

  ggsave(
    filename = file.path(out_dir_figs, forest_name),
    plot = p_forest,
    width = 7.5,
    height = max(3.5, 0.25 * nrow(forest_df) + 1),
    units = "in"
  )
}

top_genes <- primary$gene_symbol |> as.character() |> head(scatter_top_n)
scatter_df <- snp_df |> filter(gene_symbol %in% top_genes)
if (nrow(scatter_df) > 0) {
  scatter_df <- scatter_df |>
    mutate(gene_symbol = factor(gene_symbol, levels = top_genes))

  slopes <- primary |>
    filter(gene_symbol %in% top_genes) |>
    select(gene_symbol, beta) |>
    mutate(gene_symbol = factor(gene_symbol, levels = top_genes))

  plot_title <- if (nzchar(outcome_trait)) outcome_trait else outcome_id
  p_scatter <- ggplot(scatter_df, aes(x = bx, y = by)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_point(alpha = 0.75, size = 1.6, color = "grey20") +
    geom_abline(
      data = slopes,
      aes(intercept = 0, slope = beta),
      color = "#2C7FB8",
      linewidth = 0.7
    ) +
    facet_wrap(~gene_symbol, scales = "free") +
    labs(
      x = "SNP effect on expression (beta)",
      y = "SNP effect on migraine (beta)",
      title = sprintf("MR scatter (Outcome: %s)", plot_title),
      subtitle = sprintf("Line = IVW/Wald primary estimate; outcome ID: %s", outcome_id)
    ) +
    pub_theme()

  ggsave(
    filename = file.path(out_dir_figs, scatter_name),
    plot = p_scatter,
    width = 8.5,
    height = 6.5,
    units = "in"
  )
}

config <- list(
  created_at_utc = utc_now_iso(),
  script = "scripts/04_mr_causality/01_two_sample_mr_opengwas.R",
  inputs = list(
    genes = input_genes,
    outcome_id = outcome_id,
    outcome_trait_trial = outcome_trait_trial,
    outcome_trait_gwasinfo = outcome_trait_gwasinfo,
    outcome_trait_verified = outcome_trait,
    outcome_trait_source = outcome_trait_source
  ),
  exposure = list(
    eqtl_mode = eqtl_mode,
    eqtl_catalogue = if (eqtl_mode == "eqtl_catalogue") list(
      dataset_id = eqtl_catalogue_dataset_id,
      metadata_path = eqtl_catalogue_metadata_path,
      api_base = eqtl_catalogue_api_base,
      page_size = eqtl_catalogue_page_size,
      max_pages = eqtl_catalogue_max_pages,
      sleep_s = eqtl_catalogue_sleep_s
    ) else NULL
  ),
  opengwas = list(
    base = opengwas_base,
    token_fingerprint = token_fingerprint(opengwas_jwt),
    eqtl_source_pattern = eqtl_source_pattern,
    eqtl_mapping = if (nzchar(eqtl_mapping_path)) eqtl_mapping_path else NULL,
    symbol_to_ensembl = if (file.exists(symbol_to_ensembl_path)) symbol_to_ensembl_path else NULL,
    eqtl_id_prefix = eqtl_id_prefix
  ),
  parameters = list(
    pval_threshold = pval_threshold,
    clump_r2 = clump_r2,
    clump_kb = clump_kb,
    clump_pop = clump_pop,
    max_instruments_per_gene = max_instruments_per_gene,
    proxies = use_proxies,
    proxy_r2 = proxy_r2,
    maf_ambiguous = maf_ambiguous,
    f_stat_min = f_stat_min,
    scatter_top_n = scatter_top_n,
    seed = seed
  ),
  outputs = list(
    mapping = file.path(out_dir_tables, out_mapping_name),
    iv_coverage = file.path(out_dir_tables, out_iv_coverage_name),
    snps = file.path(out_dir_tables, out_snp_name),
    results = file.path(out_dir_tables, out_results_name),
    forest = file.path(out_dir_figs, forest_name),
    scatter = file.path(out_dir_figs, scatter_name)
  )
)

jsonlite::write_json(config, file.path(out_dir_tables, out_config_name), pretty = TRUE, auto_unbox = TRUE)
