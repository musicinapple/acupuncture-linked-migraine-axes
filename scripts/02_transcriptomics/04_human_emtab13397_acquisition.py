#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import logging
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


LOGGER = logging.getLogger("emtab13397_acquisition")


ASSAY_RE = re.compile(r"^(?P<subject>\d{3})(?P<session>S\d)(?P<rep>\d)(?P<suffix>[A-Z])$")
ENSG_RE = re.compile(r"\bENSG\d+\b")


@dataclass(frozen=True)
class AssayParts:
    subject_id: str
    session: str
    replicate: str
    suffix: str


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_assay_name(assay_name: str) -> AssayParts:
    m = ASSAY_RE.match(assay_name.strip())
    if not m:
        raise ValueError(f"Unexpected assay name format: {assay_name}")
    return AssayParts(
        subject_id=m.group("subject"),
        session=m.group("session"),
        replicate=m.group("rep"),
        suffix=m.group("suffix"),
    )


def extract_ensembl_id(gene_id: str) -> str:
    m = ENSG_RE.search(gene_id)
    return m.group(0) if m else ""


def require_mygene() -> None:
    try:
        import mygene  # noqa: F401
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("mygene is required. Install with: pip3 install mygene") from exc


def load_mapping_cache(cache_path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not os.path.exists(cache_path):
        return mapping
    with open(cache_path, "r", encoding="utf-8") as file_handle:
        for line in file_handle:
            line = line.rstrip("\n")
            if not line:
                continue
            key, value = line.split("\t", 1)
            mapping[key] = value
    return mapping


def save_mapping_cache(cache_path: str, mapping: Dict[str, str]) -> None:
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as file_handle:
        for key in sorted(mapping.keys()):
            file_handle.write(f"{key}\t{mapping[key]}\n")


def map_ensembl_to_symbol(ensembl_ids: List[str], cache_path: str) -> Dict[str, str]:
    require_mygene()
    import mygene

    mapping = load_mapping_cache(cache_path)
    missing = sorted({eid for eid in ensembl_ids if eid and eid not in mapping})
    if not missing:
        return mapping

    LOGGER.info("Mapping %s Ensembl IDs via mygene (species=human)", len(missing))
    mg = mygene.MyGeneInfo()
    results = mg.querymany(
        missing,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=False,
        returnall=False,
        verbose=False,
    )

    for item in results:
        query = str(item.get("query", "")).strip()
        symbol = str(item.get("symbol", "")).strip()
        if query and symbol:
            mapping[query] = symbol
        elif query and query not in mapping:
            mapping[query] = ""

    save_mapping_cache(cache_path, mapping)
    return mapping


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    x = counts.to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    lib_sizes = x.sum(axis=0)
    lib_sizes[lib_sizes == 0] = 1.0
    cpm = x / lib_sizes * 1e6
    return pd.DataFrame(np.log2(cpm + 1.0), index=counts.index, columns=counts.columns)


def save_pca_figure(expr_logcpm: pd.DataFrame, meta: pd.DataFrame, output_pdf: str) -> None:
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA

    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    meta_key = "sample_id" if "sample_id" in meta.columns else ("index" if "index" in meta.columns else None)
    if not meta_key:
        raise KeyError("Metadata must contain 'sample_id' (or 'index') column for PCA alignment.")
    meta = meta.set_index(meta_key).loc[expr_logcpm.columns, :]

    groups = meta["group"].astype(str).tolist()
    sessions = meta["session"].astype(str).tolist()
    group_colors = {"case": "#D62728", "control": "#1F77B4"}
    session_markers = {"S1": "o", "S2": "^"}

    x = expr_logcpm.T.to_numpy(dtype=float)
    x = x - x.mean(axis=0, keepdims=True)
    pcs = PCA(n_components=2, svd_solver="full", random_state=12345).fit_transform(x)

    fig, ax = plt.subplots(1, 1, figsize=(6.5, 5.0))
    for i in range(pcs.shape[0]):
        ax.scatter(
            pcs[i, 0],
            pcs[i, 1],
            s=34,
            color=group_colors.get(groups[i], "#666666"),
            marker=session_markers.get(sessions[i], "o"),
            alpha=0.85,
            edgecolor="white",
            linewidth=0.5,
        )
    ax.set_title("E-MTAB-13397 PCA (log2-CPM)", fontsize=11)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=group_colors["control"], markersize=7, label="control"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=group_colors["case"], markersize=7, label="case"),
        plt.Line2D([0], [0], marker=session_markers["S1"], color="#333333", linestyle="None", markersize=7, label="S1"),
        plt.Line2D([0], [0], marker=session_markers["S2"], color="#333333", linestyle="None", markersize=7, label="S2"),
    ]
    ax.legend(handles=handles, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(output_pdf, dpi=300)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build human transcriptomics matrix for ArrayExpress E-MTAB-13397.")
    parser.add_argument(
        "--input-dir",
        default="data/raw/arrayexpress/E-MTAB-13397/Files",
        help="Directory containing gene_count_matrix.txt and SDRF/IDF files.",
    )
    parser.add_argument("--counts", default="gene_count_matrix.txt", help="Counts matrix filename (CSV).")
    parser.add_argument("--sdrf", default="E-MTAB-13397.sdrf.txt", help="SDRF filename (tab-delimited).")
    parser.add_argument("--out-matrix", default="data/processed/human_expression_matrix.tsv", help="Output log2-CPM matrix TSV (gene symbols x samples).")
    parser.add_argument("--out-meta", default="data/processed/human_sample_metadata.tsv", help="Output sample metadata TSV.")
    parser.add_argument("--out-mapping-cache", default="data/processed/human_id_to_symbol.tsv", help="Output mapping cache TSV (Ensembl->symbol).")
    parser.add_argument("--fig-pca", default="figures/raw_plots/Fig3H_human_PCA_emtab13397.pdf", help="PCA figure PDF.")
    parser.add_argument("--min-total-n", type=int, default=30, help="Circuit breaker: stop if total N < this.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    args = parser.parse_args()

    configure_logging(args.verbose)
    retrieved_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    counts_path = os.path.join(args.input_dir, args.counts)
    sdrf_path = os.path.join(args.input_dir, args.sdrf)

    if not os.path.exists(counts_path):
        raise FileNotFoundError(counts_path)
    if not os.path.exists(sdrf_path):
        raise FileNotFoundError(sdrf_path)

    LOGGER.info("Loading SDRF: %s", sdrf_path)
    sdrf = pd.read_csv(sdrf_path, sep="\t")
    required_cols = [
        "Characteristics[analysis id]",
        "Assay Name",
        "Characteristics[disease]",
        "Characteristics[sex]",
        "Characteristics[age]",
    ]
    for col in required_cols:
        if col not in sdrf.columns:
            raise RuntimeError(f"Missing required SDRF column: {col}")

    sdrf = sdrf.copy()
    sdrf["sample_id"] = sdrf["Characteristics[analysis id]"].astype(str)
    sdrf["assay_name"] = sdrf["Assay Name"].astype(str)
    disease = sdrf["Characteristics[disease]"].astype(str).str.strip().str.lower()
    sdrf["group"] = disease.map({"normal": "control", "interictal migraine": "case"}).fillna("unknown")
    sdrf["disease"] = disease
    sdrf["sex"] = sdrf["Characteristics[sex]"].astype(str).str.strip().str.lower()
    sdrf["age_years"] = pd.to_numeric(sdrf["Characteristics[age]"].astype(str).str.replace(",", ".", regex=False), errors="coerce")

    parts: List[AssayParts] = [parse_assay_name(a) for a in sdrf["assay_name"].astype(str).tolist()]
    sdrf["subject_id"] = [p.subject_id for p in parts]
    sdrf["session"] = [p.session for p in parts]
    sdrf["replicate"] = [p.replicate for p in parts]

    keep_meta_cols = [
        "sample_id",
        "assay_name",
        "subject_id",
        "session",
        "replicate",
        "group",
        "disease",
        "sex",
        "age_years",
        "Characteristics[allergy history]",
        "Characteristics[smoking status]",
        "Comment[ENA_RUN]",
        "Comment[ENA_SAMPLE]",
        "Comment[ENA_EXPERIMENT]",
        "Characteristics[manufacturer id]",
        "Characteristics[organism part]",
    ]
    existing_meta_cols = [c for c in keep_meta_cols if c in sdrf.columns]
    meta = sdrf[existing_meta_cols].copy()
    meta["batch"] = "E-MTAB-13397"
    meta["retrieved_at"] = retrieved_at
    meta = meta.drop_duplicates(subset=["sample_id"], keep="first")

    LOGGER.info("Loading counts matrix: %s", counts_path)
    counts_raw = pd.read_csv(counts_path)
    if "gene_id" not in counts_raw.columns:
        raise RuntimeError("Counts matrix missing 'gene_id' column.")
    gene_ids = counts_raw["gene_id"].astype(str)
    counts_raw = counts_raw.drop(columns=["gene_id"])
    counts_raw.index = gene_ids

    counts_raw = counts_raw.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts_raw.index.name = "gene_id"

    # Align to SDRF samples by analysis id.
    sdrf_ids = set(meta["sample_id"].astype(str).tolist())
    original_cols = counts_raw.columns.astype(str).tolist()

    # The count-matrix column IDs sometimes use '-' while SDRF uses '.' (same tokens).
    normalized_cols = [c.replace("-", ".") for c in original_cols]
    rename_map = {orig: norm for orig, norm in zip(original_cols, normalized_cols)}
    counts_raw = counts_raw.rename(columns=rename_map)

    cols_keep = [c for c in counts_raw.columns.astype(str).tolist() if c in sdrf_ids]
    if not cols_keep:
        raise RuntimeError("No overlapping sample IDs between counts matrix columns and SDRF analysis IDs.")
    counts_raw = counts_raw.loc[:, cols_keep]

    # Reorder metadata to match counts columns.
    meta = meta.set_index("sample_id").loc[counts_raw.columns, :].reset_index()
    if "sample_id" not in meta.columns and "index" in meta.columns:
        meta = meta.rename(columns={"index": "sample_id"})

    total_n = counts_raw.shape[1]
    if total_n < args.min_total_n:
        raise RuntimeError(f"Circuit breaker: total sample size N={total_n} < {args.min_total_n}. Stop.")

    # Map Ensembl -> symbol and aggregate by symbol.
    ensembl_ids = [extract_ensembl_id(gid) for gid in counts_raw.index.astype(str).tolist()]
    mapping = map_ensembl_to_symbol([eid for eid in ensembl_ids if eid], cache_path=args.out_mapping_cache)
    symbols = pd.Series([mapping.get(eid, "") if eid else "" for eid in ensembl_ids], index=counts_raw.index, name="gene_symbol")

    counts = counts_raw.copy()
    counts.insert(0, "gene_symbol", symbols.values)
    counts = counts[counts["gene_symbol"].astype(str).str.strip() != ""]
    counts = counts.set_index("gene_symbol")
    counts.index = counts.index.astype(str)
    counts = counts.groupby(counts.index).sum(numeric_only=True)
    counts.index.name = "gene_symbol"

    expr = log_cpm(counts)

    os.makedirs(os.path.dirname(args.out_matrix), exist_ok=True)
    expr.to_csv(args.out_matrix, sep="\t")
    os.makedirs(os.path.dirname(args.out_meta), exist_ok=True)
    meta.to_csv(args.out_meta, sep="\t", index=False)
    save_pca_figure(expr, meta, args.fig_pca)

    n_cases = int((meta["group"] == "case").sum())
    n_controls = int((meta["group"] == "control").sum())
    LOGGER.info("Saved: %s (genes=%s samples=%s)", args.out_matrix, expr.shape[0], expr.shape[1])
    LOGGER.info("Metadata: %s (case=%s control=%s)", args.out_meta, n_cases, n_controls)
    LOGGER.info("Figure: %s", args.fig_pca)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
