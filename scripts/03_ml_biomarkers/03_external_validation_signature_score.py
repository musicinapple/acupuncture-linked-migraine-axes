#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import platform
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve


LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class SignatureResult:
    auc: float
    auc_ci_lower: Optional[float]
    auc_ci_upper: Optional[float]
    n_total: int
    n_cases: int
    n_controls: int
    n_signature_genes: int
    n_genes_used: int
    genes_missing: list[str]


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _run_cmd(cmd: list[str]) -> str:
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
        return out.strip()
    except Exception as exc:  # noqa: BLE001 - keep logging resilient
        return f"<unavailable: {exc}>"


def _read_expression_matrix(path: Path, gene_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={gene_col: "string"})
    if gene_col not in df.columns:
        raise ValueError(f"Expression matrix missing required gene column '{gene_col}': {path}")

    df = df.dropna(subset=[gene_col]).copy()
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df = df[df[gene_col] != ""].copy()

    sample_cols = [c for c in df.columns if c != gene_col]
    if not sample_cols:
        raise ValueError(f"Expression matrix has no sample columns: {path}")

    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = (
        df.groupby(gene_col, as_index=False)[sample_cols]
        .mean(numeric_only=True)
        .sort_values(gene_col)
        .reset_index(drop=True)
    )
    return df


def _read_signature_genes(path: Path, gene_col: str) -> list[str]:
    if path.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t", dtype={gene_col: "string"})
    else:
        df = pd.read_csv(path, dtype={gene_col: "string"})

    if gene_col not in df.columns:
        raise ValueError(f"Signature file missing gene column '{gene_col}': {path}")

    genes = df[gene_col].dropna().astype(str).str.strip()
    genes = genes[genes != ""].unique().tolist()
    return genes


def _read_deg_directions(
    path: Optional[Path],
    gene_col: str,
    logfc_col: str,
) -> dict[str, int]:
    if path is None:
        return {}

    df = pd.read_csv(path)
    if gene_col not in df.columns or logfc_col not in df.columns:
        raise ValueError(
            f"DEG table missing required columns '{gene_col}' and '{logfc_col}': {path}"
        )

    df = df.dropna(subset=[gene_col, logfc_col]).copy()
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df = df[df[gene_col] != ""].copy()
    df[logfc_col] = pd.to_numeric(df[logfc_col], errors="coerce")
    df = df.dropna(subset=[logfc_col]).copy()

    directions: dict[str, int] = {}
    for gene, logfc in zip(df[gene_col], df[logfc_col]):
        directions[gene] = 1 if float(logfc) >= 0 else -1
    return directions


def _signature_scores(
    expression_df: pd.DataFrame,
    signature_genes: Iterable[str],
    directions: dict[str, int],
    gene_col: str,
) -> tuple[pd.Series, list[str]]:
    signature_genes = [g.strip() for g in signature_genes if g and str(g).strip()]
    signature_genes = list(dict.fromkeys(signature_genes))

    expr = expression_df.set_index(gene_col)
    present = [g for g in signature_genes if g in expr.index]
    missing = [g for g in signature_genes if g not in expr.index]

    if not present:
        raise ValueError("No signature genes are present in the validation expression matrix.")

    m = expr.loc[present].to_numpy(dtype=float)
    m = np.where(np.isfinite(m), m, np.nan)

    means = np.nanmean(m, axis=1, keepdims=True)
    stds = np.nanstd(m, axis=1, keepdims=True)
    stds = np.where(stds == 0, np.nan, stds)
    z = (m - means) / stds
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)

    signs = np.array([directions.get(g, 1) for g in present], dtype=float).reshape(-1, 1)
    score = np.mean(z * signs, axis=0)
    return pd.Series(score, index=expr.columns, name="signature_score"), missing


def _bootstrap_auc_ci(
    y_true: np.ndarray,
    y_score: np.ndarray,
    n_bootstrap: int,
    seed: int,
) -> tuple[Optional[float], Optional[float]]:
    if n_bootstrap <= 0:
        return None, None

    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)

    idx_case = np.where(y_true == 1)[0]
    idx_ctrl = np.where(y_true == 0)[0]

    if len(idx_case) == 0 or len(idx_ctrl) == 0:
        return None, None

    aucs: list[float] = []
    for _ in range(n_bootstrap):
        res_case = rng.choice(idx_case, size=len(idx_case), replace=True)
        res_ctrl = rng.choice(idx_ctrl, size=len(idx_ctrl), replace=True)
        idx = np.concatenate([res_case, res_ctrl])
        try:
            aucs.append(float(roc_auc_score(y_true[idx], y_score[idx])))
        except Exception:
            continue

    if not aucs:
        return None, None

    lower = float(np.quantile(aucs, 0.025))
    upper = float(np.quantile(aucs, 0.975))
    return lower, upper


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, sort_keys=True)
        f.write("\n")


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="External validation of a fixed signature via direction-aware z-score scoring.",
    )
    p.add_argument(
        "--expression-matrix",
        required=True,
        type=Path,
        help="TSV with gene_symbol column + GSM sample columns.",
    )
    p.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="TSV metadata with sample_id + group columns.",
    )
    p.add_argument(
        "--signature",
        required=True,
        type=Path,
        help="CSV/TSV listing signature genes (default column: gene_symbol).",
    )
    p.add_argument(
        "--deg-table",
        type=Path,
        default=None,
        help="Optional DEG table to define gene direction signs using logFC.",
    )
    p.add_argument("--expr-gene-col", default="gene_symbol")
    p.add_argument("--signature-gene-col", default="gene_symbol")
    p.add_argument("--deg-gene-col", default="gene_symbol")
    p.add_argument("--deg-logfc-col", default="logFC")
    p.add_argument("--metadata-sample-col", default="sample_id")
    p.add_argument("--metadata-group-col", default="group")
    p.add_argument("--case-label", default="case")
    p.add_argument("--out-scores", required=True, type=Path)
    p.add_argument("--out-summary", required=True, type=Path)
    p.add_argument("--out-config", required=True, type=Path)
    p.add_argument("--out-roc-figure", required=True, type=Path)
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    LOG.info("Reading validation expression matrix: %s", args.expression_matrix)
    expr = _read_expression_matrix(args.expression_matrix, gene_col=args.expr_gene_col)

    LOG.info("Reading validation metadata: %s", args.metadata)
    meta = pd.read_csv(args.metadata, sep="\t", dtype="string")
    if args.metadata_sample_col not in meta.columns or args.metadata_group_col not in meta.columns:
        raise ValueError(
            f"Metadata missing required columns '{args.metadata_sample_col}' and "
            f"'{args.metadata_group_col}': {args.metadata}"
        )
    meta[args.metadata_sample_col] = meta[args.metadata_sample_col].astype(str).str.strip()
    meta[args.metadata_group_col] = meta[args.metadata_group_col].astype(str).str.strip()
    meta = meta.dropna(subset=[args.metadata_sample_col, args.metadata_group_col]).copy()
    meta = meta[meta[args.metadata_sample_col] != ""].copy()

    signature_genes = _read_signature_genes(args.signature, gene_col=args.signature_gene_col)
    LOG.info("Signature genes: %d", len(signature_genes))

    directions = _read_deg_directions(
        args.deg_table, gene_col=args.deg_gene_col, logfc_col=args.deg_logfc_col
    )
    if directions:
        LOG.info("Loaded DEG directions for %d genes from %s", len(directions), args.deg_table)
    else:
        LOG.info("No DEG directions provided; default sign=+1 for all genes.")

    score, missing_genes = _signature_scores(
        expr, signature_genes=signature_genes, directions=directions, gene_col=args.expr_gene_col
    )

    sample_ids = score.index.tolist()
    meta = meta[meta[args.metadata_sample_col].isin(sample_ids)].copy()
    if meta.empty:
        raise ValueError("No overlapping samples between metadata and expression matrix.")

    meta = meta.drop_duplicates(subset=[args.metadata_sample_col]).copy()
    meta = meta.sort_values(args.metadata_sample_col)

    score_df = (
        meta[[args.metadata_sample_col, args.metadata_group_col]]
        .rename(
            columns={
                args.metadata_sample_col: "sample_id",
                args.metadata_group_col: "group",
            }
        )
        .copy()
    )
    score_df["signature_score"] = score.loc[score_df["sample_id"]].to_numpy(dtype=float)
    score_df["y_true"] = (score_df["group"] == args.case_label).astype(int)

    y_true = score_df["y_true"].to_numpy(dtype=int)
    y_score = score_df["signature_score"].to_numpy(dtype=float)

    n_cases = int(np.sum(y_true == 1))
    n_controls = int(np.sum(y_true == 0))
    if n_cases == 0 or n_controls == 0:
        raise ValueError("Need both case and control samples to compute ROC/AUC.")

    auc = float(roc_auc_score(y_true, y_score))
    auc_inverted = float(roc_auc_score(y_true, -y_score))
    ci_low, ci_high = _bootstrap_auc_ci(y_true, y_score, args.bootstrap, seed=args.seed)

    fpr, tpr, _ = roc_curve(y_true, y_score)

    # Plot (matplotlib imported lazily to keep failures obvious if missing)
    import matplotlib.pyplot as plt  # noqa: PLC0415
    
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    args.out_roc_figure.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(5.2, 4.2))
    plt.plot(fpr, tpr, color="#1f77b4", lw=2, label="ROC")
    plt.plot([0, 1], [0, 1], color="gray", lw=1, linestyle="--", label="Chance")
    summary_text = f"AUC={auc:.3f}"
    if ci_low is not None and ci_high is not None and not (math.isnan(ci_low) or math.isnan(ci_high)):
        summary_text += f", 95% CI {ci_low:.3f}-{ci_high:.3f}"
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    ax = plt.gca()
    ax.text(
        0.95,
        0.05,
        summary_text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"),
    )
    plt.legend(loc="lower right", frameon=False)
    plt.tight_layout()
    plt.savefig(args.out_roc_figure, dpi=300)
    plt.close()

    args.out_scores.parent.mkdir(parents=True, exist_ok=True)
    score_df.to_csv(args.out_scores, index=False)

    summary = SignatureResult(
        auc=auc,
        auc_ci_lower=ci_low,
        auc_ci_upper=ci_high,
        n_total=int(len(score_df)),
        n_cases=n_cases,
        n_controls=n_controls,
        n_signature_genes=int(len(signature_genes)),
        n_genes_used=int(len(signature_genes) - len(missing_genes)),
        genes_missing=missing_genes,
    )

    summary_payload = {
        "method": "direction-aware per-gene z-score; sample score = mean(sign * zscore_gene)",
        "note": (
            "AUC is computed using the pre-defined score direction (case-label and sign definition). "
            "If AUC<0.5, the signature separates classes in the opposite direction; "
            "AUC for the inverted score is provided for transparency."
        ),
        "dataset": {
            "expression_matrix": str(args.expression_matrix),
            "metadata": str(args.metadata),
        },
        "signature": {
            "file": str(args.signature),
            "n_genes": summary.n_signature_genes,
        },
        "deg_directions": {
            "deg_table": str(args.deg_table) if args.deg_table else None,
            "deg_gene_col": args.deg_gene_col,
            "deg_logfc_col": args.deg_logfc_col,
        },
        "results": {
            "auc": summary.auc,
            "auc_inverted_score": auc_inverted,
            "auc_ci_95": [summary.auc_ci_lower, summary.auc_ci_upper],
            "n_total": summary.n_total,
            "n_cases": summary.n_cases,
            "n_controls": summary.n_controls,
            "n_genes_used": summary.n_genes_used,
            "genes_used": sorted(set(signature_genes) - set(missing_genes)),
            "genes_missing": summary.genes_missing,
        },
        "created_at": _utc_now_iso(),
        "environment": {
            "os": platform.platform(),
            "python": _run_cmd(["python3", "-V"]),
        },
    }
    _write_json(args.out_summary, summary_payload)

    config_payload = {
        "created_at": _utc_now_iso(),
        "argv": list(map(str, argv if argv is not None else [])),
        "params": {
            "expression_matrix": str(args.expression_matrix),
            "metadata": str(args.metadata),
            "signature": str(args.signature),
            "deg_table": str(args.deg_table) if args.deg_table else None,
            "expr_gene_col": args.expr_gene_col,
            "signature_gene_col": args.signature_gene_col,
            "metadata_sample_col": args.metadata_sample_col,
            "metadata_group_col": args.metadata_group_col,
            "case_label": args.case_label,
            "bootstrap": args.bootstrap,
            "seed": args.seed,
        },
    }
    _write_json(args.out_config, config_payload)

    LOG.info(
        "Done. AUC=%.3f (cases=%d, controls=%d); genes used=%d/%d; missing=%d",
        auc,
        n_cases,
        n_controls,
        summary.n_genes_used,
        summary.n_signature_genes,
        len(missing_genes),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
