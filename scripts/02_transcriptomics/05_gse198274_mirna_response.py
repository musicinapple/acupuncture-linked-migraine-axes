#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import logging
import os
import re
import tempfile
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats


LOGGER = logging.getLogger("gse198274_mirna_response")


SOFT_SAMPLE_RE = re.compile(r"^\^SAMPLE = (GSM\d+)\s*$")
CHAR_RE = re.compile(r"^!Sample_characteristics_ch1 =\s*(?P<key>[^:]+):\s*(?P<value>.*)\s*$")
TITLE_RE = re.compile(r"^!Sample_title =\s*(.*)\s*$")
REL_RE = re.compile(r"^!Sample_relation =\s*(.*)\s*$")


@dataclass(frozen=True)
class SoftSample:
    gsm: str
    title: str
    treatment: str
    tissue: str
    age: str
    biosample_url: str
    srx_url: str


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)


def _benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    out = np.full(n, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * (len(ranked) / (np.arange(1, len(ranked) + 1)))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out_idx = np.where(ok)[0][order]
    out[out_idx] = q
    return out


def parse_soft_samples(soft_gz_path: str) -> List[SoftSample]:
    samples: List[SoftSample] = []
    current: Dict[str, str] = {}
    current_rel: List[str] = []

    def flush() -> None:
        nonlocal current, current_rel
        if not current:
            return
        gsm = current.get("gsm", "").strip()
        if not gsm:
            current = {}
            current_rel = []
            return
        rel_joined = "\n".join(current_rel)
        biosample = ""
        srx = ""
        for line in current_rel:
            if line.startswith("BioSample:"):
                biosample = line.split("BioSample:", 1)[1].strip()
            if line.startswith("SRA:"):
                srx = line.split("SRA:", 1)[1].strip()
        samples.append(
            SoftSample(
                gsm=gsm,
                title=current.get("title", "").strip(),
                treatment=current.get("treatment", "").strip(),
                tissue=current.get("tissue", "").strip(),
                age=current.get("age", "").strip(),
                biosample_url=biosample,
                srx_url=srx,
            )
        )
        current = {}
        current_rel = []

    with gzip.open(soft_gz_path, "rt", encoding="utf-8", errors="replace") as file_handle:
        for line in file_handle:
            line = line.rstrip("\n")
            m = SOFT_SAMPLE_RE.match(line)
            if m:
                flush()
                current["gsm"] = m.group(1)
                continue

            m = TITLE_RE.match(line)
            if m:
                current["title"] = m.group(1)
                continue

            m = CHAR_RE.match(line)
            if m:
                key = m.group("key").strip().lower()
                value = m.group("value").strip()
                if key in {"treatment", "tissue", "age"}:
                    current[key] = value
                continue

            m = REL_RE.match(line)
            if m:
                current_rel.append(m.group(1).strip())
                continue

    flush()
    if not samples:
        raise RuntimeError(f"No samples parsed from SOFT: {soft_gz_path}")
    return samples


def load_expression_table(path_gz: str) -> pd.DataFrame:
    if not os.path.exists(path_gz):
        raise FileNotFoundError(path_gz)
    with gzip.open(path_gz, "rt", encoding="utf-8", errors="replace") as file_handle:
        df = pd.read_csv(file_handle, sep="\t")
    first_col = df.columns[0]
    if first_col != "#ID":
        # Be tolerant to variants
        LOGGER.warning("Unexpected first column name %r; treating as miRNA id.", first_col)
    df = df.rename(columns={first_col: "mirna"})
    df["mirna"] = df["mirna"].astype(str)
    numeric_cols = [c for c in df.columns if c != "mirna"]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=numeric_cols, how="all")
    return df


def classify_columns(columns: Iterable[str]) -> pd.DataFrame:
    rows = []
    for col in columns:
        name = str(col)
        group = "unknown"
        subject_id = ""
        timepoint = ""
        if name.startswith("HCS"):
            group = "control"
            timepoint = "untreated"
            subject_id = name
        elif "before acupuncture" in name.lower():
            group = "case"
            timepoint = "before"
            subject_id = name.split(" before", 1)[0].strip()
        elif "after acupuncture" in name.lower():
            group = "case"
            timepoint = "after"
            subject_id = name.split(" after", 1)[0].strip()
        rows.append({"sample_name": name, "group": group, "timepoint": timepoint, "subject_id": subject_id})
    return pd.DataFrame(rows)


def compute_contrast_unpaired(
    x: np.ndarray,  # features x samples
    group_a_idx: np.ndarray,
    group_b_idx: np.ndarray,
    label_a: str,
    label_b: str,
) -> Tuple[np.ndarray, np.ndarray]:
    a = x[:, group_a_idx]
    b = x[:, group_b_idx]
    mean_a = np.nanmean(a, axis=1)
    mean_b = np.nanmean(b, axis=1)
    log2fc = mean_a - mean_b
    _, pvals = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    return log2fc, pvals


def compute_contrast_paired(
    x: np.ndarray,  # features x samples
    idx_a: np.ndarray,
    idx_b: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    a = x[:, idx_a]
    b = x[:, idx_b]
    diff = b - a
    log2fc = np.nanmean(diff, axis=1)
    _, pvals = stats.ttest_1samp(diff, popmean=0.0, axis=1, nan_policy="omit")
    return log2fc, pvals


def volcano_plot(df: pd.DataFrame, output_pdf: str, title: str, x_label: str) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    # SOTA 2025: Clean, publication-ready style
    sns.set_style("whitegrid")
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 10,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        }
    )
    
    x = df["log2fc"].to_numpy(dtype=float)
    p = df["padj"].to_numpy(dtype=float)
    y = -np.log10(np.clip(p, 1e-300, 1.0))
    
    # Define significance status
    df = df.copy()
    df['neglog10_padj'] = y # Ensure Y column exists in DataFrame
    df['status'] = 'NS'
    df.loc[(df['padj'] < 0.05) & (df['log2fc'] >= 0.5), 'status'] = 'Up'
    df.loc[(df['padj'] < 0.05) & (df['log2fc'] <= -0.5), 'status'] = 'Down'
    
    # SOTA Colors: Navy/Firebrick (Same as Fig 3A)
    colors = {"Up": "#firebrick3", "Down": "#000080", "NS": "#grey85"} 
    # Use hex codes for safety
    hex_colors = {"Up": "#C44E52", "Down": "#4C72B0", "NS": "#D3D3D3"}

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sns.scatterplot(
        data=df, x='log2fc', y='neglog10_padj', hue='status',
        palette=hex_colors, alpha=0.7, s=20, edgecolor=None, ax=ax
    )
    
    ax.axvline(-0.5, color="#444444", lw=1, ls="--", alpha=0.4)
    ax.axvline(0.5, color="#444444", lw=1, ls="--", alpha=0.4)
    ax.axhline(-np.log10(0.05), color="#444444", lw=1, ls="--", alpha=0.4)
    
    ax.set_title(title, loc="center", fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel("-log10(FDR)", fontsize=10)
    
    # Remove legend title and place it better
    ax.legend(title=None, frameon=True, framealpha=0.9, loc='upper right', fontsize=9)
    
    sig_count = (df['status'] != 'NS').sum()
    ax.text(
        0.05,
        0.05,
        f"Significant: {sig_count}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight='bold',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
    )
    
    fig.tight_layout()
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

def key_mirna_panel(
    expr_log2: pd.DataFrame,
    sample_meta: pd.DataFrame,
    key_mirnas: List[str],
    output_pdf: str,
) -> pd.DataFrame:
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    meta = sample_meta.copy()
    if "sample_name" not in meta.columns:
        raise ValueError("sample_meta must include a 'sample_name' column.")
    meta["sample_name"] = meta["sample_name"].astype(str)
    meta["subject_id"] = meta["subject_id"].astype(str)

    rows: List[Dict[str, object]] = []
    for mir in key_mirnas:
        if mir not in expr_log2.index:
            rows.append({"mirna": mir, "present": False})
            continue
        vals = expr_log2.loc[mir]
        ctrl_names = meta.loc[meta["group"] == "control", "sample_name"].tolist()
        bef_names = meta.loc[(meta["group"] == "case") & (meta["timepoint"] == "before"), "sample_name"].tolist()
        aft_names = meta.loc[(meta["group"] == "case") & (meta["timepoint"] == "after"), "sample_name"].tolist()

        ctrl = vals[ctrl_names].astype(float)
        bef = vals[bef_names].astype(float)
        aft = vals[aft_names].astype(float)

        # paired by subject_id (Axxx)
        before_by = (
            meta[(meta["group"] == "case") & (meta["timepoint"] == "before")][["subject_id", "sample_name"]]
            .drop_duplicates()
            .set_index("subject_id")["sample_name"]
        )
        after_by = (
            meta[(meta["group"] == "case") & (meta["timepoint"] == "after")][["subject_id", "sample_name"]]
            .drop_duplicates()
            .set_index("subject_id")["sample_name"]
        )
        common = sorted(set(before_by.index) & set(after_by.index))
        diffs = []
        for sid in common:
            diffs.append(float(vals[after_by.loc[sid]] - vals[before_by.loc[sid]]))
        diffs_arr = np.asarray(diffs, dtype=float)
        diffs_arr = diffs_arr[np.isfinite(diffs_arr)]
        p_paired = (
            float(stats.ttest_1samp(diffs_arr, popmean=0.0, nan_policy="omit").pvalue) if diffs_arr.size else np.nan
        )

        rows.append(
            {
                "mirna": mir,
                "present": True,
                "mean_control_log2": float(np.nanmean(ctrl)) if len(ctrl) else np.nan,
                "mean_before_log2": float(np.nanmean(bef)) if len(bef) else np.nan,
                "mean_after_log2": float(np.nanmean(aft)) if len(aft) else np.nan,
                "delta_after_minus_before_log2": float(np.nanmean(diffs_arr)) if diffs_arr.size else np.nan,
                "pvalue_paired_after_vs_before": p_paired,
            }
        )

    key_df = pd.DataFrame(rows)

    # Plot: paired before/after for cases + control distribution (no bar charts).
    before_map = (
        meta[(meta["group"] == "case") & (meta["timepoint"] == "before")][["subject_id", "sample_name"]]
        .drop_duplicates()
        .set_index("subject_id")["sample_name"]
    )
    after_map = (
        meta[(meta["group"] == "case") & (meta["timepoint"] == "after")][["subject_id", "sample_name"]]
        .drop_duplicates()
        .set_index("subject_id")["sample_name"]
    )
    paired_subjects = sorted(set(before_map.index) & set(after_map.index))

    n_panels = max(1, len(key_mirnas))
    ncols = 3 if n_panels >= 3 else n_panels
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10.8, 3.4 * nrows), sharey=False)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.reshape(-1)

    colors = {"control": "#4C72B0", "before": "#DD8452", "after": "#55A868"}
    x_control, x_before, x_after = 0.0, 1.0, 2.0

    for i, mir in enumerate(key_mirnas):
        ax = axes[i]
        ax.set_title(mir, loc="left", fontsize=10)
        if mir not in expr_log2.index:
            ax.text(0.5, 0.5, "missing", transform=ax.transAxes, ha="center", va="center", color="#888888")
            ax.set_axis_off()
            continue
        vals = expr_log2.loc[mir]
        ctrl_names = meta.loc[meta["group"] == "control", "sample_name"].tolist()
        ctrl = vals[ctrl_names].astype(float).to_numpy()
        ctrl = ctrl[np.isfinite(ctrl)]

        # Control distribution (jittered points)
        if ctrl.size:
            rng = np.random.default_rng(12345)
            jitter = rng.normal(0, 0.06, size=ctrl.size)
            ax.scatter(
                x_control + jitter,
                ctrl,
                s=18,
                color=colors["control"],
                alpha=0.75,
                linewidths=0,
                label="Control" if i == 0 else None,
            )

        # Paired case before/after lines
        ys_bef = []
        ys_aft = []
        for sid in paired_subjects:
            yb = float(vals[before_map.loc[sid]])
            ya = float(vals[after_map.loc[sid]])
            if not (np.isfinite(yb) and np.isfinite(ya)):
                continue
            ys_bef.append(yb)
            ys_aft.append(ya)
            ax.plot(
                [x_before, x_after],
                [yb, ya],
                color="#808080",
                alpha=0.45,
                lw=1.0,
                zorder=1,
            )
        if ys_bef:
            ax.scatter([x_before] * len(ys_bef), ys_bef, s=20, color=colors["before"], alpha=0.80, linewidths=0, label="Before" if i == 0 else None)
        if ys_aft:
            ax.scatter([x_after] * len(ys_aft), ys_aft, s=20, color=colors["after"], alpha=0.80, linewidths=0, label="After" if i == 0 else None)

        ax.set_xticks([x_control, x_before, x_after])
        ax.set_xticklabels(["Control", "Before", "After"], rotation=0)
        ax.set_ylabel("log2(expression + 1)")
        ax.grid(axis="y", linestyle="--", alpha=0.18)
        ax.set_axisbelow(True)

    for j in range(len(key_mirnas), axes.size):
        axes[j].set_axis_off()

    fig.suptitle("GSE198274 key miRNAs (control vs paired before/after acupuncture)", fontsize=12, y=0.99)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.975))
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return key_df

def main() -> int:
    parser = argparse.ArgumentParser(description="GSE198274 acupuncture intervention miRNA response analysis.")
    parser.add_argument("--soft", default="data/raw/geo/GSE198274/soft/GSE198274_family.soft.gz", help="GSE198274 SOFT gz path.")
    parser.add_argument(
        "--expr",
        default="data/raw/geo/GSE198274/suppl/GSE198274_All_miRNA.DEG.xls.gz",
        help="GSE198274 miRNA expression table (gz; tab-delimited despite .xls extension).",
    )
    parser.add_argument("--out-meta", default="data/processed/gse198274_sample_metadata.tsv", help="Output sample metadata TSV.")
    parser.add_argument("--out-table", default="results/tables/T2.8_gse198274_mirna_response.csv", help="Output miRNA response table CSV.")
    parser.add_argument(
        "--fig-volcano",
        default="figures/raw_plots/Fig3D.png",
        help="Volcano plot (case-before vs control) for Figure 3D.",
    )
    parser.add_argument(
        "--fig-volcano-after",
        default="figures/raw_plots/Fig3I_gse198274_miRNA_volcano.pdf",
        help="Volcano plot (after vs before) for supplementary use.",
    )
    parser.add_argument("--out-key-table", default="results/tables/T2.8_gse198274_key_mirnas.csv", help="Output key miRNA summary CSV.")
    parser.add_argument("--fig-key", default="figures/raw_plots/Fig3J_gse198274_key_miRNAs.pdf", help="Key miRNA panel PDF.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    args = parser.parse_args()

    configure_logging(args.verbose)
    retrieved_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    LOGGER.info("Parsing SOFT: %s", args.soft)
    samples = parse_soft_samples(args.soft)
    meta_soft = pd.DataFrame([s.__dict__ for s in samples])
    meta_soft["retrieved_at"] = retrieved_at
    os.makedirs(os.path.dirname(args.out_meta), exist_ok=True)
    meta_soft.to_csv(args.out_meta, sep="\t", index=False)
    LOGGER.info("Saved SOFT-derived metadata: %s (n=%s)", args.out_meta, len(meta_soft))

    LOGGER.info("Loading miRNA expression table: %s", args.expr)
    expr_df = load_expression_table(args.expr)
    sample_cols = [c for c in expr_df.columns if c != "mirna"]
    if len(sample_cols) != 30:
        LOGGER.warning("Expected 30 sample columns; got %s", len(sample_cols))

    sample_meta = classify_columns(sample_cols)
    # Verify groups
    group_counts = sample_meta.groupby(["group", "timepoint"], dropna=False).size().reset_index(name="n")
    LOGGER.info("Parsed sample groups from expression table:\n%s", group_counts.to_string(index=False))

    x = expr_df[sample_cols].to_numpy(dtype=float)
    x = np.log2(np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0) + 1.0)
    expr_log2 = pd.DataFrame(x, index=expr_df["mirna"].astype(str).tolist(), columns=sample_cols)

    # Indices
    idx_control = sample_meta.index[(sample_meta["group"] == "control")].to_numpy()
    idx_before = sample_meta.index[(sample_meta["group"] == "case") & (sample_meta["timepoint"] == "before")].to_numpy()
    idx_after = sample_meta.index[(sample_meta["group"] == "case") & (sample_meta["timepoint"] == "after")].to_numpy()

    # Pairing by subject_id (Axxx)
    before_subjects = sample_meta.loc[idx_before, ["subject_id", "sample_name"]].copy()
    after_subjects = sample_meta.loc[idx_after, ["subject_id", "sample_name"]].copy()
    before_subjects = before_subjects.set_index("subject_id")
    after_subjects = after_subjects.set_index("subject_id")
    common_subjects = sorted(set(before_subjects.index) & set(after_subjects.index))
    if not common_subjects:
        raise RuntimeError("No common subjects between before and after; cannot run paired test.")

    idx_before_paired = np.array([sample_cols.index(before_subjects.loc[s, "sample_name"]) for s in common_subjects], dtype=int)
    idx_after_paired = np.array([sample_cols.index(after_subjects.loc[s, "sample_name"]) for s in common_subjects], dtype=int)

    rows: List[Dict[str, object]] = []

    # Contrast 1: before vs control (unpaired; disease baseline)
    log2fc, pvals = compute_contrast_unpaired(x, idx_before, idx_control, "case_before", "control")
    padj = _benjamini_hochberg(pvals)
    for mirna, lfc, pv, qv in zip(expr_df["mirna"].tolist(), log2fc, pvals, padj):
        rows.append(
            {
                "mirna": mirna,
                "contrast": "case_before_vs_control",
                "test": "welch_t",
                "n_group_a": int(len(idx_before)),
                "n_group_b": int(len(idx_control)),
                "log2fc": float(lfc),
                "pvalue": float(pv) if np.isfinite(pv) else np.nan,
                "padj": float(qv) if np.isfinite(qv) else np.nan,
                "direction": "up_in_case_before" if np.isfinite(lfc) and lfc > 0 else ("down_in_case_before" if np.isfinite(lfc) else "unknown"),
            }
        )

    # Contrast 2: after vs before (paired; acupuncture response)
    log2fc2, pvals2 = compute_contrast_paired(x, idx_before_paired, idx_after_paired)
    padj2 = _benjamini_hochberg(pvals2)
    for mirna, lfc, pv, qv in zip(expr_df["mirna"].tolist(), log2fc2, pvals2, padj2):
        rows.append(
            {
                "mirna": mirna,
                "contrast": "case_after_vs_before",
                "test": "paired_t",
                "n_group_a": int(len(idx_before_paired)),
                "n_group_b": int(len(idx_after_paired)),
                "log2fc": float(lfc),
                "pvalue": float(pv) if np.isfinite(pv) else np.nan,
                "padj": float(qv) if np.isfinite(qv) else np.nan,
                "direction": "up_after" if np.isfinite(lfc) and lfc > 0 else ("down_after" if np.isfinite(lfc) else "unknown"),
            }
        )

    # Contrast 3: after vs control (unpaired; post-treatment vs healthy)
    log2fc3, pvals3 = compute_contrast_unpaired(x, idx_after, idx_control, "case_after", "control")
    padj3 = _benjamini_hochberg(pvals3)
    for mirna, lfc, pv, qv in zip(expr_df["mirna"].tolist(), log2fc3, pvals3, padj3):
        rows.append(
            {
                "mirna": mirna,
                "contrast": "case_after_vs_control",
                "test": "welch_t",
                "n_group_a": int(len(idx_after)),
                "n_group_b": int(len(idx_control)),
                "log2fc": float(lfc),
                "pvalue": float(pv) if np.isfinite(pv) else np.nan,
                "padj": float(qv) if np.isfinite(qv) else np.nan,
                "direction": "up_in_case_after" if np.isfinite(lfc) and lfc > 0 else ("down_in_case_after" if np.isfinite(lfc) else "unknown"),
            }
        )

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out_table), exist_ok=True)
    out.to_csv(args.out_table, index=False)
    LOGGER.info("Saved miRNA response table: %s (rows=%s)", args.out_table, len(out))

    # Volcano for the disease-baseline contrast (before vs control)
    v_before = out[out["contrast"] == "case_before_vs_control"].copy()
    volcano_plot(
        v_before,
        args.fig_volcano,
        title="GSE198274: Before vs Control (Welch t-test)",
        x_label="log2FC (case_before − control; log2-expression)",
    )
    LOGGER.info("Saved volcano: %s", args.fig_volcano)

    # Volcano for the acupuncture response contrast (after vs before)
    v_after = out[out["contrast"] == "case_after_vs_before"].copy()
    volcano_plot(
        v_after,
        args.fig_volcano_after,
        title="GSE198274: After vs Before acupuncture (paired t-test)",
        x_label="log2FC (after − before; log2-expression)",
    )
    LOGGER.info("Saved volcano (after vs before): %s", args.fig_volcano_after)

    key_mirnas = ["hsa-miR-369-5p", "hsa-miR-1268b", "hsa-miR-145-5p", "hsa-miR-222-5p", "hsa-miR-4488"]
    key_df = key_mirna_panel(expr_log2, sample_meta, key_mirnas=key_mirnas, output_pdf=args.fig_key)
    os.makedirs(os.path.dirname(args.out_key_table), exist_ok=True)
    key_df.to_csv(args.out_key_table, index=False)
    LOGGER.info("Saved key miRNA summary: %s, %s", args.out_key_table, args.fig_key)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
