#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import os
import textwrap
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from scipy import stats


LOGGER = logging.getLogger("cross_omics_integration")
GPROFILER_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


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

def _apply_mpl_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
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


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.size, np.nan, dtype=float)
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


def read_mirna_response(path_csv: str) -> pd.DataFrame:
    df = pd.read_csv(path_csv)
    required = {"mirna", "contrast", "log2fc", "pvalue", "padj", "direction"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"miRNA response table missing columns: {sorted(missing)}")
    df["mirna"] = df["mirna"].astype(str)
    df["contrast"] = df["contrast"].astype(str)
    return df


def select_mirnas(
    df: pd.DataFrame,
    contrast: str,
    top_n: int,
    fdr: Optional[float],
    min_abs_log2fc: float,
) -> pd.DataFrame:
    sub = df[df["contrast"] == contrast].copy()
    if sub.empty:
        raise ValueError(f"Contrast not found in miRNA response table: {contrast}")

    sub = sub.dropna(subset=["pvalue", "log2fc"])
    sub = sub[np.isfinite(sub["pvalue"].to_numpy(dtype=float))]
    sub = sub[np.isfinite(sub["log2fc"].to_numpy(dtype=float))]
    if min_abs_log2fc > 0:
        sub = sub[sub["log2fc"].abs() >= float(min_abs_log2fc)]

    if sub.empty:
        raise ValueError("No miRNAs remain after filtering (check min_abs_log2fc / missing values).")

    if fdr is not None:
        sig = sub.dropna(subset=["padj"]).copy()
        sig = sig[np.isfinite(sig["padj"].to_numpy(dtype=float)) & (sig["padj"] <= float(fdr))]
        if not sig.empty:
            out = sig.sort_values(["padj", "pvalue", "mirna"]).head(int(top_n)).reset_index(drop=True)
            LOGGER.info("Selected %s miRNAs at FDR<=%s for %s.", len(out), fdr, contrast)
            return out

    out = sub.sort_values(["pvalue", "mirna"]).head(int(top_n)).reset_index(drop=True)
    LOGGER.info("Selected top %s miRNAs by p-value (fallback) for %s.", len(out), contrast)
    return out


def read_human_deg(path_csv: str) -> pd.DataFrame:
    df = pd.read_csv(path_csv)
    required = {"gene_symbol", "logFC", "P.Value", "adj.P.Val"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Human DEG table missing columns: {sorted(missing)}")
    df["gene_symbol"] = df["gene_symbol"].astype(str)
    return df


def select_disease_gene_set(df: pd.DataFrame, fdr: float, top_n_fallback: int) -> Tuple[Set[str], str]:
    df = df.dropna(subset=["gene_symbol", "P.Value", "adj.P.Val"])
    df = df[np.isfinite(df["P.Value"].to_numpy(dtype=float))]
    df = df[np.isfinite(df["adj.P.Val"].to_numpy(dtype=float))]

    sig = df[df["adj.P.Val"] <= float(fdr)]
    if not sig.empty:
        return set(sig["gene_symbol"].tolist()), f"FDR<= {fdr:g}"

    top = df.sort_values(["P.Value", "gene_symbol"]).head(int(top_n_fallback))
    return set(top["gene_symbol"].tolist()), f"top_{top_n_fallback}_by_P.Value (no FDR<= {fdr:g})"


@dataclass(frozen=True)
class MirtarBaseRow:
    mirna: str
    mirna_species: str
    target_gene: str
    target_species: str
    experiments: str
    support_type: str
    references_pmid: str


def read_mirtarbase_strong_csv(path_csv: str) -> pd.DataFrame:
    df = pd.read_csv(path_csv)
    expected = {
        "miRNA",
        "Species (miRNA)",
        "Target Gene",
        "Species (Target Gene)",
        "Experiments",
        "Support Type",
        "References (PMID)",
    }
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"miRTarBase CSV missing columns: {sorted(missing)}")
    out = df.rename(
        columns={
            "miRNA": "mirna",
            "Species (miRNA)": "mirna_species",
            "Target Gene": "target_gene",
            "Species (Target Gene)": "target_species",
            "Experiments": "experiments",
            "Support Type": "support_type",
            "References (PMID)": "pmid",
        }
    ).copy()
    out["mirna"] = out["mirna"].astype(str)
    out["mirna_species"] = out["mirna_species"].astype(str)
    out["target_gene"] = out["target_gene"].astype(str)
    out["target_species"] = out["target_species"].astype(str)
    out["experiments"] = out["experiments"].astype(str)
    out["support_type"] = out["support_type"].astype(str)
    out["pmid"] = out["pmid"].astype(str)
    return out


def build_mirna_target_edges(
    mirtar_df: pd.DataFrame,
    mirnas: Iterable[str],
    species: str,
) -> pd.DataFrame:
    mirna_set = set(str(m) for m in mirnas)
    df = mirtar_df[
        (mirtar_df["mirna_species"] == species)
        & (mirtar_df["target_species"] == species)
        & (mirtar_df["mirna"].isin(mirna_set))
    ].copy()

    if df.empty:
        raise RuntimeError(f"No miRTarBase interactions found for selected miRNAs (species={species}).")

    def _collapse(series: pd.Series) -> str:
        vals = [v for v in series.astype(str).tolist() if v and v.lower() != "nan"]
        uniq = sorted(set(vals))
        return " // ".join(uniq)

    agg = (
        df.groupby(["mirna", "target_gene"], as_index=False)
        .agg(
            n_records=("pmid", "size"),
            pmids=("pmid", _collapse),
            support_types=("support_type", _collapse),
            experiments=("experiments", _collapse),
        )
        .sort_values(["mirna", "n_records"], ascending=[True, False])
        .reset_index(drop=True)
    )
    return agg


def fisher_enrichment(
    universe: Set[str],
    disease_set: Set[str],
    target_set: Set[str],
) -> Tuple[float, int, int]:
    u = universe
    d = disease_set & u
    t = target_set & u
    a = len(d & t)
    b = len(t - d)
    c = len(d - t)
    d0 = len(u - (d | t))
    _, p = stats.fisher_exact([[a, b], [c, d0]], alternative="greater")
    return float(p), a, len(t)


def draw_bipartite_network(
    edges: pd.DataFrame,
    mirna_direction: Dict[str, str],
    human_deg: pd.DataFrame,
    output_pdf: str,
    max_edges: int,
) -> None:
    import matplotlib.pyplot as plt
    import networkx as nx

    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    _apply_mpl_style()

    if edges.empty:
        raise ValueError("No edges to plot.")

    plot_edges = edges.copy()
    if len(plot_edges) > int(max_edges):
        plot_edges = plot_edges.sort_values(["n_records", "target_in_disease_set"], ascending=[False, False]).head(int(max_edges))

    g = nx.DiGraph()
    for _, row in plot_edges.iterrows():
        g.add_edge(str(row["mirna"]), str(row["target_gene"]))

    mirnas = sorted({u for u, _ in g.edges()})
    genes = sorted({v for _, v in g.edges()})

    pos: Dict[str, Tuple[float, float]] = {}
    for i, m in enumerate(mirnas):
        pos[m] = (0.0, float(i))
    for j, ge in enumerate(genes):
        pos[ge] = (1.0, float(j))

    fig_h = max(4.0, 0.18 * max(len(mirnas), len(genes)))
    fig, ax = plt.subplots(1, 1, figsize=(9.6, fig_h))

    mirna_colors = []
    for m in mirnas:
        d = mirna_direction.get(m, "")
        if "up" in d:
            mirna_colors.append("#2CA02C")
        elif "down" in d:
            mirna_colors.append("#D62728")
        else:
            mirna_colors.append("#7F7F7F")

    gene_df = human_deg.set_index("gene_symbol")
    gene_colors = []
    for ge in genes:
        if ge in gene_df.index:
            lfc = float(gene_df.loc[ge, "logFC"])
            gene_colors.append("#C44E52" if lfc > 0 else "#4C72B0")
        else:
            gene_colors.append("#9AA0A6")

    nx.draw_networkx_edges(g, pos=pos, ax=ax, alpha=0.25, width=0.8, arrows=False)
    nx.draw_networkx_nodes(
        g,
        pos=pos,
        nodelist=mirnas,
        node_color=mirna_colors,
        node_size=340,
        ax=ax,
        linewidths=0.8,
        edgecolors="white",
    )
    nx.draw_networkx_nodes(
        g,
        pos=pos,
        nodelist=genes,
        node_color=gene_colors,
        node_size=230,
        ax=ax,
        linewidths=0.8,
        edgecolors="white",
    )

    for m in mirnas:
        ax.text(pos[m][0] - 0.03, pos[m][1], m, ha="right", va="center", fontsize=9, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.2))
    for ge in genes:
        ax.text(pos[ge][0] + 0.03, pos[ge][1], ge, ha="left", va="center", fontsize=9, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.2))

    ax.set_title("Candidate miRNA→gene axes (miRTarBase + disease background)", loc="left", fontsize=12, fontweight="bold")
    ax.set_axis_off()
    fig.tight_layout(pad=1.0)
    # Ensure png output if requested
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _require_requests() -> "object":
    try:
        import requests  # type: ignore

        return requests
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "Missing dependency: requests. Install with `python3 -m pip install requests` and re-run."
        ) from exc


def gprofiler_profile(
    genes: Sequence[str],
    organism: str,
    sources: Sequence[str],
    user_threshold: float,
    significance_threshold_method: str,
    no_iea: bool,
    timeout_s: int,
) -> dict:
    requests = _require_requests()
    payload = {
        "organism": str(organism),
        "query": list(map(str, genes)),
        "sources": list(map(str, sources)),
        "user_threshold": float(user_threshold),
        # g:Profiler API uses this name (defaults to g_SCS if not provided).
        "significance_threshold_method": str(significance_threshold_method),
        "no_iea": bool(no_iea),
    }
    resp = requests.post(GPROFILER_ENDPOINT, json=payload, timeout=int(timeout_s))
    resp.raise_for_status()
    return resp.json()


def normalise_gprofiler_results(raw: dict) -> pd.DataFrame:
    # g:Profiler returns { result: [ {...}, ... ], meta: {...} }
    rows = raw.get("result", []) if isinstance(raw, dict) else []
    if not rows:
        return pd.DataFrame(
            columns=[
                "source",
                "term_id",
                "term_name",
                "p_value",
                "p_value_adjusted",
                "intersection_size",
                "term_size",
                "query_size",
                "effective_domain_size",
                "intersection",
            ]
        )
    df = pd.DataFrame(rows).copy()

    query_genes: List[str] = []
    try:
        queries = raw.get("meta", {}).get("query_metadata", {}).get("queries", {})
        if isinstance(queries, dict) and queries:
            first_key = sorted(queries.keys())[0]
            q = queries.get(first_key, [])
            if isinstance(q, list):
                query_genes = [str(x) for x in q]
    except Exception:
        query_genes = []

    keep = {
        "source",
        "native",
        "name",
        "p_value",
        "p_value_adjusted",
        "intersection_size",
        "term_size",
        "query_size",
        "effective_domain_size",
        "intersection",
        "intersections",
    }
    df = df[[c for c in df.columns if c in keep]].copy()
    df = df.rename(columns={"native": "term_id", "name": "term_name"})
    for col in ["p_value", "p_value_adjusted"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "p_value_adjusted" not in df.columns and "p_value" in df.columns:
        # Some g:Profiler deployments return only p_value which is already corrected
        # according to meta.query_metadata.significance_threshold_method (default g_SCS).
        df["p_value_adjusted"] = pd.to_numeric(df["p_value"], errors="coerce")

    for col in ["intersection_size", "term_size", "query_size", "effective_domain_size"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Derive overlap gene names from per-query "intersections" evidence arrays.
    if "intersections" in df.columns and query_genes:
        def _overlap_genes(x: object) -> str:
            if not isinstance(x, list):
                return ""
            genes = []
            for i, ev in enumerate(x):
                if i >= len(query_genes):
                    break
                if isinstance(ev, list) and len(ev) > 0:
                    genes.append(query_genes[i])
            return ";".join(genes)

        df["intersection"] = df["intersections"].apply(_overlap_genes)
    elif "intersection" in df.columns:
        df["intersection"] = df["intersection"].apply(lambda x: ";".join(map(str, x)) if isinstance(x, list) else str(x))
    else:
        df["intersection"] = ""

    if "intersections" in df.columns:
        df = df.drop(columns=["intersections"])
    df = df.sort_values(["source", "p_value_adjusted", "p_value", "term_id"], na_position="last").reset_index(drop=True)
    return df


def plot_enrichment_bar(
    df: pd.DataFrame,
    output_pdf: str,
    max_terms_per_source: int,
) -> None:
    import matplotlib.pyplot as plt

    if df.empty:
        raise RuntimeError("No enrichment terms to plot (empty g:Profiler results).")

    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    _apply_mpl_style()

    plot_df = df.dropna(subset=["p_value_adjusted"]).copy()
    plot_df = plot_df[np.isfinite(plot_df["p_value_adjusted"].to_numpy(dtype=float))].copy()
    plot_df = plot_df.sort_values(["source", "p_value_adjusted", "p_value", "term_id"], na_position="last").copy()
    plot_df = plot_df.groupby("source", as_index=False).head(int(max_terms_per_source)).copy()
    plot_df["neglog10_fdr"] = -np.log10(np.clip(plot_df["p_value_adjusted"].astype(float), 1e-300, 1.0))
    plot_df["gene_ratio"] = plot_df["intersection_size"].astype(float) / np.clip(
        plot_df["query_size"].astype(float), 1.0, np.inf
    )

    def _wrap_label(v: str, width: int = 34, max_chars: int = 120) -> str:
        s = str(v or "")
        if len(s) > max_chars:
            s = s[: max_chars - 1] + "…"
        return textwrap.fill(s, width=width)

    plot_df["label"] = plot_df["term_name"].map(_wrap_label)

    sources = list(plot_df["source"].dropna().unique())
    sources = sorted(sources)
    n = len(sources)
    total_terms = int(plot_df.shape[0])
    fig_h = max(6.5, 0.30 * total_terms + 0.9 * n)
    fig, axes = plt.subplots(n, 1, figsize=(12.2, fig_h), sharex=False)
    if n == 1:
        axes = [axes]

    cmap = plt.get_cmap("viridis")
    vmin = float(plot_df["neglog10_fdr"].min()) if not plot_df.empty else 0.0
    vmax = float(plot_df["neglog10_fdr"].max()) if not plot_df.empty else 1.0

    for ax, src in zip(axes, sources):
        sub = plot_df[plot_df["source"] == src].copy()
        sub = sub.sort_values(["neglog10_fdr", "gene_ratio"], ascending=[True, True]).reset_index(drop=True)
        y = np.arange(sub.shape[0])
        colors = cmap((sub["neglog10_fdr"].to_numpy(dtype=float) - vmin) / max(vmax - vmin, 1e-9))
        sizes = 40 + 90 * (sub["intersection_size"].to_numpy(dtype=float) / max(sub["intersection_size"].max(), 1.0))
        ax.scatter(sub["gene_ratio"], y, s=sizes, c=colors, edgecolors="white", linewidths=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["label"].tolist(), fontsize=7)
        ax.set_xlabel("GeneRatio (overlap / query)")
        ax.set_title(f"{src} (top {max_terms_per_source} by FDR)", loc="left", fontsize=10)
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        ax.set_axisbelow(True)

    fig.suptitle("Pathway enrichment for miRNA→mRNA mechanism axes (g:Profiler ORA)", fontsize=12)
    fig.subplots_adjust(left=0.40, right=0.96, top=0.93, bottom=0.06, hspace=0.70)
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.965, 0.20, 0.015, 0.55])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("-log10(FDR)")

    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_mechanism_note(
    axes: pd.DataFrame,
    enr: pd.DataFrame,
    out_md: str,
    run_at: str,
    disease_deg_path: str,
    mirna_response_path: str,
    mirtarbase_path: str,
    max_terms_per_source: int,
) -> None:
    os.makedirs(os.path.dirname(out_md), exist_ok=True)
    genes = sorted(set(axes["target_gene"].astype(str).tolist())) if not axes.empty else []
    lines: List[str] = []
    lines.append("# Mechanism Note (minimal)\n")
    lines.append(f"- Run at (UTC): {run_at}")
    lines.append(f"- Axes input (miRNA→gene): {len(axes)} edges; {len(genes)} unique target genes")
    lines.append(f"- Disease background DEG: `{disease_deg_path}`")
    lines.append(f"- Acupuncture miRNA response: `{mirna_response_path}`")
    lines.append(f"- miRNA→target mapping: `{mirtarbase_path}` (miRTarBase strong evidence)\n")

    if not enr.empty:
        lines.append("## Enrichment Summary (g:Profiler ORA)\n")
        top = (
            enr.dropna(subset=["p_value_adjusted"])
            .sort_values(["source", "p_value_adjusted", "p_value"])
            .groupby("source", as_index=False)
            .head(int(max_terms_per_source))
        )
        for src in sorted(top["source"].unique()):
            sub = top[top["source"] == src].copy()
            if sub.empty:
                continue
            lines.append(f"### {src}")
            for _, r in sub.iterrows():
                term = str(r.get("term_name", ""))
                tid = str(r.get("term_id", ""))
                fdr = r.get("p_value_adjusted", float("nan"))
                inter = str(r.get("intersection", ""))
                inter_genes = [g for g in inter.split(";") if g]
                inter_genes = inter_genes[:15]
                lines.append(f"- {term} ({tid}) — FDR={fdr:.3g}; overlap_genes={', '.join(inter_genes)}")
            lines.append("")
    else:
        lines.append("## Enrichment Summary\n")
        lines.append("- No enrichment results returned (empty).")

    lines.append("## Candidate Axes (miRNA → gene)\n")
    if axes.empty:
        lines.append("- No axes were supported by the disease gene set in this run.\n")
    else:
        preview = axes.sort_values(["mirna", "target_fdr", "target_pvalue"], na_position="last").copy()
        for _, r in preview.head(30).iterrows():
            lines.append(f"- {r['mirna']} → {r['target_gene']} (target_FDR={r.get('target_fdr', float('nan')):.3g})")
        if len(preview) > 30:
            lines.append(f"- … ({len(preview) - 30} more)\n")
        else:
            lines.append("")

    content = "\n".join(lines).rstrip() + "\n"
    Path = __import__("pathlib").Path
    Path(out_md).write_text(content, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Task 2.9: Cross-omics integration (miRNA → targets → human mRNA).")
    parser.add_argument(
        "--mirna-response",
        default="results/tables/T2.8_gse198274_mirna_response.csv",
        help="miRNA response table (Task 2.8).",
    )
    parser.add_argument("--mirna-contrast", default="case_after_vs_before", help="Which miRNA contrast to treat as acupuncture response.")
    parser.add_argument("--top-mirnas", type=int, default=25, help="Number of top miRNAs to carry into integration.")
    parser.add_argument("--mirna-fdr", type=float, default=None, help="If provided, prefer miRNAs with padj<=FDR; fallback to top by p-value.")
    parser.add_argument("--min-abs-log2fc", type=float, default=0.0, help="Minimum |log2fc| for selected miRNAs.")
    parser.add_argument(
        "--mirna-prefix",
        default="hsa-",
        help="If set (default: hsa-), restrict miRNA selection to names starting with this prefix to ensure target mapping.",
    )
    parser.add_argument(
        "--require-mirtarbase-mapping",
        action="store_true",
        help="If set, restrict miRNA selection to those present in the miRTarBase file after species filtering.",
    )

    parser.add_argument(
        "--mirtarbase",
        default="data/raw/mirtarbase/10.0/miRTarBase_SE_WR.csv",
        help="miRTarBase mapping CSV (strong evidence set).",
    )
    parser.add_argument("--species", default="hsa", help="Species code to filter (default: hsa).")

    parser.add_argument("--human-deg", default="results/tables/T2.6_human_deg_emtab13397.csv", help="Human DEG table (Task 2.6).")
    parser.add_argument("--human-fdr", type=float, default=0.05, help="Human gene-set FDR threshold (fallback to top by p-value).")
    parser.add_argument(
        "--human-pvalue",
        type=float,
        default=0.05,
        help="If no genes pass FDR, define the disease-background gene set as P.Value<=threshold before falling back to top-N.",
    )
    parser.add_argument("--human-top-genes", type=int, default=500, help="Fallback human disease-gene set size (top by p-value).")

    parser.add_argument("--out-edges", default="results/tables/T2.9_mirna_target_edges.tsv", help="Output edge list TSV (miRNA→gene).")
    parser.add_argument("--out-axes", default="results/tables/T2.9_mechanism_axes.tsv", help="Output filtered axes TSV (targets supported by human DEG set).")
    parser.add_argument("--out-enrichment", default="results/tables/T2.9_mirna_target_enrichment.csv", help="Per-miRNA enrichment CSV.")
    parser.add_argument("--fig-network", default="figures/raw_plots/Fig3E.png", help="Draft network figure (PNG).")
    parser.add_argument("--max-edges-plot", type=int, default=180, help="Maximum edges to plot in the draft network.")
    parser.add_argument("--out-config", default="results/tables/T2.9_integration_config.json", help="Run config JSON output.")
    parser.add_argument("--run-pathway-enrichment", action="store_true", help="If set, run g:Profiler ORA on the mechanism axes target genes.")
    parser.add_argument(
        "--enrichment-sources",
        default="GO:BP,GO:CC,GO:MF,REAC,KEGG",
        help="Comma-separated g:Profiler sources (e.g., GO:BP,REAC,KEGG).",
    )
    parser.add_argument("--enrichment-organism", default="hsapiens", help="g:Profiler organism (default: hsapiens).")
    parser.add_argument("--enrichment-user-threshold", type=float, default=0.05, help="g:Profiler user_threshold (default: 0.05).")
    parser.add_argument(
        "--enrichment-significance-threshold-method",
        "--enrichment-correction-method",
        dest="enrichment_significance_threshold_method",
        default="g_SCS",
        choices=["g_SCS", "fdr", "bonferroni"],
        help="Multiple testing correction method (g:Profiler API name: significance_threshold_method).",
    )
    parser.add_argument("--enrichment-no-iea", action="store_true", help="If set, exclude GO IEA annotations (no_iea=true).")
    parser.add_argument("--enrichment-timeout", type=int, default=60, help="HTTP timeout seconds for g:Profiler API.")
    parser.add_argument("--out-pathway-terms", default="", help="Output CSV for pathway enrichment terms (axes gene set).")
    parser.add_argument("--out-pathway-raw-json", default="", help="Output raw g:Profiler JSON for audit.")
    parser.add_argument("--fig-pathway-enrichment", default="", help="Output PDF barplots for pathway enrichment.")
    parser.add_argument("--out-mechanism-note", default="", help="Output markdown mechanism note summary.")
    parser.add_argument("--max-terms-per-source", type=int, default=10, help="Max terms per source in plots/narrative.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    args = parser.parse_args()

    configure_logging(args.verbose)
    run_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    LOGGER.info("Loading miRTarBase mapping: %s", args.mirtarbase)
    mtb = read_mirtarbase_strong_csv(args.mirtarbase)
    available_mirnas = set(mtb.loc[mtb["mirna_species"] == args.species, "mirna"].astype(str).tolist())

    LOGGER.info("Loading miRNA response: %s", args.mirna_response)
    mir_df = read_mirna_response(args.mirna_response)
    if args.mirna_prefix:
        mir_df = mir_df[mir_df["mirna"].astype(str).str.startswith(str(args.mirna_prefix))].copy()
    if args.require_mirtarbase_mapping:
        mir_df = mir_df[mir_df["mirna"].isin(available_mirnas)].copy()
    if mir_df.empty:
        raise RuntimeError("No miRNAs remain after applying mirna-prefix / require-mirtarbase-mapping filters.")

    sel_mir = select_mirnas(
        mir_df,
        contrast=args.mirna_contrast,
        top_n=args.top_mirnas,
        fdr=args.mirna_fdr,
        min_abs_log2fc=args.min_abs_log2fc,
    )
    selected_mirnas = sel_mir["mirna"].tolist()
    mirna_direction = dict(zip(sel_mir["mirna"], sel_mir["direction"]))
    LOGGER.info("Selected miRNAs (n=%s): %s", len(selected_mirnas), ", ".join(selected_mirnas))

    edges = build_mirna_target_edges(mtb, selected_mirnas, species=args.species)
    LOGGER.info("Raw edges from miRTarBase (selected miRNAs): %s", len(edges))

    LOGGER.info("Loading human DEG table: %s", args.human_deg)
    human = read_human_deg(args.human_deg)
    universe = set(human["gene_symbol"].astype(str).tolist())

    disease_set, disease_rule = select_disease_gene_set(human, fdr=args.human_fdr, top_n_fallback=args.human_top_genes)
    if disease_rule.startswith("top_"):
        pv = float(args.human_pvalue)
        pv_hits = human.dropna(subset=["gene_symbol", "P.Value"]).copy()
        pv_hits = pv_hits[np.isfinite(pv_hits["P.Value"].to_numpy(dtype=float)) & (pv_hits["P.Value"] <= pv)]
        if not pv_hits.empty:
            disease_set = set(pv_hits["gene_symbol"].astype(str).tolist())
            disease_rule = f"P.Value<= {pv:g} (no FDR<= {args.human_fdr:g})"

    LOGGER.info("Disease gene-set rule: %s (n=%s; universe=%s)", disease_rule, len(disease_set), len(universe))

    # Annotate edges with human DEG stats and disease-set membership
    human_idx = human.set_index("gene_symbol")
    edges["target_in_universe"] = edges["target_gene"].isin(universe)
    edges = edges[edges["target_in_universe"]].copy()
    edges["target_in_disease_set"] = edges["target_gene"].isin(disease_set)

    def _lookup(col: str, gene: str) -> float:
        if gene not in human_idx.index:
            return float("nan")
        v = human_idx.loc[gene, col]
        try:
            return float(v)
        except Exception:
            return float("nan")

    edges["target_logFC"] = edges["target_gene"].map(lambda g: _lookup("logFC", g))
    edges["target_pvalue"] = edges["target_gene"].map(lambda g: _lookup("P.Value", g))
    edges["target_fdr"] = edges["target_gene"].map(lambda g: _lookup("adj.P.Val", g))

    # Export full edge list
    os.makedirs(os.path.dirname(args.out_edges), exist_ok=True)
    edges.to_csv(args.out_edges, sep="\t", index=False)
    LOGGER.info("Saved edges: %s", args.out_edges)

    # Filtered “axes” (only targets supported by human disease gene set)
    axes = edges[edges["target_in_disease_set"]].copy()
    axes = axes.sort_values(["mirna", "target_pvalue", "n_records"], ascending=[True, True, False]).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_axes), exist_ok=True)
    axes.to_csv(args.out_axes, sep="\t", index=False)
    LOGGER.info("Saved mechanism axes: %s (edges=%s)", args.out_axes, len(axes))

    # Per-miRNA enrichment of disease-set genes among targets
    enr_rows: List[Dict[str, object]] = []
    for mir in selected_mirnas:
        tset = set(edges.loc[edges["mirna"] == mir, "target_gene"].astype(str).tolist())
        p, overlap, n_targets = fisher_enrichment(universe, disease_set, tset)
        enr_rows.append(
            {
                "mirna": mir,
                "direction": mirna_direction.get(mir, ""),
                "n_targets_in_universe": int(n_targets),
                "n_targets_in_disease_set": int(overlap),
                "disease_set_rule": disease_rule,
                "pvalue_fisher_greater": p,
            }
        )
    enr = pd.DataFrame(enr_rows)
    enr["padj_fisher"] = benjamini_hochberg(enr["pvalue_fisher_greater"].to_numpy(dtype=float))
    enr = enr.sort_values(["padj_fisher", "pvalue_fisher_greater", "mirna"]).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_enrichment), exist_ok=True)
    enr.to_csv(args.out_enrichment, index=False)
    LOGGER.info("Saved enrichment: %s", args.out_enrichment)

    # Draft network figure (plot only axes if available; else plot the unfiltered edges)
    plot_edges = axes if not axes.empty else edges.copy()
    draw_bipartite_network(
        plot_edges,
        mirna_direction=mirna_direction,
        human_deg=human,
        output_pdf=args.fig_network,
        max_edges=args.max_edges_plot,
    )
    LOGGER.info("Saved network figure: %s", args.fig_network)

    pathway_terms_df = pd.DataFrame()
    if args.run_pathway_enrichment:
        if axes.empty:
            raise RuntimeError("Requested --run-pathway-enrichment but mechanism axes are empty (no supported targets).")

        sources = [s.strip() for s in str(args.enrichment_sources).split(",") if s.strip()]
        if not sources:
            raise ValueError("No enrichment sources parsed from --enrichment-sources.")

        axis_genes = sorted(set(axes["target_gene"].astype(str).tolist()))
        LOGGER.info("Running g:Profiler ORA for axes genes: n=%d sources=%s", len(axis_genes), ",".join(sources))
        raw = gprofiler_profile(
            axis_genes,
            organism=args.enrichment_organism,
            sources=sources,
            user_threshold=args.enrichment_user_threshold,
            significance_threshold_method=args.enrichment_significance_threshold_method,
            no_iea=args.enrichment_no_iea,
            timeout_s=args.enrichment_timeout,
        )
        pathway_terms_df = normalise_gprofiler_results(raw)

        out_terms = args.out_pathway_terms or os.path.join(os.path.dirname(args.out_edges), "T2.9_mechanism_axes_pathway_enrichment.csv")
        out_raw = args.out_pathway_raw_json or os.path.join(os.path.dirname(args.out_edges), "T2.9_mechanism_axes_pathway_enrichment_raw.json")
        out_fig = args.fig_pathway_enrichment or os.path.join(os.path.dirname(args.fig_network), "Fig3O_mechanism_axes_pathway_enrichment.pdf")
        out_note = args.out_mechanism_note or os.path.join(os.path.dirname(args.out_edges), "T2.9_mechanism_note.md")

        os.makedirs(os.path.dirname(out_terms), exist_ok=True)
        pathway_terms_df.to_csv(out_terms, index=False)
        LOGGER.info("Saved pathway enrichment terms: %s", out_terms)

        os.makedirs(os.path.dirname(out_raw), exist_ok=True)
        with open(out_raw, "w", encoding="utf-8") as file_handle:
            json.dump(raw, file_handle, ensure_ascii=False, indent=2)
            file_handle.write("\n")
        LOGGER.info("Saved pathway enrichment raw JSON: %s", out_raw)

        if not pathway_terms_df.empty:
            plot_enrichment_bar(pathway_terms_df, out_fig, max_terms_per_source=args.max_terms_per_source)
            LOGGER.info("Saved pathway enrichment figure: %s", out_fig)

        write_mechanism_note(
            axes=axes,
            enr=pathway_terms_df,
            out_md=out_note,
            run_at=run_at,
            disease_deg_path=args.human_deg,
            mirna_response_path=args.mirna_response,
            mirtarbase_path=args.mirtarbase,
            max_terms_per_source=args.max_terms_per_source,
        )
        LOGGER.info("Saved mechanism note: %s", out_note)

    config = {
        "run_at_utc": run_at,
        "inputs": {
            "mirna_response": args.mirna_response,
            "mirna_contrast": args.mirna_contrast,
            "mirtarbase": args.mirtarbase,
            "species": args.species,
            "human_deg": args.human_deg,
        },
        "parameters": {
            "top_mirnas": args.top_mirnas,
            "mirna_fdr_prefer": args.mirna_fdr,
            "min_abs_log2fc": args.min_abs_log2fc,
            "human_fdr": args.human_fdr,
            "human_top_genes_fallback": args.human_top_genes,
            "max_edges_plot": args.max_edges_plot,
            "run_pathway_enrichment": bool(args.run_pathway_enrichment),
            "enrichment_sources": str(args.enrichment_sources),
            "enrichment_organism": str(args.enrichment_organism),
            "enrichment_user_threshold": float(args.enrichment_user_threshold),
            "enrichment_significance_threshold_method": str(args.enrichment_significance_threshold_method),
            "enrichment_no_iea": bool(args.enrichment_no_iea),
            "enrichment_timeout": int(args.enrichment_timeout),
            "max_terms_per_source": int(args.max_terms_per_source),
        },
        "derived": {
            "selected_mirnas": selected_mirnas,
            "disease_set_rule": disease_rule,
            "n_edges_total": int(len(edges)),
            "n_edges_axes": int(len(axes)),
            "n_axes_genes": int(axes["target_gene"].nunique()) if not axes.empty else 0,
            "n_pathway_terms": int(len(pathway_terms_df)) if isinstance(pathway_terms_df, pd.DataFrame) else 0,
        },
        "outputs": {
            "edges_tsv": args.out_edges,
            "axes_tsv": args.out_axes,
            "enrichment_csv": args.out_enrichment,
            "network_pdf": args.fig_network,
        },
    }
    if args.run_pathway_enrichment:
        config["outputs"]["pathway_terms_csv"] = args.out_pathway_terms or os.path.join(os.path.dirname(args.out_edges), "T2.9_mechanism_axes_pathway_enrichment.csv")
        config["outputs"]["pathway_raw_json"] = args.out_pathway_raw_json or os.path.join(os.path.dirname(args.out_edges), "T2.9_mechanism_axes_pathway_enrichment_raw.json")
        config["outputs"]["pathway_fig_pdf"] = args.fig_pathway_enrichment or os.path.join(os.path.dirname(args.fig_network), "Fig3O_mechanism_axes_pathway_enrichment.pdf")
        config["outputs"]["mechanism_note_md"] = args.out_mechanism_note or os.path.join(os.path.dirname(args.out_edges), "T2.9_mechanism_note.md")
    os.makedirs(os.path.dirname(args.out_config), exist_ok=True)
    with open(args.out_config, "w", encoding="utf-8") as file_handle:
        json.dump(config, file_handle, ensure_ascii=False, indent=2)
        file_handle.write("\n")

    LOGGER.info("Saved run config: %s", args.out_config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
