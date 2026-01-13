#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import gzip
import logging
import os
import re
import sys
import textwrap
import urllib.request
from collections import defaultdict
from typing import Iterable


LOGGER = logging.getLogger("target_intersection")


DEFAULT_CTD_CURATED_URL = "https://ctdbase.org/reports/CTD_curated_genes_diseases.tsv.gz"


GENE_MAPPINGS: dict[str, list[str]] = {
    "CGRP": ["CALCA", "CALCB"],
    "PACAP": ["ADCYAP1"],
    "Substance P": ["TAC1"],
    "NPY": ["NPY"],
    "TNF": ["TNF"],
    "IL-1β": ["IL1B"],
    "IL-6": ["IL6"],
    "IL-10": ["IL10"],
    "IL-17": ["IL17A", "IL17F"],
    "ET-1": ["EDN1"],
    "TRPV1": ["TRPV1"],
    "BDNF": ["BDNF"],
    "Nitric oxide": ["NOS1", "NOS2", "NOS3"],
}


@dataclasses.dataclass(frozen=True)
class CtdGeneDiseaseRow:
    gene_symbol: str
    gene_id: str
    disease_name: str
    disease_id: str
    direct_evidence: str
    omim_ids: str
    pubmed_ids: str


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    # Keep matplotlib/font subsetting noise out of logs.
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)


def download_file(url: str, output_path: str) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    LOGGER.info("Downloading %s -> %s", url, output_path)
    with urllib.request.urlopen(url) as response, open(output_path, "wb") as file_handle:
        file_handle.write(response.read())


def read_ctd_curated(path: str) -> Iterable[CtdGeneDiseaseRow]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as file_handle:
        for line in file_handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            yield CtdGeneDiseaseRow(
                gene_symbol=parts[0],
                gene_id=parts[1],
                disease_name=parts[2],
                disease_id=parts[3],
                direct_evidence=parts[4],
                omim_ids=parts[5],
                pubmed_ids=parts[6],
            )


def load_molecules_from_literature(path: str) -> set[str]:
    molecules: set[str] = set()
    with open(path, "r", newline="", encoding="utf-8") as file_handle:
        reader = csv.DictReader(file_handle)
        for row in reader:
            molecule = (row.get("molecule") or "").strip()
            if molecule:
                molecules.add(molecule)
    return molecules


def build_acupuncture_target_genes(molecules_present: set[str]) -> dict[str, set[str]]:
    gene_to_sources: dict[str, set[str]] = defaultdict(set)
    for molecule in sorted(molecules_present):
        mapped_genes = GENE_MAPPINGS.get(molecule)
        if not mapped_genes:
            continue
        for gene in mapped_genes:
            gene_to_sources[gene].add(molecule)
    return gene_to_sources


def write_acupuncture_targets_csv(output_path: str, gene_to_sources: dict[str, set[str]]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=["gene_symbol", "source_molecules"])
        writer.writeheader()
        for gene_symbol in sorted(gene_to_sources.keys()):
            writer.writerow(
                {
                    "gene_symbol": gene_symbol,
                    "source_molecules": "|".join(sorted(gene_to_sources[gene_symbol])),
                }
            )


def write_migraine_targets_csv(output_path: str, rows: list[CtdGeneDiseaseRow]) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(
            file_handle,
            fieldnames=[
                "gene_symbol",
                "gene_id",
                "disease_name",
                "disease_id",
                "direct_evidence",
                "omim_ids",
                "pubmed_ids",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(dataclasses.asdict(row))


def write_intersection_csv(
    output_path: str,
    intersection_genes: set[str],
    gene_to_sources: dict[str, set[str]],
    ctd_by_gene: dict[str, list[CtdGeneDiseaseRow]],
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(
            file_handle,
            fieldnames=["gene_symbol", "source_molecules", "ctd_diseases", "ctd_pubmed_ids"],
        )
        writer.writeheader()
        for gene in sorted(intersection_genes):
            diseases = sorted({row.disease_name for row in ctd_by_gene.get(gene, [])})
            pubmed = sorted(
                {pmid for row in ctd_by_gene.get(gene, []) for pmid in (row.pubmed_ids.split("|") if row.pubmed_ids else [])}
            )
            writer.writerow(
                {
                    "gene_symbol": gene,
                    "source_molecules": "|".join(sorted(gene_to_sources.get(gene, set()))),
                    "ctd_diseases": "|".join(diseases),
                    "ctd_pubmed_ids": "|".join(pubmed),
                }
            )


def save_venn_plot(output_path: str, *, a: int, b: int, ab: int, label_a: str, label_b: str) -> None:
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # SOTA 2025: Clean, flat style
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 14,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        }
    )
    fig, ax = plt.subplots(figsize=(6, 5))

    # Pastel colors (Blue/Orange -> classic scientific contrast)
    v = venn2(subsets=(a, b, ab), set_labels=(label_a, label_b), set_colors=('#a6cee3', '#fdbf6f'), alpha=0.6)
    
    # Improve text visibility
    for text in v.set_labels:
        if text:
            text.set_fontsize(14)
            text.set_fontweight('bold')
    for text in v.subset_labels:
        if text:
            text.set_fontsize(16)
            text.set_color('black') # Black numbers for visibility on pastel
            text.set_fontweight('bold')

    # Panel title is added during figure board assembly; avoid redundant titles inside the panel.
    fig.tight_layout(pad=1.0)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_network_plot(
    output_path: str,
    *,
    disease_node: str,
    intersection_genes: set[str],
    gene_to_sources: dict[str, set[str]],
) -> None:
    import matplotlib.pyplot as plt
    import networkx as nx
    import math

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        }
    )
    graph = nx.Graph()

    graph.add_node(disease_node, kind="disease")
    for gene in sorted(intersection_genes):
        graph.add_node(gene, kind="gene")
        graph.add_edge(disease_node, gene, kind="gene_disease")
        for molecule in sorted(gene_to_sources.get(gene, set())):
            mol_node = f"MED::{molecule}"
            graph.add_node(mol_node, kind="molecule", display_name=molecule)
            graph.add_edge(mol_node, gene, kind="molecule_gene")

    if graph.number_of_nodes() == 0:
        raise RuntimeError("Network is empty; cannot plot.")

    # SOTA Layout: Concentric/Shell but with force adjustment
    # Use spring layout with fixed seed but high k for spacing
    pos = nx.spring_layout(graph, k=0.5, iterations=50, seed=42)

    fig, ax = plt.subplots(figsize=(8.0, 8.0))
    ax.set_axis_off()

    # Color Palette: SOTA 2025 (Flat UI colors)
    color_map = {
        "disease": "#e74c3c", # Red
        "gene": "#3498db",    # Blue
        "molecule": "#2ecc71" # Green
    }
    
    node_colors = [color_map.get(graph.nodes[n].get("kind"), "gray") for n in graph.nodes]
    node_sizes = [1500 if graph.nodes[n].get("kind")=="disease" else (1000 if graph.nodes[n].get("kind")=="gene" else 600) for n in graph.nodes]

    nx.draw_networkx_edges(graph, pos, ax=ax, alpha=0.2, width=1.5, edge_color="#95a5a6")
    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_color=node_colors,
        node_size=node_sizes,
        linewidths=1.5,
        edgecolors="white",
    )

    labels = {}
    for node in graph.nodes:
        kind = graph.nodes[node].get("kind")
        if kind == "molecule":
             labels[node] = graph.nodes[node].get("display_name")
        else:
             labels[node] = node

    # Labels with white outline for readability
    text_items = nx.draw_networkx_labels(
        graph,
        pos,
        labels=labels,
        font_size=10,
        font_weight="bold",
        font_family="sans-serif",
        font_color="#2c3e50",
        ax=ax,
    )
    # Add white outline to text
    import matplotlib.patheffects as path_effects
    for _, t in text_items.items():
        t.set_path_effects([path_effects.withStroke(linewidth=3, foreground='white')])

    # Custom Legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Disease', markerfacecolor=color_map['disease'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', label='Target Gene', markerfacecolor=color_map['gene'], markersize=12),
        plt.Line2D([0], [0], marker='o', color='w', label='Active Component', markerfacecolor=color_map['molecule'], markersize=10)
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, framealpha=0.9, fontsize=10)

    # Panel title is added during figure board assembly; avoid redundant titles inside the panel.

    fig.tight_layout(pad=1.0)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_ctd_evidence_bar_plot(
    output_path: str,
    *,
    migraine_rows: list[CtdGeneDiseaseRow],
    intersection_genes: set[str],
    top_n: int = 15,
) -> None:
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        }
    )

    by_gene: dict[str, list[CtdGeneDiseaseRow]] = defaultdict(list)
    for row in migraine_rows:
        if row.gene_symbol:
            by_gene[row.gene_symbol].append(row)

    records = []
    for gene, rows in by_gene.items():
        pmids = {pmid for r in rows for pmid in (r.pubmed_ids.split("|") if r.pubmed_ids else [])}
        records.append(
            {
                "gene": gene,
                "n_rows": len(rows),
                "n_pmids": len(pmids),
                "is_intersection": gene in intersection_genes,
            }
        )
    if not records:
        raise RuntimeError("No CTD evidence rows available to plot.")

    records = sorted(records, key=lambda r: (r["n_rows"], r["n_pmids"]), reverse=True)[: int(top_n)]
    genes = [r["gene"] for r in records][::-1]
    n_rows = [r["n_rows"] for r in records][::-1]
    colors = ["#C44E52" if r["is_intersection"] else "#4C72B0" for r in records][::-1]

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    # Lollipop (cleaner than solid bars for dense labels).
    y = list(range(len(genes)))
    ax.hlines(y=y, xmin=0, xmax=n_rows, color="#bdbdbd", lw=2.0, alpha=0.6)
    ax.scatter(n_rows, y, s=80, c=colors, edgecolors="white", linewidths=1.0, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlabel("CTD curated gene–disease evidence records (count)")
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    # Legend-free but explicit note: red bars are intersection genes.
    ax.text(
        0.99,
        0.02,
        "Red = intersection genes",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        color="#C44E52",
    )

    fig.tight_layout(pad=1.0)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_intersection_table_plot(
    output_path: str,
    *,
    migraine_rows: list[CtdGeneDiseaseRow],
    intersection_genes: set[str],
    gene_to_sources: dict[str, set[str]],
    max_rows: int = 12,
) -> None:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 10,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        }
    )

    by_gene: dict[str, list[CtdGeneDiseaseRow]] = defaultdict(list)
    for row in migraine_rows:
        if row.gene_symbol:
            by_gene[row.gene_symbol].append(row)

    rows = []
    for gene in sorted(intersection_genes):
        pmids = sorted({pmid for r in by_gene.get(gene, []) for pmid in (r.pubmed_ids.split("|") if r.pubmed_ids else [])})
        rows.append(
            [
                gene,
                str(len(gene_to_sources.get(gene, set()))),
                str(len(pmids)),
                pmids
            ]
        )

    if not rows:
        raise RuntimeError("Intersection is empty; cannot write summary table panel.")

    rows = rows[: int(max_rows)]
    # Create DataFrame but keep pmids as list for now
    summary_data = []
    for r in rows:
        summary_data.append({"gene": r[0], "n_mediators": int(r[1]), "n_pmids": int(r[2]), "pmids_list": r[3]})
        
    summary = pd.DataFrame(summary_data)

    genes = summary["gene"].tolist()[::-1]
    mediators = summary["n_mediators"].tolist()[::-1]
    n_pmids = summary["n_pmids"].tolist()[::-1]
    y = np.arange(len(genes))
    bar_h = 0.34

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    ax.barh(y - bar_h / 2, mediators, height=bar_h, color="#4C72B0", label="#Mediators")
    ax.barh(y + bar_h / 2, n_pmids, height=bar_h, color="#C44E52", label="#CTD PMIDs")
    
    for i, (m, p) in enumerate(zip(mediators, n_pmids)):
        ax.text(m + 0.1, i - bar_h / 2, str(m), va="center", ha="left", fontsize=9, color="#2c2c2c")
        ax.text(p + 0.1, i + bar_h / 2, str(p), va="center", ha="left", fontsize=9, color="#2c2c2c")
        
        # Add truncated PMID list if space permits or just keep it simple with counts? 
        # The prompt mentioned "PMID=0" issue and "text overflow".
        # Let's display a few PMIDs if count is small, otherwise "..."
        # Actually, the original code didn't print PMIDs, it printed counts.
        # "Fig2D.png the TNF CTD PMID=0 is real..." 
        # Ah, the user said "CTD curated table TNF pubmed_ids is empty". 
        # The chart plots counts. The TABLE refers to the CSV probably?
        # Re-reading: "Fig2D.png... TNF CTD PMID=0... not bug".
        # The original code only plotted bars of COUNTS.
        # So "text overflow" likely referred to my audit finding "If PubMed ID list over long..." 
        # Wait, I (the assistant) wrote the audit. I said "If PubMed ID list over long... text will overflow".
        # But the *original* code `save_intersection_table_pdf` only plotted BARS of counts. 
        # It did NOT plot the actual IDs.
        # Let's double check `save_intersection_table_pdf` in `old_string`.
        # Yes, it plots bars. `ax.barh`.
        # So where did I see text overflow?
        # Maybe I hallucinated that it prints text? 
        # Or maybe I meant the *Table* in the manuscript?
        # No, "Fig 2D (Intersection Table) ... ax.text ... If PubMed ID list...".
        # Use `read_file` earlier showed `save_intersection_table_pdf` creates a BAR CHART, not a text table.
        # It is called "Intersection Summary Table Panel" but it makes a BAR CHART.
        # Okay, I will stick to the bar chart but ensure it looks good (PNG, fonts).
        # I will also add a note about TNF PMIDs being 0 if needed, or just let it be 0.
        pass

    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlabel("Count", fontsize=10)
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right", fontsize=9)

    fig.tight_layout(pad=1.0)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_edge_list_tsv(
    output_path: str,
    *,
    disease_node: str,
    migraine_rows: list[CtdGeneDiseaseRow],
    intersection_genes: set[str],
    gene_to_sources: dict[str, set[str]],
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    by_gene: dict[str, list[CtdGeneDiseaseRow]] = defaultdict(list)
    for row in migraine_rows:
        by_gene[row.gene_symbol].append(row)

    with open(output_path, "w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(
            file_handle,
            fieldnames=["source", "target", "interaction", "evidence", "source_db"],
            delimiter="\t",
        )
        writer.writeheader()

        for gene in sorted(intersection_genes):
            ctd_evidence = sorted({row.direct_evidence for row in by_gene.get(gene, []) if row.direct_evidence})
            ctd_pmids = sorted(
                {pmid for row in by_gene.get(gene, []) for pmid in (row.pubmed_ids.split("|") if row.pubmed_ids else [])}
            )
            writer.writerow(
                {
                    "source": gene,
                    "target": disease_node,
                    "interaction": "ctd_gene_disease",
                    "evidence": f"direct_evidence={ '|'.join(ctd_evidence) };pmids={ '|'.join(ctd_pmids) }",
                    "source_db": "CTD_curated_genes_diseases",
                }
            )
            for molecule in sorted(gene_to_sources.get(gene, set())):
                writer.writerow(
                    {
                        "source": molecule,
                        "target": gene,
                        "interaction": "literature_mediator_gene_mapping",
                        "evidence": "mapping=curated_dictionary",
                        "source_db": "Task1.1",
                    }
                )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Collect migraine target genes (CTD curated) and intersect with acupuncture targets derived from Task 1.1 mediators.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """
            Outputs:
              - data/processed/migraine_targets_ctd.csv
              - data/processed/acupuncture_targets_from_literature.csv
              - data/processed/targets_intersection.csv
              - data/processed/task1.2_network_edges.tsv (Cytoscape-compatible edge list)
              - figures/raw_plots/Fig2A_venn.pdf
              - figures/raw_plots/Fig2B_network.pdf
              - figures/raw_plots/Fig2C_ctd_evidence_bar.pdf
              - figures/raw_plots/Fig2D_intersection_summary_table.pdf
            """
        ).strip(),
    )
    parser.add_argument(
        "--literature-csv",
        default="data/processed/literature_compounds.csv",
        help="Task 1.1 output CSV (default: data/processed/literature_compounds.csv)",
    )
    parser.add_argument(
        "--ctd-curated",
        default="data/raw/ctd/CTD_curated_genes_diseases.tsv.gz",
        help="CTD curated gene-disease report path (default: data/raw/ctd/CTD_curated_genes_diseases.tsv.gz)",
    )
    parser.add_argument("--ctd-url", default=DEFAULT_CTD_CURATED_URL, help="CTD curated report download URL.")
    parser.add_argument("--download-ctd", action="store_true", help="Download CTD curated report if missing.")
    parser.add_argument(
        "--disease-regex",
        default="migraine",
        help='Case-insensitive regex applied to CTD DiseaseName (default: "migraine").',
    )
    parser.add_argument(
        "--disease-node",
        default="Migraine",
        help='Disease node label used in outputs (default: "Migraine").',
    )
    parser.add_argument(
        "--out-migraine-targets",
        default="data/processed/migraine_targets_ctd.csv",
        help="Output CSV for CTD migraine targets.",
    )
    parser.add_argument(
        "--out-acu-targets",
        default="data/processed/acupuncture_targets_from_literature.csv",
        help="Output CSV for acupuncture targets derived from Task 1.1.",
    )
    parser.add_argument(
        "--out-intersection",
        default="data/processed/targets_intersection.csv",
        help="Output CSV for the intersection gene set.",
    )
    parser.add_argument(
        "--out-edges",
        default="data/processed/task1.2_network_edges.tsv",
        help="Output TSV edge list (Cytoscape compatible).",
    )
    parser.add_argument(
        "--fig-venn",
        default="figures/raw_plots/Fig2A.png",
        help="Output Venn figure (PNG/PDF).",
    )
    parser.add_argument(
        "--fig-network",
        default="figures/raw_plots/Fig2B.png",
        help="Output network figure (PNG/PDF).",
    )
    parser.add_argument(
        "--fig-evidence-bar",
        default="figures/raw_plots/Fig2C.png",
        help="Output CTD evidence bar plot (PNG/PDF).",
    )
    parser.add_argument(
        "--fig-intersection-table",
        default="figures/raw_plots/Fig2D.png",
        help="Output intersection summary table panel (PNG/PDF).",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    configure_logging(args.verbose)

    if not os.path.exists(args.literature_csv):
        LOGGER.error("Missing literature CSV: %s", args.literature_csv)
        return 2

    if not os.path.exists(args.ctd_curated):
        if not args.download_ctd:
            LOGGER.error("Missing CTD curated report: %s (use --download-ctd)", args.ctd_curated)
            return 2
        download_file(args.ctd_url, args.ctd_curated)

    disease_pattern = re.compile(args.disease_regex, re.I)
    migraine_rows = [row for row in read_ctd_curated(args.ctd_curated) if disease_pattern.search(row.disease_name)]
    if not migraine_rows:
        LOGGER.error("No CTD rows matched disease regex: %s", args.disease_regex)
        return 1

    write_migraine_targets_csv(args.out_migraine_targets, migraine_rows)
    LOGGER.info("CTD migraine rows: %s", len(migraine_rows))

    molecules_present = load_molecules_from_literature(args.literature_csv)
    gene_to_sources = build_acupuncture_target_genes(molecules_present)
    write_acupuncture_targets_csv(args.out_acu_targets, gene_to_sources)

    migraine_genes = {row.gene_symbol for row in migraine_rows if row.gene_symbol}
    acupuncture_genes = set(gene_to_sources.keys())
    intersection = migraine_genes & acupuncture_genes
    LOGGER.info("Migraine genes: %s | Acupuncture genes: %s | Intersection: %s", len(migraine_genes), len(acupuncture_genes), len(intersection))

    ctd_by_gene: dict[str, list[CtdGeneDiseaseRow]] = defaultdict(list)
    for row in migraine_rows:
        ctd_by_gene[row.gene_symbol].append(row)

    write_intersection_csv(args.out_intersection, intersection, gene_to_sources, ctd_by_gene)
    write_edge_list_tsv(
        args.out_edges,
        disease_node=args.disease_node,
        migraine_rows=migraine_rows,
        intersection_genes=intersection,
        gene_to_sources=gene_to_sources,
    )

    retrieved_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    save_venn_plot(
        args.fig_venn,
        a=len(acupuncture_genes),
        b=len(migraine_genes),
        ab=len(intersection),
        label_a="Acupuncture",
        label_b="Migraine",
    )
    save_network_plot(
        args.fig_network,
        disease_node=args.disease_node,
        intersection_genes=intersection,
        gene_to_sources=gene_to_sources,
    )
    save_ctd_evidence_bar_plot(
        args.fig_evidence_bar,
        migraine_rows=migraine_rows,
        intersection_genes=intersection,
        top_n=15,
    )
    save_intersection_table_plot(
        args.fig_intersection_table,
        migraine_rows=migraine_rows,
        intersection_genes=intersection,
        gene_to_sources=gene_to_sources,
        max_rows=12,
    )
    LOGGER.info(
        "Generated figures (retrieved_at=%s): %s , %s , %s , %s",
        retrieved_at,
        args.fig_venn,
        args.fig_network,
        args.fig_evidence_bar,
        args.fig_intersection_table,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
