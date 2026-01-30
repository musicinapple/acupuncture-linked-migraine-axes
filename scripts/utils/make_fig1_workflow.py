#!/usr/bin/env python3
"""
Generate Figure 1 workflow schematic (matplotlib).

Outputs are written to:
  - figures/final_figures/Fig1_Study_Workflow.(png|jpeg|pdf|svg)
"""

from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"]

BLUE = "#2C5F8A"
ORANGE = "#D97B29"
GREEN = "#4E8F58"
TEXT = "#1A1A1A"


def _wrap(text: str, width: int) -> str:
    lines: list[str] = []
    for part in text.split("\n"):
        part = part.strip()
        if not part:
            lines.append("")
            continue
        lines.append(fill(part, width=width, break_long_words=False, break_on_hyphens=False))
    return "\n".join(lines)


def header(ax, y_top: float, label: str, color: str) -> float:
    h = 0.055
    rect = FancyBboxPatch(
        (0.05, y_top - h),
        0.90,
        h,
        boxstyle="round,pad=0.005,rounding_size=0.01",
        facecolor=color,
        edgecolor=color,
        linewidth=0.0,
        zorder=2,
    )
    ax.add_patch(rect)
    ax.text(0.5, y_top - h / 2 + 0.003, label, ha="center", va="center", color="white", fontsize=11, fontweight="bold", zorder=3)
    return y_top - h - 0.02


def box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str | None,
    *,
    edge: str,
    face: str = "white",
    title_fs: float = 9.6,
    body_fs: float = 8.4,
    wrap_title: int | None = None,
    wrap_body: int | None = None,
    body_offset: float = 0.060,
    zorder: int = 10,
    title_align: str = "center",
    title_x_offset: float = 0.0,
) -> None:
    rect = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.01,rounding_size=0.015",
        facecolor=face,
        edgecolor=edge,
        linewidth=1.2,
        zorder=zorder,
    )
    ax.add_patch(rect)

    title_text = _wrap(title, wrap_title) if wrap_title else title

    if title_align == "center":
        tx = x + w / 2
        ha = "center"
    else:
        tx = x + 0.02 + title_x_offset
        ha = "left"

    ax.text(tx, y + h - 0.018, title_text, ha=ha, va="top", fontsize=title_fs, fontweight="bold", color=TEXT, zorder=zorder + 1)

    if body:
        body_text = _wrap(body, wrap_body) if wrap_body else body
        ax.text(x + 0.02, y + h - body_offset, body_text, ha="left", va="top", fontsize=body_fs, color=TEXT, linespacing=1.4, zorder=zorder + 1)


def arrow(ax, x1: float, y1: float, x2: float, y2: float, *, color: str) -> None:
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="->", lw=1.5, color=color, shrinkA=0, shrinkB=0, mutation_scale=15),
        zorder=5,
    )


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    out_dir = root / "figures" / "final_figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_png = out_dir / "Fig1_Study_Workflow.png"
    out_jpg = out_dir / "Fig1_Study_Workflow.jpeg"
    out_pdf = out_dir / "Fig1_Study_Workflow.pdf"
    out_svg = out_dir / "Fig1_Study_Workflow.svg"

    fig = plt.figure(figsize=(7.5, 10.5), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()

    # STAGE 1
    y1 = header(ax, 0.99, "STAGE 1 — KNOWLEDGE & DATA HARVESTING", BLUE)
    box(
        ax,
        0.30,
        y1 - 0.07,
        0.40,
        0.05,
        "Knowledge-Driven Target Discovery",
        None,
        edge=BLUE,
        face="#EBF2F8",
        wrap_title=40,
        title_fs=10,
    )

    pub_x, pub_y = 0.05, 0.67
    pub_w, pub_h = 0.40, 0.16
    box(
        ax,
        pub_x,
        pub_y,
        pub_w,
        pub_h,
        "PubMed mining (endogenous mediators)\n+ CTD (migraine targets)",
        "Intersection: 4 shared targets\nCALCA, CALCB, TAC1, TNF",
        edge=BLUE,
        wrap_body=35,
        title_fs=9.5,
        body_fs=9.0,
        body_offset=0.07,
    )

    stack_x = 0.50
    stack_w = 0.45
    card_h = 0.085
    gap = 0.010
    y_start = pub_y + pub_h - card_h

    box(
        ax,
        stack_x,
        y_start - 2 * (card_h + gap),
        stack_w,
        card_h + 0.01,
        "C. GSE198274 — Serum exosomal miRNA",
        "N=30 (10 Ctrl, 10 Pre, 10 Post)\nPaired cases",
        edge=BLUE,
        title_fs=9,
        body_fs=8.5,
        body_offset=0.05,
        title_align="left",
    )
    box(
        ax,
        stack_x,
        y_start - 1 * (card_h + gap),
        stack_w,
        card_h,
        "B. PRJEB40032 — PBMC RNA-seq",
        "N=42 (32 subj); Interictal/Ictal vs Basal",
        edge=BLUE,
        title_fs=9,
        body_fs=8.5,
        body_offset=0.05,
        title_align="left",
    )
    box(
        ax,
        stack_x,
        y_start,
        stack_w,
        card_h,
        "A. E-MTAB-13397 — Whole blood RNA-seq",
        "N=104 (52 subj); Interictal vs Control\nRepeated measures",
        edge=BLUE,
        title_fs=9,
        body_fs=8.5,
        body_offset=0.05,
        title_align="left",
    )

    arrow(ax, 0.50, y1 - 0.07, 0.25, pub_y + pub_h + 0.005, color=BLUE)
    arrow(ax, 0.50, y1 - 0.07, 0.725, y_start + card_h + 0.005, color=BLUE)
    arrow(ax, pub_x + pub_w, pub_y + pub_h / 2, stack_x - 0.01, y_start - card_h / 2, color=BLUE)

    # STAGE 2
    header(ax, 0.55, "STAGE 2 — ANALYSIS & MULTI-OMICS INTEGRATION", ORANGE)
    s2_y = 0.38
    s2_h = 0.12
    box(
        ax,
        0.05,
        s2_y,
        0.26,
        s2_h,
        "Subject-aware ML\nfeature selection",
        "LASSO / RF / SVM-RFE\nInternal AUC = 0.531\nExternal AUC = 0.483\nInv. AUC = 0.517",
        edge=ORANGE,
        wrap_body=22,
        title_fs=9,
        body_fs=8.5,
        body_offset=0.05,
    )
    box(
        ax,
        0.35,
        s2_y,
        0.30,
        s2_h,
        "miRNA–mRNA integration\n& pathway profiling",
        "miRNA baseline: 38\n(FDR<0.05)\nPost-tx delta: 0\n(FDR<0.05)",
        edge=ORANGE,
        wrap_body=24,
        title_fs=9,
        body_fs=8.5,
        body_offset=0.05,
    )
    box(
        ax,
        0.69,
        s2_y,
        0.26,
        s2_h,
        "Two-sample\nMendelian randomization",
        "Genetic triangulation\nNominal signals\nFDR>0.05\n(no robust drivers)",
        edge=ORANGE,
        wrap_body=22,
        title_fs=9,
        body_fs=8.5,
        body_offset=0.05,
    )

    arrow(ax, 0.31, s2_y + s2_h / 2, 0.35, s2_y + s2_h / 2, color=ORANGE)
    arrow(ax, 0.69, s2_y + s2_h / 2, 0.65, s2_y + s2_h / 2, color=ORANGE)

    # STAGE 3
    header(ax, 0.35, "STAGE 3 — RESULTS & INTERPRETATION", GREEN)
    s3_y = 0.18
    s3_h = 0.10
    box(
        ax,
        0.08,
        s3_y,
        0.38,
        s3_h,
        "32-gene migraine signature",
        "Intersection:\nML-selected ∩\nAcupuncture-linked",
        edge=GREEN,
        wrap_body=30,
        title_fs=9.5,
        body_fs=8.5,
        body_offset=0.055,
    )
    box(
        ax,
        0.54,
        s3_y,
        0.38,
        s3_h,
        "Molecular subtyping (k=3)",
        "General: C1=7, C2=15, C3=8\nAcupuncture: C1=9, C2=16, C3=5",
        edge=GREEN,
        wrap_body=35,
        title_fs=9.5,
        body_fs=8.5,
        body_offset=0.055,
    )
    arrow(ax, 0.46, s3_y + s3_h / 2, 0.54, s3_y + s3_h / 2, color=GREEN)

    # Conclusion
    c_x, c_y, c_w, c_h = 0.05, 0.05, 0.90, 0.10
    rect = FancyBboxPatch(
        (c_x, c_y),
        c_w,
        c_h,
        boxstyle="round,pad=0.01,rounding_size=0.015",
        facecolor="white",
        edgecolor=GREEN,
        linewidth=1.2,
        zorder=10,
    )
    ax.add_patch(rect)
    ax.text(
        c_x + c_w / 2,
        c_y + c_h - 0.018,
        "Conclusion (conservative)",
        ha="center",
        va="top",
        fontsize=9.6,
        fontweight="bold",
        color=TEXT,
        zorder=11,
    )
    ax.text(
        c_x + c_w / 2,
        c_y + 0.035,
        "Findings are consistent with state-dependent molecular heterogeneity and downstream response pathways\nrather than primary inherited drivers.",
        ha="center",
        va="bottom",
        fontsize=8.4,
        color=TEXT,
        zorder=11,
    )

    ax.text(
        0.50,
        0.020,
        "Legend: Blue=Input/Data  Orange=Analysis  Green=Interpretation/Output  Gray=Limitations",
        fontsize=7.6,
        color="#444444",
        ha="center",
        va="bottom",
    )

    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    fig.savefig(out_svg)
    plt.close(fig)

    from PIL import Image

    Image.open(out_png).convert("RGB").save(out_jpg, quality=95)
    print(f"Wrote: {out_png}")
    print(f"Wrote: {out_jpg}")
    print(f"Wrote: {out_pdf}")
    print(f"Wrote: {out_svg}")


if __name__ == "__main__":
    main()

